# Reindeer case study - thinning experiment: summarise 
# load fits from reindeer_fit.R, score UDs, average over clock phases, summarise, write csvs
# Betas and scores are first averaged over the clock phases in each leave-one-out 
# CV fold, then summarised over the folds:
#   * betas: jackknife SE (folds share all but one animal, so their sd
#     understates the uncertainty; see jack_se())
#   * scores: SE over withheld animals (each fold scores a different animal),
#     and paired differences from BBIS (removes how predictable each animal is)

# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))                             # ud_from_beta(), eval_ud(), jack_se(), ...
source(here("case_studies/reindeer/R/reindeer_data.R"))  # load_reindeer_tracks()
source(here("case_studies/reindeer/R/reindeer_covariate_smoothing.R"))
# Load necessary packages
load_lib(terra, dplyr, tidyr, lubridate,  # data
         ggplot2)                         # visualisation

# settings ---------------------------------------------------------------- ####
exp_name       <- "thinning_fixed_hour"   # as in reindeer_fit.R
fits_dir       <- here("case_studies/reindeer/fitted_estimates/experiments", exp_name)
score_dir      <- file.path(fits_dir, "scores")  # UD scores cached here between runs
results_dir    <- here("case_studies/reindeer/results")
metrics        <- c("logdens", "kstar", "logdens_kstar",  # magnitude-sensitive (main)
                    "auc", "cbi", "tss")#,                  # rank-only (supplementary)
                    # "UDcor")                              # ud overlap
ref_interval   <- 3                       # native GPS interval: reference for drift
anchor_legacy  <- 12                      # as in reindeer_fit.R (fits without anchor_h)
converged_only <- TRUE                    # drop non-converged fits from summaries
make_path(results_dir)
make_path(score_dir)

# helper functions
## generate output file paths
out <- function(suffix) file.path(results_dir, paste0(exp_name, "_", suffix))

## calculate overlap between two UD rasters 
ud_overlap <- function(r1, r2) {
  v  <- na.omit(cbind(values(r1, mat = FALSE), values(r2, mat = FALSE)))
  a  <- pmax(v[, 1], 0); b <- pmax(v[, 2], 0)      # clamp smoothing negatives
  p  <- a / sum(a);      q <- b / sum(b)           # normalise to prob mass
  c(spearman = cor(a, b, method = "spearman"),
    BA       = sum(sqrt(p * q)),                   # Bhattacharyya affinity
    D        = 1 - 0.5 * sum(abs(p - q)))          # Schoener's D  ( = VI )
}

## restrict a table to converged fits
use <- function(df) if (converged_only) df[df$converged, , drop = FALSE] else df

# load data --------------------------------------------------------------- ####
covs_raw    <- rast(here("case_studies/reindeer/data/reindeer_covariates.tiff"))
covs_smooth <- rast(here("case_studies/reindeer/data/smoothed/reindeer_covariates_smooth_0.5km.tif")) |>
  setNames(names(covs_raw))
cov_fields  <- list(`0` = covs_raw, `0.5` = covs_smooth)  # as in reindeer_fit.R
tracks <- load_reindeer_tracks(covs_raw)

# load checkpoints -------------------------------------------------------- ####
# one record per fit: the data set (method, interval, anchor_h, fold, ...) + the fit
recs <- collect_checkpoints(fits_dir)
if (!length(recs)) stop("no checkpoints in ", fits_dir)
message(sprintf("%s: %d fits read", exp_name, length(recs)))
# fits run before anchors were added carry no anchor_h: they used 12:00 UTC
recs <- lapply(recs, function(r) {
  if (is.null(r$anchor_h)) r$anchor_h <- anchor_legacy %% r$interval
  r
})

# per-fit table: data set, fit scalars and bridge-ESS summary
fits <- bind_rows(lapply(recs, function(r) {
  data.frame(r[c("method", "interval", "anchor_h", "fold", "smooth_km", "M", "label", "n_train",
                 "sigma", "loglik", "convergence", "converged", "n_attempts", "time_secs")],
             as.list(r$ess))
}))
# save
write.csv(fits, out("fits.csv"), row.names = FALSE)

# estimated coefficients: one row per fit x covariate
betas <- bind_rows(lapply(recs, function(r) {
  data.frame(method = r$method, interval = r$interval, anchor_h = r$anchor_h, fold = r$fold,
             converged = r$converged, coef = names(r$beta), beta = unname(r$beta))
}))
# save
write.csv(betas, out("betas.csv"), row.names = FALSE)

# score UDs --------------------------------------------------------------- ####
# stationary UD exp(beta . c) on the method's own covariate field, scored
# against the withheld animal's full 3 h fixes (out-of-sample)
message("scoring UDs (cached in ", score_dir, ") ...")
scores <- bind_rows(lapply(recs, function(rec_i) {
  file <- file.path(score_dir, paste0(rec_i$label, ".rds"))
  score_i <- if (file.exists(file)) readRDS(file)
  # score if not cached, or cached before a metric was added
  if (!all(metrics %in% names(score_i))) {
    score_i <- setNames(rep(NA_real_, length(metrics)), metrics)
    if (all(is.finite(rec_i$beta))) {
      ud   <- ud_from_beta(rec_i$beta, cov_fields[[as.character(rec_i$smooth_km)]])
      test <- tracks[tracks$ID == rec_i$fold, ]
      score_i    <- eval_ud(ud, vect(as.matrix(test[, c("x", "y")]), crs = crs(ud)),
                      metrics = metrics)
    }
    saveRDS(score_i, file)
  }
  score_i <- score_i[metrics]
  data.frame(method = rec_i$method, interval = rec_i$interval, 
             anchor_h = rec_i$anchor_h, fold = rec_i$fold,
             converged = rec_i$converged, metric = names(score_i), 
             value = unname(score_i))
}))

# UD correlation score
## individual KDE
KDE <- list()
for (i in unique(tracks$ID)) {
  KDE[[i]] <-
    # pull track & convert to spatial
    filter(tracks, ID == i) %>% 
    dplyr::select(x, y) %>% 
    as.matrix() %>% 
    vect(crs = crs(cov_fields[[1]], T)) %>% 
    rasterize(covs[[1]], fun = "length", background = 0) %>% 
    potential_GPU(shape = "gaus", 
                  alpha = dist_to_alpha(1, "gaus"))
}
KDE <- rast(KDE) %>% 
  setNames(as.character(unique(tracks$ID)))

# calculate UD correlation score
UDCor_scores <- bind_rows(lapply(recs, function(rec_i) {
    if (all(is.finite(rec_i$beta))) {
      ud <- ud_from_beta(rec_i$beta, cov_fields[[as.character(rec_i$smooth_km)]])
      score_i <- ud_overlap(ud,
                        KDE[[as.character(rec_i$fold)]])
      # score_i <- ud_overlap(values(ud, ID = F),
      #                   values(KDE[[as.character(rec_i$fold)]], ID = F))[[1]]
    }
  data.frame(method = rec_i$method, interval = rec_i$interval, 
             anchor_h = rec_i$anchor_h, fold = rec_i$fold,
             converged = rec_i$converged,
             metric = sprintf("UDCor_%s", names(score_i)), 
             value = unname(score_i))
}))

scores <- bind_rows(scores, UDCor_scores)

# save
write.csv(scores, out("scores.csv"), row.names = FALSE)


# summarise betas --------------------------------------------------------- ####
# phase-averaged estimate per fold: mean over the clock phases (anchors) of a schedule
beta_fold <- use(betas) %>%
  group_by(method, interval, fold, coef) %>%
  summarise(beta = mean(beta), n_anchor = n(), .groups = "drop")

# mean +/- jackknife SE over leave-one-out folds, per method x interval
beta_sm <- beta_fold %>%
  group_by(method, interval, coef) %>%
  summarise(mean = mean(beta), sd = sd(beta), n = n(), se_jack = jack_se(beta),
            n_anchor = mean(n_anchor), .groups = "drop")
# save
write.csv(beta_sm, out("beta_summary.csv"), row.names = FALSE)

# time-of-day effect: betas per clock phase (UTC hour kept), before averaging
beta_phase <- use(betas) %>%
  group_by(method, interval, anchor_h, coef) %>%
  summarise(mean = mean(beta), n = n(), se_jack = jack_se(beta), .groups = "drop")
# save
write.csv(beta_phase, out("beta_phase.csv"), row.names = FALSE)

# drift from the native 3 h, paired by fold (same withheld animal):
#   d     = beta - beta_3h       signed change
#   d_abs = |beta| - |beta_3h|   > 0: further from 0 (stronger selection/avoidance);
#                                misleading where the mean sign switches (sign_switch)
beta_drift <- beta_fold %>%
  group_by(method, fold, coef) %>%
  mutate(b_ref = beta[interval == ref_interval][1]) %>%
  ungroup() %>%
  filter(interval != ref_interval) %>%
  mutate(d = beta - b_ref, d_abs = abs(beta) - abs(b_ref)) %>%
  group_by(method, interval, coef) %>%
  summarise(sign_switch = sign(mean(beta)) != sign(mean(b_ref, na.rm = TRUE)),
            across(c(d, d_abs), list(mean = ~ mean(.x, na.rm = TRUE), se = jack_se)),
            .groups = "drop") %>%
  mutate(z = d_mean / d_se, z_abs = d_abs_mean / d_abs_se)
# save
write.csv(beta_drift, out("beta_drift.csv"), row.names = FALSE)

# whole-vector shift from the native 3 h betas: L2 distance and angle.
# (mean |beta| is not used: it cancels betas moving in opposite directions)
beta_ref <- beta_sm %>%
  filter(interval == ref_interval) %>%
  select(method, coef, ref = mean)
beta_shift <- beta_sm %>%
  left_join(beta_ref, by = c("method", "coef")) %>%
  group_by(method, interval) %>%
  summarise(L2 = sqrt(sum((mean - ref)^2)),
            angle_deg = acos(pmin(1, sum(mean * ref) / sqrt(sum(mean^2) * sum(ref^2)))) * 180 / pi,
            .groups = "drop")
# save
write.csv(beta_shift, out("beta_shift.csv"), row.names = FALSE)

# summarise scores -------------------------------------------------------- ####
# per fold (= withheld animal): mean over the clock phases
score_fold <- scores %>%
  group_by(method, interval, fold, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(value = ifelse(is.nan(value), NA, value))

# mean +/- SE over withheld animals (each fold scores a different animal)
score_sm <- score_fold %>%
  group_by(method, interval, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            n = sum(is.finite(value)), .groups = "drop") %>%
  mutate(se = sd / sqrt(n))
# save
write.csv(score_sm, out("score_summary.csv"), row.names = FALSE)

# paired by withheld animal: each method minus BBIS (removes how predictable each
# animal is, which dominates the spread of the raw scores)
score_rel <- score_fold %>%
  group_by(interval, fold, metric) %>%
  mutate(value = value - value[method == "bbis"][1]) %>%
  group_by(method, interval, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            n = sum(is.finite(value)), .groups = "drop") %>%
  mutate(se = sd / sqrt(n), z = mean / se)
# save
write.csv(score_rel, out("score_vs_bbis.csv"), row.names = FALSE)

# fit diagnostics --------------------------------------------------------- ####
# per level, pooled over clock phases; beta_se = jackknife SE (mean over covariates)
beta_spread <- beta_sm %>%
  group_by(method, interval) %>%
  summarise(beta_se = mean(se_jack, na.rm = TRUE), .groups = "drop")
diag_sm <- fits %>%
  group_by(method, interval) %>%
  summarise(n_fits = n(),
            converged_rate = mean(converged),
            n_attempts = mean(n_attempts, na.rm = TRUE),
            n_train = mean(n_train),
            sigma = mean(sigma[converged], na.rm = TRUE),
            time_secs = mean(time_secs, na.rm = TRUE),
            ess_frac_mean = mean(ess_frac_mean, na.rm = TRUE),
            frac_ess_lt10 = mean(frac_ess_lt10, na.rm = TRUE),
            N_mean = mean(N_mean, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(beta_spread, by = c("method", "interval")) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA, .x)))
# save
write.csv(diag_sm, out("diagnostic_summary.csv"), row.names = FALSE)

# predict UD maps --------------------------------------------------------- ####
# predicted UD averaged over folds and clock phases, per method x interval,
# on a shared colour scale
keys <- distinct(use(fits), method, interval, smooth_km) %>% arrange(method, interval)
uds <- lapply(seq_len(nrow(keys)), function(i) {
  rs <- Filter(function(r) r$method == keys$method[i] && r$interval == keys$interval[i] &&
                 (!converged_only || isTRUE(r$converged)), recs)
  mean_ud(lapply(rs, `[[`, "beta"), cov_fields[[as.character(keys$smooth_km[i])]])
})
keep <- !vapply(uds, is.null, logical(1))
ud_stack <- rast(uds[keep])
names(ud_stack) <- sprintf("%s | %g h", keys$method[keep], keys$interval[keep])
writeRaster(ud_stack, out("mean_UDs.tif"), overwrite = TRUE)
