# Reindeer case study - SSF vs Langevin (Euler & brownian bridge)

# prep workspace ---------------------------------------------------------- ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, tidyr, mvnfast, parallel, terra, ggplot2, viridis,
         RColorBrewer, sf, oneimpact, survival, lubridate, mgcv)

# case study dir
cs_dir <- here("case_studies/reindeer/")

# import standardised covs (built by reindeer_spatial_covariates.R)
covs <- rast(here(cs_dir, "data/reindeer_covariates.tiff"))

# import tracks
data("reindeer")
reindeer <- sf::st_as_sf(reindeer, coords = c("x", "y"), crs = 25833)
# write.csv2(reindeer, here(cs_dir, "data", "reindeer_tracks.csv"))
# st_write(reindeer, here(cs_dir, "data", "reindeer_tracks.shp"))
# filter individuals to model 
unique(reindeer$animal_year_id)
subset_ids <- NULL

if (!is.null(subset_ids)) {
  reindeer <- reindeer |> 
    filter(animal_year_id %in% subset_ids)
}
# project tracks to covs crs 
tracks <- reindeer |>
  rename(time = t, ID = animal_year_id) |> 
  sf::st_transform(crs(covs)) |>
  terra::vect() |> 
  as.data.frame(geom = "XY") |> 
  dplyr::arrange(ID, time) |> 
  mutate(time = lubridate::round_date(time, unit = "hour"),
         dt = ifelse(ID == lead(ID),
                     lead(time) - time, NA),
         step = ifelse(ID == lead(ID),
                       sqrt((lead(x) - x)^2 + (lead(y) - y)^2), NA))

# generate smoothed covs for langevin gradient
covs_1km <- potential_GPU(covs, 
                          shape = "gaus",
                          alpha = dist_to_alpha(1, shape = "gaus"))
covs_2km <- potential_GPU(covs, 
                          shape = "gaus",
                          alpha = dist_to_alpha(2, shape = "gaus"))
covs_4km <- potential_GPU(covs, 
                          shape = "gaus",
                          alpha = dist_to_alpha(4, shape = "gaus"))
covs_5km <- potential_GPU(covs, 
                          shape = "gaus",
                          alpha = dist_to_alpha(5, shape = "gaus"))

# fit ssf based on top model in Panzacchi et al 2016. --------------------- ####
tracks2 <- read.csv(here(cs_dir, "data/dataframe.csv"))
tracks2$maxslope <- tracks2$maxslope^2

hab_terms <- c("maxslope", "maxTrailK", "sun", "prop2_lc", 
               "prop4_lc", "prop5_lc", "prop8_lc", "prop9_lc", "x_road", "maxRoK")
mov_terms <- c("realSL")
ssf_form  <- reformulate(c(hab_terms, mov_terms, "strata(strat)"),
                         response = "use")
ssf_fit   <- survival::clogit(ssf_form, data = tracks2)
coef(ssf_fit)

# comparable (with panzacchi et al 2016 covs)
tracks3 <- transmute(tracks2, strat = strat, use = use, step = realSL,
                     slope2 = maxslope, trail_density = maxTrailK,
                     road_density = maxRoK, solar_radiation = sun,
                     bog = prop2_lc, mountain_rock = prop4_lc,
                     mountain_veg = prop5_lc, 
                     lakes = prop8_lc, resevoirs = prop9_lc)
ssf_form3  <- reformulate(c(names(covs), "step", "strata(strat)"),
                          response = "use")
ssf_fit3   <- survival::clogit(ssf_form3, data = tracks3)
coef(ssf_fit3)

# conventional ssf (same data, only covs from end of step) ---------------- ####
# Build steps + random steps per animal-year, extract covariates at the step
# end point, and fit a conditional logistic model with a movement kernel.
trk <- amt::make_track(tracks, x, y, time, crs = terra::crs(covs), 
                       all_cols = TRUE) 
ssf <- trk |>
  tidyr::nest(.by = ID) |>
  dplyr::mutate(steps = purrr::map(data, function(d) {
    d |>
      amt::track_resample(rate = hours(3), tolerance = minutes(30)) |>
      amt::filter_min_n_burst(min_n = 3) |>
      amt::steps_by_burst(keep_cols = "start") |>
      amt::random_steps(n_control = 10)
  })) |>
  dplyr::select(-data) |>
  tidyr::unnest(steps) |>
  dplyr::mutate(step_id = paste(ID, burst_, step_id_, sep = "_"))

# covariates at the step end point (used = case_ TRUE, available = FALSE)
end_v <- terra::vect(as.matrix(ssf[, c("x2_", "y2_")]),
                     type = "points", crs = terra::crs(covs))
ecov  <- terra::extract(covs, end_v, ID = FALSE)
ssf   <- dplyr::bind_cols(ssf, ecov)

# movement-kernel terms (first step of each burst has NA turn angle -> that
# stratum is dropped by clogit; standard in iSSA)
ssf <- ssf |>
  dplyr::mutate(log_sl_ = log(sl_),
                cos_ta_ = cos(ta_))

mov_terms <- c("sl_", "log_sl_", "cos_ta_")
ssf_form  <- reformulate(c(names(covs), mov_terms, "strata(step_id)"),
                         response = "case_")
ssf_fit_conv <- survival::clogit(ssf_form, data = ssf)

# fit Euler Langevin ------------------------------------------------------ ####
covlist_km <- spatRast_to_covlist(covs)
locs       <- as.matrix(tracks[, c("x", "y")])
times      <- as.numeric(tracks$time) / 3600          # hours (within-ID lags used)
grad_array <- bilinearGradArray(locs, covlist_km)      # dim [n, 2, J]
# fit
euler      <- langevinUD(locs, times, ID = tracks$ID, grad_array = grad_array)
euler_beta <- setNames(euler$betaHat, names(covs))

# fit with smoothed covs
# 1km
grad_array <- bilinearGradArray(locs, spatRast_to_covlist(covs_1km))      # dim [n, 2, J]
# fit
euler1km      <- langevinUD(locs, times, ID = tracks$ID, grad_array = grad_array)
euler_beta1km <- setNames(euler1km$betaHat, names(covs))
# 2km
grad_array <- bilinearGradArray(locs, spatRast_to_covlist(covs_2km))      # dim [n, 2, J]
# fit
euler2kmkm       <- langevinUD(locs, times, ID = tracks$ID, grad_array = grad_array)
euler_beta2km <- setNames(euler2km$betaHat, names(covs))

# fit Langevin BBIS - single fit ------------------------------------------ ####
# reference BBIS 
ncores  <- 10
M       <- 20   # bridges for the single fit
dt_max  <- 3  # hours (reindeer fixes are ~3-hourly); single-fit value
bbis_M20_N3h  <- fit_langevin_bbis(tracks, covs,
                                   M = 2,
                                   dt_max = dt_max,
                                   dt_units = "hours",
                                   ncores = ncores,
                                   fixed_sampling = FALSE)

bbis_M20_N3h_1km  <- fit_langevin_bbis(tracks, 
                                      covs_1km,
                                      M = M,
                                      dt_max = dt_max,
                                      dt_units = "hours",
                                      ncores = ncores,
                                      fixed_sampling = FALSE)
bbis_M20_N3h_2km  <- fit_langevin_bbis(tracks, 
                                      covs_2km,
                                      M = M,
                                      dt_max = dt_max,
                                      dt_units = "hours",
                                      ncores = ncores,
                                      fixed_sampling = FALSE)
# 1 hr
dt_max  <- 1
bbis_M20_N1h  <- fit_langevin_bbis(tracks, covs,
                                   M = M,
                                   dt_max = dt_max,
                                   dt_units = "hours",
                                   ncores = ncores,
                                   fixed_sampling = FALSE)
# 30 min
dt_max  <- 30
bbis_M20_N30m  <- fit_langevin_bbis(tracks, covs,
                                    M = M,
                                    dt_max = dt_max,
                                    dt_units = "mins",
                                    ncores = ncores,
                                    fixed_sampling = FALSE)
# 15 min
dt_max  <- 15
bbis_M20_N15m  <- fit_langevin_bbis(tracks, covs,
                                    M = M,
                                    dt_max = dt_max,
                                    dt_units = "mins",
                                    ncores = ncores,
                                    fixed_sampling = FALSE)
# 10 min
dt_max  <- 10
bbis_M20_N10m  <- fit_langevin_bbis(tracks, covs,
                                    M = M,
                                    dt_max = dt_max,
                                    dt_units = "mins",
                                    ncores = ncores,
                                    fixed_sampling = FALSE)

# combine coefficient estimates from all models --------------------------- ####
      
model_names <- c("ssf_panzacchi", 
                 "ssf_conv",
                 "euler", "euler1km","euler2km", 
                 "bbis_3h", "bbis_3h_1km", "bbis_3h_2km",
                 "bbis_1h", "bbis_30m", "bbis_15m",
                 "bbis_10m")
coef_df <- data.frame(coef(ssf_fit3)[1:9],
                      coef(ssf_fit_conv)[1:9],
                      euler_beta, euler_beta1km, euler_beta2km,
                      bbis_M20_N3h$beta, bbis_M20_N3h_1km$beta, bbis_M20_N3h_2km$beta,
                      bbis_M20_N1h$beta, bbis_M20_N30m$beta, 
                      bbis_M20_N15m$beta, bbis_M20_N10m$beta) |> 
  setNames(model_names)
coef_df$coef <- rownames(coef_df)
rownames(coef_df) <- NULL
coef_df_long <- coef_df |> 
  pivot_longer(
    cols = all_of(model_names),
    names_to = "model",      # Name of the new column containing old headers
    values_to = "est"   # Name of the new column containing the values
  )
coef_df_long$model <- factor(coef_df_long$model, levels = model_names)

# save fit estimates
make_path(here("case_studies/reindeer/fitted_estimates"))
saveRDS(coef_df_long, here("case_studies/reindeer/fitted_estimates/reindeer_fits.rds"))

# summary plot of estimated coefficients ---------------------------------- ####
coef_df_long |> 
  ggplot(aes(x = coef, y = est, fill = model, col = model)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_col(width = 0.1, position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = NULL, y = "Estimated coefficient") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# predicted UD from estimated coefficients -------------------------------- ####
ud_ls <- list()
for (i in seq(model_names)) {
  lin_pred  <- sum(covs * coef_df[,model_names[i]])                          # linear predictor
  ud        <- exp(lin_pred)
  ud       <- ud / terra::global(ud, "sum", na.rm = TRUE)[[1]] # normalise to sum 1
  names(ud) <- "UD"
  ud_ls[[i]] <- ud
}

ud_stack <- rast(ud_ls) |> setNames(model_names)

# add marginal/empirical ud
obs_v  <- terra::vect(as.matrix(tracks[, c("x", "y")]),
                      type = "points", crs = terra::crs(covs))
cnt    <- terra::rasterize(obs_v, covs[[1]], fun = "length", background = 0)
cnt    <- terra::mask(cnt, covs[[1]])
bw_km  <- 0.25                                              # KDE bandwidth (km)
gk     <- terra::focalMat(cnt, d = bw_km, type = "Gauss")
obs_ud <- terra::focal(cnt, w = gk, na.rm = TRUE, na.policy = "omit") |> 
  terra::mask(covs[[1]])
obs_ud <- obs_ud / terra::global(obs_ud, "sum", na.rm = TRUE)[[1]]
ud_stack$obs_ud <- obs_ud

# save the predicted UD raster
writeRaster(ud_stack, here("case_studies/reindeer", "reindeer_predicted_UDs.tif"),
            overwrite = TRUE)

# Performance: Spearman corr. + Boyce indices (binned, continuous, smoothed) #
# Spearman rank correlation between predicted and observed UD across grid cells.
spear <- function(pred_r) {
  v <- terra::values(c(pred_r, obs_ud))
  v <- v[stats::complete.cases(v), , drop = FALSE]
  suppressWarnings(stats::cor(v[, 1], v[, 2], method = "spearman"))
}

# Original binned Boyce index (OBI; Boyce et al. 2002) 
# split the predicted UD into `nclass` fixed classes, compare the frequency of
# used points (E) to the frequency of available cells (P) per class, and take
# the Spearman correlation of the predicted-to-expected ratio F = E/P against
# class value.
boyce_index <- function(pred_r, used_pts, nclass = 20) {
  # UD value of used cells
  pau <- terra::extract(pred_r, used_pts, ID = FALSE)[, 1]
  pau <- pau[is.finite(pau)]
  # UD value of available cells
  av  <- terra::values(pred_r); av <- av[is.finite(av)]
  # define break values for UD distribution
  brks <- seq(min(av), max(av), length.out = nclass + 1)
  mids <- (brks[-1] + brks[-length(brks)]) / 2  # bin midpoints 
  # distribution of available UD (binned)
  Pi   <- as.numeric(table(cut(av,  brks, include.lowest = TRUE))) / length(av)
  # distribution of used UD (binned)
  Ei   <- as.numeric(table(cut(pau, brks, include.lowest = TRUE))) / length(pau)
  # ratio of used/available (binned)
  Fi   <- Ei / Pi
  # bins with data
  keep <- is.finite(Fi) & Pi > 0
  # return correlation midpoint
  suppressWarnings(stats::cor(mids[keep], Fi[keep], method = "spearman"))
}

# Continuous Boyce index (CBI; Hirzel et al. 2006)
# Moving-window version of the Boyce index. Instead of a small number of fixed
# bins, slide a window of width `window` (a fraction of the predicted-value
# range) across that range in `res` steps. For each window compute the P/E
# ratio F = (used pts in window / n_used) / (available cells in window /
# n_cells), then take the Spearman correlation of F against the window
# midpoints. This is the algorithm implemented by ecospat::ecospat.boyce()
# (defaults window.w = 10% of range, res = 100) and removes the sensitivity of
# the OBI to the number and starting location of bins.
boyce_continuous <- function(pred_r, used_pts, window = 0.1, res = 100) {
  # predicted UD at used (presence) points
  obs <- terra::extract(pred_r, used_pts, ID = FALSE)[, 1]
  obs <- obs[is.finite(obs)]
  # predicted UD across all available cells
  fit <- terra::values(pred_r); fit <- fit[is.finite(fit)]
  mini <- min(fit); maxi <- max(fit)
  window.w <- window * (maxi - mini)                 # moving-window width
  vec.mov  <- seq(mini, maxi - window.w, length.out = res + 1)  # window lows
  intervals <- cbind(vec.mov, vec.mov + window.w)
  # P/E ratio in each moving window
  Fi <- apply(intervals, 1, function(w) {
    Pi <- mean(obs >= w[1] & obs <= w[2])            # used freq. in window
    Ei <- mean(fit >= w[1] & fit <= w[2])            # available freq. in window
    Pi / Ei
  })
  mids <- vec.mov + window.w / 2                      # window midpoints
  keep <- is.finite(Fi)                              # drop empty windows (Ei=0)
  if (sum(keep) < 2) return(NA_real_)
  suppressWarnings(stats::cor(Fi[keep], mids[keep], method = "spearman"))
}

# calc spearman correlation and Boyce indices
spearman <- list()
for (i in seq(model_names)) {
  spearman[[i]] <- spear(ud_stack[[i]])
}

boyce <- list()          # original binned Boyce index (OBI)
for (i in seq(model_names)) {
  boyce[[i]] <- boyce_index(ud_stack[[i]], obs_v, nclass = 100)
}

boyce_cbi <- list()      # continuous Boyce index (Hirzel et al. 2006)
for (i in seq(model_names)) {
  boyce_cbi[[i]] <- boyce_continuous(ud_stack[[i]], obs_v,
                                     window = 0.025, res = 100)
}

## combine into df
perf_df <- data.frame(
  method    = model_names,
  spearman  = unlist(spearman),
  # pearson   = unlist(pearson),
  boyce_obi = unlist(boyce),
  boyce_cbi = unlist(boyce_cbi)
) |>
  dplyr::mutate(method = factor(method,
                                levels = model_names))
## piv long
perf_long <- tidyr::pivot_longer(perf_df, c(spearman,
                                            # pearson,
                                            boyce_obi, boyce_cbi),
                                 names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric,
                                levels = c("spearman", "boyce_obi",
                                           "boyce_cbi")))
# plot
perf_plot <- ggplot(perf_long, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.3, size = 3) +
  facet_wrap(~ metric, scales = "free_y",
             labeller = as_labeller(
               c(spearman  = "Spearman corr. (predicted vs observed UD)",
                 # pearson   = "pearson corr. (predicted vs observed UD)",
                 boyce_obi = "Boyce index, binned (Boyce et al. 2002)",
                 boyce_cbi = "Continuous Boyce index (Hirzel et al. 2006)",
                 boyce_sbi = "Smoothing Boyce index, tp (Liu et al. 2024)"))) +
  labs(x = NULL, y = NULL,
       title = "UD-prediction performance by method") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
perf_plot

# Out-of-sample validation: leave-one-animal-out cross-validation ---------- ####
# refit each model on all but one animal-year and score its predicted UD on 
# the withheld animal's fixes. This is the k-fold Boyce cross-validation of Boyce et al. (2002),
# applied to a predicted UD raster so it is identical for SSFs and Langevin
# models. We evaluate the withheld *fixes* directly (used-vs-available), which
# needs no empirical/KDE reference surface.
#
# Metrics per fold:
#   cbi     - continuous Boyce index (Hirzel 2006) on withheld fixes  [discrim.]
#   obi     - binned Boyce index on withheld fixes                    [discrim.]
#   logdens - mean log predicted UD at withheld fixes                 [magnitude]
#             (the UD is a normalised density on a common grid, so this proper
#              score is comparable across model types; higher = better)
#
cv_M   <- 20  # bridges per BBIS refit (as single fit)
set.seed(1)

# build a normalised predicted UD from a habitat-beta vector -------------- ####
ud_from_beta <- function(beta) {
  b  <- setNames(as.numeric(beta), names(covs))[names(covs)]
  ud <- exp(sum(covs * b))
  ud / terra::global(ud, "sum", na.rm = TRUE)[[1]]
}

# model registry: each refits on a training track set `tr` and returns ---- ####
#     a habitat-beta vector named over names(covs). Comment out to skip.
cv_fitters <- list(
  ssf_conv = function(tr) {
    st <- amt::make_track(tr, x, y, time, crs = terra::crs(covs),
                          all_cols = TRUE) |>
      tidyr::nest(.by = ID) |>
      dplyr::mutate(steps = purrr::map(data, function(d) {
        d |>
          amt::track_resample(rate = hours(3), tolerance = minutes(30)) |>
          amt::filter_min_n_burst(min_n = 3) |>
          amt::steps_by_burst(keep_cols = "start") |>
          amt::random_steps(n_control = 10)
      })) |>
      dplyr::select(-data) |>
      tidyr::unnest(steps) |>
      dplyr::mutate(step_id = paste(ID, burst_, step_id_, sep = "_"))
    ev <- terra::vect(as.matrix(st[, c("x2_", "y2_")]),
                      type = "points", crs = terra::crs(covs))
    st <- dplyr::bind_cols(st, terra::extract(covs, ev, ID = FALSE))
    st <- dplyr::mutate(st, log_sl_ = log(sl_), cos_ta_ = cos(ta_))
    form <- reformulate(c(names(covs), c("sl_", "log_sl_", "cos_ta_"),
                          "strata(step_id)"), response = "case_")
    coef(survival::clogit(form, data = st))[names(covs)]
  },
  euler = function(tr) {
    locs <- as.matrix(tr[, c("x", "y")])
    ga   <- bilinearGradArray(locs, covlist_km)
    e    <- langevinUD(locs, as.numeric(tr$time) / 3600, ID = tr$ID,
                       grad_array = ga)
    setNames(e$betaHat, names(covs))
  },
  bbis_3h  = function(tr) fit_langevin_bbis(tr, covs, M = cv_M, dt_max = 3,
                  dt_units = "hours", ncores = ncores, fixed_sampling = FALSE)$beta,
  bbis_1h  = function(tr) fit_langevin_bbis(tr, covs, M = cv_M, dt_max = 1,
                  dt_units = "hours", ncores = ncores, fixed_sampling = FALSE)$beta,
  bbis_30m = function(tr) fit_langevin_bbis(tr, covs, M = cv_M, dt_max = 30,
                  dt_units = "mins",  ncores = ncores, fixed_sampling = FALSE)$beta,
  bbis_15m = function(tr) fit_langevin_bbis(tr, covs, M = cv_M, dt_max = 15,
                  dt_units = "mins",  ncores = ncores, fixed_sampling = FALSE)$beta,
  bbis_10m = function(tr) fit_langevin_bbis(tr, covs, M = cv_M, dt_max = 10,
                  dt_units = "mins",  ncores = ncores, fixed_sampling = FALSE)$beta,
  bbis_3h_1km = function(tr) fit_langevin_bbis(tr, covs_1km, M = cv_M, dt_max = 3,
                  dt_units = "hours", ncores = ncores, fixed_sampling = FALSE)$beta
)
cv_models <- names(cv_fitters)            # <- subset here for a quicker run

# score one predicted UD against a set of withheld fixes ------------------ ####
eval_ud <- function(ud, test_pts) {
  c(cbi     = boyce_continuous(ud, test_pts, window = 0.025, res = 100),
    obi     = boyce_index(ud, test_pts, nclass = 100),
    logdens = mean(log(terra::extract(ud, test_pts, ID = FALSE)[, 1]),
                   na.rm = TRUE))
}

# leave-one-animal-out loop ----------------------------------------------- ####
folds   <- sort(unique(tracks$ID))
cv_rows <- list()
for (f in folds) {
  tr      <- tracks[tracks$ID != f, ]                 # training animals
  test_pts <- terra::vect(as.matrix(tracks[tracks$ID == f, c("x", "y")]),
                          type = "points", crs = terra::crs(covs))
  for (m in cv_models) {
    message(sprintf("CV fold %s | model %s", f, m))
    beta <- tryCatch(cv_fitters[[m]](tr), error = function(e) {
      warning(sprintf("  %s failed on fold %s: %s", m, f, conditionMessage(e)))
      setNames(rep(NA_real_, length(names(covs))), names(covs))
    })
    ev <- eval_ud(ud_from_beta(beta), test_pts)
    cv_rows[[length(cv_rows) + 1L]] <- data.frame(
      fold = f, model = m, metric = names(ev), value = as.numeric(ev))
  }
}
cv_df <- dplyr::bind_rows(cv_rows) |>
  dplyr::mutate(model  = factor(model,  levels = model_names),
                metric = factor(metric, levels = c("cbi", "obi", "logdens")))

# save CV results
saveRDS(cv_df, here("case_studies/reindeer/fitted_estimates/reindeer_cv.rds"))

# summarise (mean +/- SE across folds) and plot --------------------------- ####
cv_summary <- cv_df |>
  dplyr::group_by(model, metric) |>
  dplyr::summarise(mean = mean(value, na.rm = TRUE),
                   se   = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value))),
                   .groups = "drop")

cv_plot <- ggplot(cv_summary, aes(x = model, y = mean, colour = model)) +
  geom_point(data = cv_df, aes(y = value), alpha = 0.25,
             position = position_jitter(width = 0.12), show.legend = FALSE) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se), show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y",
             labeller = as_labeller(
               c(cbi     = "CBI, withheld fixes (Hirzel 2006)",
                 obi     = "Binned Boyce, withheld fixes (Boyce 2002)",
                 logdens = "Mean log predicted UD at withheld fixes"))) +
  labs(x = NULL, y = NULL,
       title = "Leave-one-animal-out cross-validated UD prediction",
       subtitle = "points = per-animal folds; range = mean +/- SE") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
cv_plot


# # bar plot of the habitat-selection coefficients from the single fit:
# # covariates as categories on x, estimates as bar heights
# model <- model_names[1]
# 
# model_coef_df <- data.frame(
#   covariate = factor(names(covs), levels = names(covs)),
#   estimate  = coef_df[,model]
# )
# 
# coef_plot <- ggplot(model_coef_df, aes(x = covariate, y = estimate,
#                                  fill = estimate > 0)) +
#   geom_col(width = 0.7, show.legend = FALSE) +
#   geom_hline(yintercept = 0, linewidth = 0.3) +
#   scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02")) +
#   labs(x = NULL, y = "Estimated coefficient",
#        title = sprintf("Langevin BBIS coefficients (M = %s, dt_max = %s h)",
#                        M, dt_max)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1))
# 
# ggsave(here("case_studies/reindeer", "reindeer_coef_estimates.png"), coef_plot,
#        width = 7, height = 5)
# 
# 
# 
# 
# 
# # map of the predicted UD
# ud_plot <- ggplot() +
#   tidyterra::geom_spatraster(data = ud) +
#   scale_fill_viridis_c(option = "inferno", na.value = "transparent", name = "UD") +
#   coord_sf(expand = FALSE) +
#   labs(title = "Predicted utilisation distribution (single fit)",
#        x = NULL, y = NULL) +
#   theme_minimal()
# 
# ggsave(here("case_studies/reindeer", "reindeer_predicted_UD.png"), ud_plot,
#        width = 7, height = 6)

# #### response curves for polynomial covariates (single fit) ####
# # Marginal relative selection for each quadratic covariate: how exp(linear
# # predictor) varies as that covariate moves across its observed (standardised)
# # range while the others are held at their mean (0 on the z-score scale).
# # A negative squared coefficient gives an interior optimum at -b1 / (2 b2).
# beta_named <- setNames(out$beta, names(covs))
# 
# rng_of <- function(nm) c(terra::global(covs[[nm]], "min", na.rm = TRUE)[[1]],
#                          terra::global(covs[[nm]], "max", na.rm = TRUE)[[1]])
# 
# resp_df <- dplyr::bind_rows(lapply(poly_covs, function(nm) {
#   rng <- rng_of(nm)
#   xs  <- seq(rng[1], rng[2], length.out = 200)
#   b1  <- beta_named[[nm]]
#   b2  <- beta_named[[paste0(nm, "_sq")]]
#   lp  <- b1 * xs + b2 * xs^2
#   data.frame(covariate = nm, x = xs, rel_sel = exp(lp - max(lp)))
# }))
# 
# opt_df <- dplyr::bind_rows(lapply(poly_covs, function(nm) {
#   rng   <- rng_of(nm)
#   b1    <- beta_named[[nm]]; b2 <- beta_named[[paste0(nm, "_sq")]]
#   xstar <- -b1 / (2 * b2)
#   data.frame(covariate = nm, xstar = xstar,
#              has_opt = is.finite(xstar) && b2 < 0 &&
#                xstar >= rng[1] && xstar <= rng[2])
# }))
# 
# resp_plot <- ggplot(resp_df, aes(x = x, y = rel_sel)) +
#   geom_line(linewidth = 0.8, colour = "#1b9e77") +
#   geom_vline(data = subset(opt_df, has_opt), aes(xintercept = xstar),
#              linetype = 2, colour = "grey40") +
#   facet_wrap(~ covariate, scales = "free", ncol = 3) +
#   labs(x = "Standardised covariate value (z-score)",
#        y = "Relative selection  exp(linear predictor), peak = 1",
#        title = "Langevin BBIS response curves (quadratic terms)") +
#   theme_bw()
# 
# ggsave(here("case_studies/reindeer", "reindeer_response_curves.png"),
#        resp_plot, width = 10, height = 4)
# 
# #### fit Langevin BBIS - dt_max & M refits ####
# # define fitting criteria
# ncores  <- 10
# Ms      <- c(25, 50, 100)
# deltas  <- exp(seq(log(0.5), log(24), length.out = 30))  # delta_max grid (hours)
# nrefits <- 10
# 
# # number of pars
# npar <- nlyr(covs) + 1
# # add 1 column to track dt_max
# params <- matrix(NA, ncol = npar + 1, nrow = nrefits * length(deltas))
# 
# for (M in Ms) {                       # for each number of bridges
#   for (k in seq_along(deltas)) {      # for each delta_max
#     for (i in 1:nrefits) {
#       delta_max <- deltas[k]
#       
#       print(sprintf("Fitting M = %s, delta_max = %.4f, refit = %s",
#                     M, delta_max, i))
#       # fit
#       out <- fit_langevin_bbis(tracks, covs,
#                                M = M,
#                                dt_max = delta_max,
#                                dt_units = "hours",
#                                ncores = ncores,
#                                fixed_sampling = FALSE)
#       
#       # store outputs (par + delta_max)
#       params[(k - 1L) * nrefits + i, ] <- c(out$par, delta_max)
#     }
#   }
#   # convert to data.frame, name it dfM (df25 / df50 / df100), and save
#   df <- as.data.frame(params) |>
#     setNames(c(paste0("beta", seq_len(npar - 1L)), "sigma", "delta_max"))
#   assign(sprintf("df%s", M), df)
#   save(list = sprintf("df%s", M),
#        file = sprintf("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=%s.Rda",
#                       M))
# }
# 
# #### generate summary plots ####
# # import estimates
# load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=25.Rda")
# load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=50.Rda")
# load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=100.Rda")
# 
# # parameter columns (all betas + sigma), labelled for facet strips
# beta_names <- paste0("beta", seq_len(nlyr(covs)))
# cov_labels <- names(covs)
# par_levels <- c(beta_names, "sigma")
# par_labels <- c(setNames(sprintf("beta[%d]~(%s)", seq_along(cov_labels), cov_labels),
#                          beta_names),
#                 sigma = "sigma")
# 
# # combine all data
# df_all <- bind_rows(mutate(df25,  M = "M=25"),
#                     mutate(df50,  M = "M=50"),
#                     mutate(df100, M = "M=100")) |>
#   pivot_longer(cols = all_of(par_levels),
#                names_to = "par", values_to = "mu") |>
#   mutate(par = factor(par, levels = par_levels))
# 
# # summarise estimates (median, sd, & confidence intervals)
# z <- qnorm(0.975)
# sum_all <- df_all |>
#   dplyr::group_by(par, delta_max, M) |>
#   dplyr::summarise(sd = sd(mu),
#                    mu = median(mu),
#                    .groups = "drop") |>
#   dplyr::mutate(lo = mu - z * sd,
#                 hi = mu + z * sd)
# 
# # generate plot
# plot <- ggplot(sum_all, aes(x = delta_max, y = mu,
#                             color = factor(M, levels = c("M=25", "M=50", "M=100")))) +
#   # BBIS estimates
#   geom_point(data = df_all, alpha = 0.15, stroke = NA) +
#   geom_line(aes(linetype = factor(M, levels = c("M=25", "M=50", "M=100"))),
#             linewidth = 0.7) +
#   # design
#   facet_wrap(~ par, scales = "free",
#              labeller = labeller(par = as_labeller(par_labels, label_parsed))) +
#   scale_x_log10() +
#   scale_color_brewer(palette = "Dark2") +
#   labs(x = expression(Delta[max]), y = expression(Estimate),
#        color = NULL, linetype = NULL) +
#   theme_bw()
# 
# # save plot
# ggsave(here("case_studies/reindeer", "reindeer_par_est_dmax_M.png"), plot,
#        width = 9, height = 6)
