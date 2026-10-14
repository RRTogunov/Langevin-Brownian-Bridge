# Reindeer case study - thinning experiment to test beta consistency at different sampling resolutions
# Leave-one-out (one animal-year withheld) for cross validation (fit to all but one, predict on left out individual)
# To remove effect of time-of-day effect, training data are thinned (3, 6, 12 or 24 h)
# at every clock phase the 3 h GPS fix rate allows. i.e.,
# 3 h -> 1 phase, reference model
# 6 h -> 2 phases (00 & 03 h)
# 12 h -> 4 phases (00, 03, 06, 09 h)
# 24 h -> 8 phases (00, 03, 06, 09, 12, 15, 18, 21 h)
# fit iSSF, Euler, and BBIS to all thinned data sets 
# each fit is defined by fold (individual), sampling interval, and anchor hour (phase)
#
# Usage (from the repo root):
#   Rscript case_studies/reindeer/reindeer_fit.R <cell | n | todo | all>
#     <cell> a row of `grid` (SLURM_ARRAY_TASK_ID); fits all three methods
#     n      print nrow(grid) (for sbatch --array=1-<n>)
#     todo   print the rows with a method still to fit (for sbatch --array=<list>)
#     all    fit every row serially (local runs)
# Re-running skips methods already on disk
# set BBIS_OVERWRITE=TRUE to refit.

# setup workspace --------------------------------------------------------- ####
library(here)
# load custom functions
source(here("functions/utility_functions.R"))
sourceDir(here("functions"))                             # thin_track(), fit_*(), checkpoints
source(here("case_studies/reindeer/R/reindeer_data.R"))  # load_reindeer_tracks()

# Load necessary packages
load_lib(mvnfast, parallel, Rcpp,      # BBIS likelihood
         terra, dplyr, lubridate,      # data
         amt, survival, tidyr, purrr)  # iSSF

# experiment settings ----------------------------------------------------- ####
exp_name      <- "thinning_fixed_hour"
fits_dir      <- here("case_studies/reindeer/fitted_estimates/experiments", exp_name)
intervals     <- c(3, 6, 12, 24)  # training fix interval (h); 3 h = native GPS schedule
native_h      <- 3                # GPS fix interval (h): clock phases step by this
anchor_legacy <- 12               # clock hour (UTC) of the fits run before anchors were added
methods       <- c("ssf", "euler", "bbis")
smooth_km     <- c(bbis = 0.5, euler = 0.5, ssf = 0)  # covariate smoothing per method (km)

# BBIS
M           <- 200  # bridges
dt_max      <- 0.5  # max hours between nodes
max_retries <- 4    # refits with fresh bridges if the optimiser does not converge
lower_sigma <- 1e-4 # floor for the Euler gamma2 used as init_sigma

# iSSF
ssf_n_control <- 10     # random steps per used step
ssf_tol_frac  <- 1 / 6  # resampling tolerance, fraction of the interval

# compute
ncores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "10"))
overwrite <- as.logical(Sys.getenv("BBIS_OVERWRITE", unset = "FALSE"))
oneimpact_lib <- "/cluster/projects/nn11017k/BioDivMapping/R"  # oneimpact on Saga
if (!dir.exists(oneimpact_lib)) oneimpact_lib <- NULL

# load data --------------------------------------------------------------- ####
# raw covariates (iSSF) and the 0.5 km smoothed field (BBIS, Euler)
covs_raw    <- rast(here("case_studies/reindeer/data/reindeer_covariates.tiff"))
covs_smooth <- rast(here("case_studies/reindeer/data/smoothed/reindeer_covariates_smooth_0.5km.tif")) |>
  setNames(names(covs_raw))
# covariate field per smoothing distance (km); must cover the values in smooth_km
cov_fields <- list(`0` = covs_raw, `0.5` = covs_smooth)

# tracks on the covariate grid (km), native 3 h fixes; one fold per animal-year
tracks <- load_reindeer_tracks(covs_raw, lib = oneimpact_lib)
folds  <- sort(unique(tracks$ID))

# build grid: one row per training data set ------------------------------ ####
# every clock phase of each schedule: anchors 0, 3, ... below the interval
# (anchor_h and anchor_h + interval keep the same fixes)
phases <- bind_rows(lapply(intervals, function(i) {
  data.frame(interval = i, anchor_h = seq(0, i - native_h, by = native_h))
}))
grid <- crossing(phases, fold = folds) %>%
  mutate(seed   = 1000 * match(fold, folds) + 10 * interval + anchor_h / native_h,  # bridges / random steps
         data   = sprintf("int%g_a%02d_fold%s", interval, anchor_h, fold),
         legacy = anchor_h == anchor_legacy %% interval) %>%  # same data as an old 12:00 UTC fit
  as.data.frame()

#' Checkpoint label of one method's fit to one data set (the old label if a
#' fit from before anchors were added is on disk for the same data)
#' @param row grid row; @param method "ssf", "euler" or "bbis"
#' @return label
fit_label <- function(row, method) {
  old <- sprintf("%s_int%g_fold%s", method, row$interval, row$fold)
  if (row$legacy && checkpoint_exists(fits_dir, old)) old else paste(method, row$data, sep = "_")
}

#' Methods still to fit for one data set
#' @param row grid row
#' @return subset of methods
to_fit <- function(row) {
  if (overwrite) return(methods)
  methods[!vapply(methods, function(m) checkpoint_exists(fits_dir, fit_label(row, m)), logical(1))]
}

#' Fit one method to one data set and write its checkpoint
#'
#' @param row    grid row (the data set)
#' @param method "ssf", "euler" or "bbis"
#' @param train  thinned training track
#' @param expr   fitting code returning list(beta, sigma, loglik, convergence,
#'               converged, n_attempts, time_secs, ess); evaluated here, so an
#'               error is caught and stored as a failed fit
save_fit <- function(row, method, train, expr) {
  label <- fit_label(row, method)
  fit <- tryCatch(expr, error = function(err) {
    warning(sprintf("%s failed: %s", label, conditionMessage(err)))
    list(beta = setNames(rep(NA_real_, nlyr(covs_raw)), names(covs_raw)),
         sigma = NA_real_, loglik = NA_real_, convergence = NA_real_, converged = FALSE,
         n_attempts = NA_real_, time_secs = NA_real_, ess = ess_summary(NULL),
         error = conditionMessage(err))
  })
  # save the data set's description with the fit, so each checkpoint describes itself
  fit$beta <- fit$beta[names(covs_raw)]
  write_checkpoint(c(list(method = method, interval = row$interval, anchor_h = row$anchor_h,
                          fold = row$fold, smooth_km = unname(smooth_km[method]),
                          M = if (method == "bbis") M else NA, seed = row$seed,
                          label = label, n_train = nrow(train)),
                     fit),
                   fits_dir, label)
  message(sprintf("  %s written", label))
}

# pick rows to fit -------------------------------------------------------- ####
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: reindeer_fit.R <cell | n | todo | all>")
if (args[1] == "n") {
  cat(nrow(grid), "\n")
  quit(save = "no")
}
if (args[1] == "todo") {
  todo_rows <- which(vapply(seq_len(nrow(grid)), function(i) length(to_fit(grid[i, ])) > 0, logical(1)))
  # compress consecutive rows to ranges for sbatch, e.g. 25-36,49-132
  runs <- split(todo_rows, cumsum(c(1, diff(todo_rows) != 1)))
  cat(paste(vapply(runs, function(r) if (length(r) > 1) sprintf("%d-%d", min(r), max(r))
                   else as.character(r), ""), collapse = ","), "\n", sep = "")
  quit(save = "no")
}
cells <- if (args[1] == "all") seq_len(nrow(grid)) else as.integer(args[1])
if (anyNA(cells) || any(cells < 1 | cells > nrow(grid)))
  stop(sprintf("cell must be in 1..%d (got '%s')", nrow(grid), args[1]))

# fit --------------------------------------------------------------------- ####
for (i in cells) {
  row  <- grid[i, ]
  todo <- to_fit(row)
  if (!length(todo)) {
    message(sprintf("skip %s (all methods on disk)", row$data))
    next
  }
  message(sprintf("== %s | cell %d/%d | %s | %d cores ==", row$data, i, nrow(grid),
                  paste(todo, collapse = ", "), ncores))

  # withhold one animal-year, then thin the training tracks to the fixed-hour schedule
  train <- thin_track(tracks[tracks$ID != row$fold, ], row$interval, row$anchor_h)

  # iSSF on the raw covariates; step duration = the thinned interval
  message("fitting SSF")
  if ("ssf" %in% todo) save_fit(row, "ssf", train, {
    s <- fit_ssf(train, cov_fields[[as.character(smooth_km[["ssf"]])]],
                 resample_rate = row$interval, tolerance = row$interval * ssf_tol_frac,
                 n_control = ssf_n_control, seed = row$seed)
    list(beta = s$beta, sigma = NA_real_, loglik = s$loglik, convergence = 0,
         converged = all(is.finite(s$beta)), n_attempts = 1, time_secs = s$time_secs,
         ess = ess_summary(NULL))
  })

  # Euler; its gamma2 (hours) also seeds init_sigma for BBIS, so fit it for either
  message("fitting euler")
  euler <- if (any(c("euler", "bbis") %in% todo))
    tryCatch(fit_euler(train, cov_fields[[as.character(smooth_km[["euler"]])]]),
             error = function(err) err)
  if ("euler" %in% todo) save_fit(row, "euler", train, {
    if (inherits(euler, "error")) stop(euler)
    list(beta = euler$beta, sigma = euler$gamma2, loglik = euler$loglik, convergence = 0,
         converged = TRUE, n_attempts = 1, time_secs = euler$time_secs,
         ess = ess_summary(NULL))
  })

  # BBIS
  message("fitting BBIS")
  if ("bbis" %in% todo) save_fit(row, "bbis", train, {
    init_sigma <- if (!inherits(euler, "error") && is.finite(euler$gamma2) &&
                      euler$gamma2 > lower_sigma) euler$gamma2 else 1
    r <- fit_langevin_bbis_retry(train, cov_fields[[as.character(smooth_km[["bbis"]])]],
                                 M = M, dt_max = dt_max, dt_units = "hours", seed = row$seed,
                                 ncores = ncores, init_sigma = init_sigma,
                                 max_retries = max_retries, diagnostics = TRUE)
    if (is.null(r$fit)) stop("all BBIS attempts failed")
    list(beta = r$fit$beta, sigma = r$fit$sigma, loglik = r$fit$loglik,
         convergence = r$fit$convergence, converged = r$converged,
         n_attempts = r$n_attempts, time_secs = as.numeric(r$fit$time, units = "secs"),
         ess = ess_summary(r$fit$diagnostics))
  })
}
message("== done ==")

