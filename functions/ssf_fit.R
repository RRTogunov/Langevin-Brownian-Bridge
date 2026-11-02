#' Fit an iSSF to a track set
#' Per ID: resample to a regular rate (amt), keep bursts with >= min_n fixes,
#' build steps, draw n_control random steps, extract covariates at step ends and
#' fit a conditional logistic model with habitat + movement terms.
#'
#' @param track         data.frame with x, y, time (POSIXct), ID
#' @param covs          SpatRaster of habitat covariates (evaluated at step ends)
#' @param resample_rate target step duration in hours
#' @param tolerance     resampling tolerance in hours
#' @param n_control     random (available) steps per used step
#' @param hab_terms     habitat covariate names (default all layers)
#' @param move_terms    movement-kernel terms (sl_, log_sl_, cos_ta_)
#' @param min_n         minimum fixes per burst
#' @param seed          optional RNG seed for random_steps()
#' @return list(beta (habitat, named by hab_terms), move (movement coefs),
#'              sigma (NA), loglik, convergence (0), n_steps, time_secs)
fit_ssf <- function(track, covs, resample_rate, tolerance, n_control = 10,
                    hab_terms = names(covs),
                    move_terms = c("sl_", "log_sl_", "cos_ta_"),
                    min_n = 3, seed = NULL) {
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)
  # Hours (numeric) to a lubridate Period, rounded to whole minutes
  rate <- lubridate::minutes(round(resample_rate * 60))
  tol <- lubridate::minutes(round(tolerance * 60))
  # prep data with amt, make tracks
  st <- amt::make_track(track, x, y, time, crs = terra::crs(covs), all_cols = TRUE) %>%
    tidyr::nest(.by = ID) %>%   # group by ID
    dplyr::mutate(steps = purrr::map(data, function(d) {
      d %>%
        amt::track_resample(rate = rate, tolerance = tol) %>%  # subsample to tol
        amt::filter_min_n_burst(min_n = min_n) %>%
        amt::steps_by_burst(keep_cols = "start") %>%
        amt::random_steps(n_control = n_control)  # generate control locs
    })) %>%
    dplyr::select(-data) %>% 
    tidyr::unnest(steps) %>%
    # add step_id for clogit strata
    dplyr::mutate(step_id = paste(ID, burst_, step_id_, sep = "_"))
  # print number of steps that will be modelled
  cat(length(unique(st$step_id)))
  # convert tracks to spatVect & extract covs
  ev <- terra::vect(as.matrix(st[, c("x2_", "y2_")]), type = "points",
                    crs = terra::crs(covs))
  st <- dplyr::bind_cols(st, terra::extract(covs[[hab_terms]], ev, ID = FALSE))
  # define step & turn angle vars for model
  st <- dplyr::mutate(st, log_sl_ = log(sl_), cos_ta_ = cos(ta_))
  # generate issf formula for clogit
  form <- stats::reformulate(c(hab_terms, move_terms, "strata(step_id)"),
                             response = "case_")
  m  <- survival::clogit(form, data = st)  # fit model
  cf <- stats::coef(m)  # extract covs 
  # combine outputs
  list(beta = setNames(as.numeric(cf[hab_terms]), hab_terms),  # named beta estimates
       move = cf[intersect(move_terms, names(cf))],  # movement coefs
       # diagnostics
       sigma = NA, loglik = as.numeric(m$loglik[2]), convergence = 0L,  
       n_steps = sum(st$case_),
       time_secs = as.numeric(Sys.time() - t0, units = "secs"))
}
