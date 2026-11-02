#' Fit the Euler-discretised Langevin model to a track
#' Wraps langevinUD() / nllkEuler() from Rhabit_functions.R with the same 
#' track/covariate interface as fit_langevin_bbis()
#'
#' @param track      data.frame with x, y, time (POSIXct), ID; consecutive rows
#'                   of the same ID form steps (sorted by ID, time)
#' @param covlist    covariate list (x, y, z) or SpatRaster
#' @param time_units units for time increments (and hence gamma2): "secs",
#'                   "mins", "hours", "days". Use the same as BBIS dt_units.
#' @return list(beta (named), gamma2, sigma (= gamma2), loglik, convergence (0),
#'              time_secs)
fit_euler <- function(track, covlist, time_units = "hours") {
  if (inherits(covlist, "SpatRaster")) covlist <- spatRast_to_covlist(covlist)
  t0    <- Sys.time()
  locs  <- as.matrix(track[, c("x", "y")])
  times <- as.numeric(difftime(track$time, min(track$time), units = time_units))
  ID    <- if ("ID" %in% names(track)) track$ID else rep(1, nrow(track))
  ga    <- bilinearGradArray(locs, covlist)
  e     <- langevinUD(locs, times, ID = ID, grad_array = ga)
  beta  <- setNames(as.numeric(e$betaHat), names(covlist))
  ll    <- -nllkEuler(beta = beta, gamma2 = e$gamma2Hat, locs = locs, times = times,
                      ID = ID, grad_array = ga)
  list(beta = beta, gamma2 = e$gamma2Hat, sigma = e$gamma2Hat, loglik = ll,
       convergence = 0L, time_secs = as.numeric(Sys.time() - t0, units = "secs"))
}
