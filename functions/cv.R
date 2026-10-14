#' Jackknife standard error from leave-one-out estimates
#'
#' Leave-one-out fits share all but one group, so their spread understates the
#' uncertainty of the full-data estimate; the jackknife scales it back up:
#' SE = sqrt((n - 1) / n * sum((x_i - mean(x))^2)) = sd(x) * (n - 1) / sqrt(n)
#' @param x one estimate per fold (non-finite values dropped)
#' @return scalar SE (NA if fewer than 2 finite estimates)
jack_se <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  sqrt((n - 1) / n * sum((x - mean(x))^2))
}

#' Sort a track by ID/time and (re)compute per-ID step duration and length
#' @param track data.frame with x, y, time, ID
#' @param dt_units units for `dt` (default hours)
#' @return track with `dt` (time to next fix) and `step` (distance to next
#'         fix); NA on the last fix of each ID
track_steps <- function(track, dt_units = "hours") {
  track <- track[order(track$ID, track$time), , drop = FALSE]
  n    <- nrow(track)
  same <- c(track$ID[-1] == track$ID[-n], FALSE)
  nx   <- c(track$x[-1], NA); ny <- c(track$y[-1], NA)
  nt   <- c(track$time[-1], track$time[n])
  track$dt   <- ifelse(same, as.numeric(difftime(nt, track$time, units = dt_units)), NA)
  track$step <- ifelse(same, sqrt((nx - track$x)^2 + (ny - track$y)^2), NA)
  rownames(track) <- NULL
  track
}

#' Thin a track set to a fixed-clock-hour schedule, per ID
#'
#' Mimics a tag programmed to fix every `interval_h` hours: keeps only fixes at
#' clock time anchor_h + k * interval_h, so every fix of a 24 h schedule is at
#' the same time of day (e.g. 24 h -> 12:00 only; 12 h -> 00:00 & 12:00;
#' 6 h -> 00, 06, 12, 18). Coarser schedules are nested in finer ones. A missed
#' fix stays a gap, as it would on a real tag. Clock time is read in the time
#' zone stored on `track$time` (UTC for the reindeer).
#' @param track      data.frame with x, y, time (POSIXct), ID
#' @param interval_h target fix interval in hours (must divide 24)
#' @param anchor_h   clock hour every schedule includes (default 12)
#' @return thinned track with recomputed dt (hours) and step
thin_track <- function(track, interval_h, anchor_h = 12) {
  if (24 %% interval_h != 0)
    stop(sprintf("interval_h = %g does not divide 24 h", interval_h))
  # minutes since midnight of each fix
  lt  <- as.POSIXlt(track$time)
  mod <- lt$hour * 60 + lt$min
  # keep fixes that fall exactly on the schedule
  keep <- lt$sec == 0 & (mod - anchor_h * 60) %% (interval_h * 60) == 0
  track_steps(track[keep, , drop = FALSE])
}
