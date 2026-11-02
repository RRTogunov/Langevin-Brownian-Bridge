#' Load the reindeer GPS tracks (oneimpact::reindeer) on the covariate grid (km)
#'
#' @param covs SpatRaster whose CRS the tracks are projected to
#' @param lib  library holding the oneimpact package (NULL = default library)
#' @param crs  CRS of the raw track coordinates
#' @return data.frame x, y (km), time (POSIXct, UTC, rounded to the hour), ID
#'         (animal_year_id), dt (h), step (km); sorted by ID, time
load_reindeer_tracks <- function(covs, lib = NULL, crs = "EPSG:25833") {
  e <- new.env()
  utils::data("reindeer", package = "oneimpact", lib.loc = lib, envir = e)
  rd <- as.data.frame(e$reindeer)
  # project to the covariate grid (km)
  v  <- terra::vect(rd, geom = c("x", "y"), crs = crs) |>
    terra::project(terra::crs(covs))
  xy <- terra::crds(v)
  tr <- data.frame(x = xy[, 1], y = xy[, 2],
                   time = lubridate::round_date(rd$t, unit = "hour"),
                   ID = rd$animal_year_id)
  track_steps(tr, dt_units = "hours")
}
