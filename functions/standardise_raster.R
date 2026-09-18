# Standardise SpatRaster by layer (mean 0, sd 1)
standardise_raster <- function(r) {
  if (!inherits(r, "SpatRaster")) {
    stop("Input must be a SpatRaster.")
  }
  # apply mean and sd to each layer
  m <- unlist(terra::global(r, "mean", na.rm = TRUE))
  s <- unlist(terra::global(r, "sd",   na.rm = TRUE))
  s[s == 0 | is.na(s)] <- 1
  (r - m) / s
}