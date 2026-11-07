# custom function crop and maintain aspect ratio
crop_aspect <- function(r, ymin = NULL, ymax = NULL,
                        anchor = c("right", "left", "center")) {
  anchor <- match.arg(anchor)
  ext   <- as.vector(ext(r))
  asp <- (ext[["ymax"]] - ext[["ymin"]]) / (ext[["xmax"]] - ext[["xmin"]])
  ymin <- if (is.null(ymin)) ext[["ymin"]] else ymin
  ymax <- if (is.null(ymax)) ext[["ymax"]] else ymax
  dx  <- (ymax - ymin) / asp
  xr  <- switch(anchor,
                right  = c(ext[["xmax"]] - dx, ext[["xmax"]]),
                left   = c(ext[["xmin"]], ext[["xmin"]] + dx),
                center = mean(ext[c("xmin", "xmax")]) + c(-1, 1) * dx / 2)
  crop(r, ext(xr[1], xr[2], ymin, ymax))
}
