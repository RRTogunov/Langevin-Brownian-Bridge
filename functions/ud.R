#' Normalise a raster to sum to 1
#' @param r SpatRaster (non-negative)
#' @return SpatRaster summing to 1 over finite cells
normalise_ud <- function(r) {
  r / terra::global(r, "sum", na.rm = TRUE)[[1]]
}

#' Stationary UD from habitat-selection coefficients
#'
#' @param beta      numeric vector named by covariate (names must match layers;
#'                  unnamed = taken in layer order)
#' @param cov_stack SpatRaster of covariates the betas were estimated on
#' @return single-layer SpatRaster, normalised to sum 1
ud_from_beta <- function(beta, cov_stack) {
  beta <- unlist(beta)
  b <- if (is.null(names(beta))) setNames(as.numeric(beta), names(cov_stack))
       else setNames(as.numeric(beta), names(beta))[names(cov_stack)]
  if (anyNA(b)) stop("ud_from_beta: betas missing for layer(s): ",
                     paste(names(cov_stack)[is.na(b)], collapse = ", "))
  normalise_ud(exp(sum(cov_stack * b)))
}

#' Mean UD across a set of fits (e.g. CV folds), renormalised
#'
#' @param beta_list  list of named beta vectors (non-finite ones are dropped)
#' @param cov_stacks one SpatRaster used for all, or a list parallel to beta_list
#' @return SpatRaster (sum 1), or NULL if no finite betas
mean_ud <- function(beta_list, cov_stacks) {
  if (inherits(cov_stacks, "SpatRaster"))
    cov_stacks <- rep(list(cov_stacks), length(beta_list))
  ok <- vapply(beta_list, function(b) all(is.finite(unlist(b))), logical(1))
  if (!any(ok)) return(NULL)
  uds <- terra::rast(Map(ud_from_beta, beta_list[ok], cov_stacks[ok]))
  normalise_ud(terra::mean(uds))
}
