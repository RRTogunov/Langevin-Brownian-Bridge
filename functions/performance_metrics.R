# Performance metrics for predicted utilisation distributions (UDs)
# Each scoring metric has
#   * a vector core  `*_score(u, a)`  (u = UD at used fixes, a = UD at available
#     cells), so eval_ud() extracts the raster values once and reuses them; and
#   * a raster wrapper (boyce_continuous, boyce_index, auc_ud, tss_ud,
#     logdens_ud) with the signature used by the old reindeer scripts.
# AUC, TSS and (log-scale) CBI only use the ranking of cells, so they cannot see
# the size of the betas (UD^k ranks cells identically for any k > 0); logdens
# and kstar can

#' Extract UD values at used fixes and at available locations
#'
#' @param pred_r    SpatRaster (single layer) predicted UD
#' @param used_pts  SpatVector of used fixes
#' @param available NULL (all finite cells of pred_r), a SpatVector of points,
#'                  or a numeric vector of UD values
#' @return list(u = finite used values, a = finite available values)
ud_values <- function(pred_r, used_pts, available = NULL) {
  u <- terra::extract(pred_r, used_pts, ID = FALSE)[, 1]
  a <- if (is.null(available)) {
    terra::values(pred_r, mat = FALSE)
  } else if (inherits(available, "SpatVector")) {
    terra::extract(pred_r, available, ID = FALSE)[, 1]
  } else {
    as.numeric(available)
  }
  list(u = u[is.finite(u)], a = a[is.finite(a)])
}

#' AUC (Mann-Whitney): P(UD at a used fix > UD at an available cell)
#' @param u used values; @param a available values
#' @return scalar AUC in [0, 1]
auc_score <- function(u, a) {
  if (!length(u) || !length(a)) return(NA_real_)
  r  <- rank(c(u, a))
  n1 <- as.numeric(length(u)); n0 <- as.numeric(length(a))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' TSS: max over thresholds of sensitivity + specificity - 1 (used vs available)
#' @param u used values; @param a available values
#' @param nthr number of thresholds (quantiles of the pooled values)
#' @return scalar TSS in [0, 1]
tss_score <- function(u, a, nthr = 200) {
  if (!length(u) || !length(a)) return(NA_real_)
  thr  <- stats::quantile(c(u, a), probs = seq(0, 1, length.out = nthr), na.rm = TRUE)
  sens <- vapply(thr, function(t) mean(u >= t), numeric(1))   # true-positive rate
  spec <- vapply(thr, function(t) mean(a <  t), numeric(1))   # true-negative rate
  max(sens + spec - 1, na.rm = TRUE)
}

#' Continuous Boyce index (CBI; Hirzel et al. 2006)
#'
#' Moving window of width `window` (fraction of the predicted range) slid in
#' `res` steps; Spearman correlation of the P/E ratio against window midpoints
#' (as ecospat::ecospat.boyce()).
#' @param u used values; @param a available values
#' @param window window width as fraction of range; @param res number of steps
#' @return scalar CBI in [-1, 1]
cbi_score <- function(u, a, window = 0.1, res = 100) {
  if (!length(u) || !length(a)) return(NA_real_)
  mini <- min(a); maxi <- max(a)
  window.w <- window * (maxi - mini)
  vec.mov  <- seq(mini, maxi - window.w, length.out = res + 1)
  Fi <- vapply(vec.mov, function(lo) {
    mean(u >= lo & u <= lo + window.w) / mean(a >= lo & a <= lo + window.w)
  }, numeric(1))
  mids <- vec.mov + window.w / 2
  keep <- is.finite(Fi)
  if (sum(keep) < 2) return(NA_real_)
  suppressWarnings(stats::cor(Fi[keep], mids[keep], method = "spearman"))
}

#' Original binned Boyce index (OBI; Boyce et al. 2002)
#' @param u used values; @param a available values; @param nclass number of bins
#' @return scalar OBI in [-1, 1]
obi_score <- function(u, a, nclass = 20) {
  if (!length(u) || !length(a)) return(NA_real_)
  brks <- seq(min(a), max(a), length.out = nclass + 1)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  Pi   <- as.numeric(table(cut(a, brks, include.lowest = TRUE))) / length(a)
  Ei   <- as.numeric(table(cut(u, brks, include.lowest = TRUE))) / length(u)
  Fi   <- Ei / Pi
  keep <- is.finite(Fi) & Pi > 0
  suppressWarnings(stats::cor(mids[keep], Fi[keep], method = "spearman"))
}

#' Mean log predicted UD at used fixes (magnitude-sensitive; higher = better)
#' @param u used values
#' @return scalar
logdens_score <- function(u) {
  u <- u[u > 0]
  if (!length(u)) return(NA_real_)
  mean(log(u))
}

#' Calibration slope of a UD (k*) and the log score it attains
#'
#' Rescales the UD to UD^k (= exp(k * beta . c): every beta times k), renormalised
#' over the available cells, and finds the k maximising the mean log UD at the
#' used fixes. k* = 1: betas the right size; k* > 1: too small (UD too flat);
#' k* < 1: too large. The log score at k* depends only on the direction of beta.
#' At k = 1 it equals logdens_score() when the UD sums to 1 over `a`.
#' @param u used values; @param a available values (UD^k is renormalised over these)
#' @param k_range search interval for k
#' @return c(k = k*, logdens = mean log UD^k* at the used fixes)
kstar_fit <- function(u, a, k_range = c(0, 50)) {
  u <- u[u > 0]; a <- a[a > 0]
  if (!length(u) || !length(a)) return(c(k = NA_real_, logdens = NA_real_))
  lu <- log(u); la <- log(a)
  # log score of UD^k; log-sum-exp keeps the normalising constant finite
  ls  <- function(k) { m <- max(k * la); k * mean(lu) - (m + log(sum(exp(k * la - m)))) }
  opt <- stats::optimize(ls, k_range, maximum = TRUE)
  c(k = opt$maximum, logdens = opt$objective)
}

#' Score one predicted UD against used fixes (used-vs-available)
#'
#' @param ud        SpatRaster predicted UD (normalised or not; rank metrics
#'                  are scale-invariant, logdens assumes normalised)
#' @param used_pts  SpatVector of used (e.g. withheld) fixes
#' @param available see ud_values(); NULL = all finite cells of `ud`
#' @param metrics   any of "auc", "cbi", "tss", "logdens", "obi", "kstar"
#'                  (calibration slope) and "logdens_kstar" (log score at k*)
#' @param cbi_window,cbi_res CBI window settings (defaults = old eval_ud())
#' @param obi_nclass OBI bins; @param tss_nthr TSS thresholds
#' @return named numeric vector, one element per metric (NA if not computable)
eval_ud <- function(ud, used_pts, available = NULL,
                    metrics = c("auc", "cbi", "tss", "logdens"),
                    cbi_window = 0.025, cbi_res = 100, obi_nclass = 100, tss_nthr = 200) {
  v  <- ud_values(ud, used_pts, available)
  # k* and its log score come from one optimisation
  ks <- if (any(c("kstar", "logdens_kstar") %in% metrics)) kstar_fit(v$u, v$a)
  vapply(metrics, function(m) switch(m,
                                     auc     = auc_score(v$u, v$a),
                                     cbi     = cbi_score(v$u, v$a, cbi_window, cbi_res),
                                     tss     = tss_score(v$u, v$a, tss_nthr),
                                     logdens = logdens_score(v$u),
                                     obi     = obi_score(v$u, v$a, obi_nclass),
                                     kstar         = ks[["k"]],
                                     logdens_kstar = ks[["logdens"]],
                                     stop(sprintf("unknown metric '%s'", m))), numeric(1))
}