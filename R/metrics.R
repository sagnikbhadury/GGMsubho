# Accuracy of an estimated network against a known one.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

.upper <- function(x) if (is.matrix(x)) x[upper.tri(x)] else x

.confusion <- function(est, truth) {
  e <- est != 0
  t <- truth != 0
  tp <- sum(e & t); tn <- sum(!e & !t); fp <- sum(e & !t); fn <- sum(!e & t)
  denom <- sqrt((tp + fp) * (tp + fn)) * sqrt((tn + fp) * (tn + fn))
  c(true_edges = tp + fn, edges = tp + fp, tp = tp, fp = fp, fn = fn, tn = tn,
    fdr = if (tp + fp > 0) fp / (tp + fp) else 0,
    fpr = if (fp + tn > 0) fp / (fp + tn) else 0,
    tpr = if (tp + fn > 0) tp / (tp + fn) else 0,
    mcc = if (denom > 0) (tp * tn - fp * fn) / denom else 0)
}

#' Accuracy of an estimated graph
#'
#' Compares two graphs over the pairs of nodes in the upper triangle.
#'
#' An FDR of 0 is returned when no edge is reported, and an MCC of 0 when its denominator
#' is 0.
#'
#' @param estimate Estimated adjacency matrix, or the upper triangle as a vector.
#' @param truth True adjacency matrix, or the upper triangle as a vector.
#' @return A named numeric vector: `true_edges`, `edges`, `tp`, `fp`, `fn`, `tn`, `fdr`,
#'   `fpr`, `tpr`, `mcc`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' truth <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
#' est <- rbind(c(0, 1, 1), c(1, 0, 0), c(1, 0, 0))
#' graph_metrics(est, truth)
#' @export
graph_metrics <- function(estimate, truth) {
  e <- .upper(estimate)
  t <- .upper(truth)
  if (length(e) != length(t)) stop("`estimate` and `truth` differ in size.", call. = FALSE)
  .confusion(e, t)
}

#' Error of estimated partial correlations
#'
#' Sums the squared errors in the upper triangle, split by whether the true value is zero.
#'
#' @param estimate,truth Matrices, or upper triangles as vectors.
#' @return A named vector: `error_zero`, `error_nonzero` and `error_total`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' partial_cor_error(c(0.1, 0.5, 0), c(0, 0.4, 0))
#' @export
partial_cor_error <- function(estimate, truth) {
  e <- .upper(estimate)
  t <- .upper(truth)
  if (length(e) != length(t)) stop("`estimate` and `truth` differ in size.", call. = FALSE)
  d <- (t - e)^2
  c(error_zero = sum(d[t == 0]), error_nonzero = sum(d[t != 0]), error_total = sum(d))
}

#' Accuracy of the largest partial correlations
#'
#' Keeps the `top` pairs with the largest absolute values and scores them as a graph.
#' Tied values are all kept.
#'
#' @param estimate Matrix of estimated values, or the upper triangle as a vector.
#' @param truth True adjacency matrix, or the upper triangle as a vector.
#' @param top Number of pairs to keep.
#' @return The result of [graph_metrics()], plus `threshold`, the smallest kept value.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' top_edge_metrics(c(0.9, 0.05, 0.4), c(1, 0, 0), top = 1)
#' @export
top_edge_metrics <- function(estimate, truth, top = 100) {
  e <- abs(.upper(estimate))
  t <- .upper(truth)
  sorted <- sort(e, decreasing = TRUE)
  top <- min(top, length(sorted))
  thr <- sorted[top]
  if (thr == 0) {
    nz <- sorted[sorted != 0]
    thr <- if (length(nz)) nz[length(nz)] else Inf
  }
  c(.confusion(1 * (e >= thr), t), threshold = thr)
}

#' Area under the ROC curve for edge detection
#'
#' Requires the suggested package 'pROC'.
#'
#' @param score Matrix of scores, or the upper triangle as a vector.
#' @param truth True adjacency matrix, or the upper triangle as a vector.
#' @return A single number.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examplesIf requireNamespace("pROC", quietly = TRUE)
#' edge_auc(c(0.9, 0.1, 0.4, 0.2), c(1, 0, 1, 0))
#' @export
edge_auc <- function(score, truth) {
  .require_suggested("pROC", "edge_auc()")
  as.numeric(pROC::auc(pROC::roc(response = .upper(truth), predictor = abs(.upper(score)),
                                 levels = c(0, 1), direction = "<", quiet = TRUE)))
}

#' Distribution of true partial and marginal correlations
#'
#' Counts the pairs in bins of absolute partial correlation, and the true edges in bins of
#' absolute marginal correlation.
#'
#' @param partial_cor True partial correlation matrix.
#' @param sigma True correlation matrix.
#' @return A list of two data frames, `partial` and `marginal`, with `bin` and `count`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 20, p = 10, structure = "ar1", seed = 1)
#' correlation_summary(sim$partial_cor, sim$sigma)
#' @export
correlation_summary <- function(partial_cor, sigma) {
  bins <- c(0, 0.1, 0.4, 0.7, Inf)
  labels <- c("(0, 0.1]", "(0.1, 0.4]", "(0.4, 0.7]", "> 0.7")
  count <- function(v) {
    data.frame(bin = c("0", labels),
               count = c(sum(v == 0), as.vector(table(cut(v[v != 0], bins, labels = labels)))))
  }
  pc <- abs(.upper(partial_cor))
  mc <- abs(.upper(sigma))
  list(partial = count(pc), marginal = count(mc * (pc != 0)))
}
