#' GGMsubho: sparse Gaussian graphical models with step-up horseshoe selection
#'
#' Implements Graphical-SUBHO (Step-Up Bayesian HOrseshoe). The main function is
#' [subho()]. [plot_network()] draws the estimated network, and
#' [simulate_ggm()] generates data with a known graph.
#'
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' @keywords internal
#' @useDynLib GGMsubho, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

utils::globalVariables(c(
  "x", "y", "name", "strength", "sign_label", "magnitude", "edge_group",
  "row_node", "col_node", "value", "method_x", "method_y", "overlap"
))
