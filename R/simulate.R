# Data simulated from a known sparse Gaussian graphical model.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

#' Simulate data from a sparse Gaussian graphical model
#'
#' Generates a zero-mean multivariate normal sample whose precision matrix has a known
#' sparsity pattern.
#'
#' * `"ar1"`: \eqn{\Sigma_{ab} = \rho^{|a-b|}}, whose graph is a chain. One of the two
#'   designs in Bhadury et al. (2023).
#' * `"gwishart"`: a Bernoulli(`edge_prob`) graph with
#'   \eqn{\Omega \sim W_G(3, I)}, rescaled so that \eqn{\Sigma} has a unit diagonal. The
#'   other design in Bhadury et al. (2023). Requires the suggested package 'BDgraph'. Most
#'   partial correlations are small under this design.
#' * `"random"`: a convenience design for examples, not used in the article. A
#'   Bernoulli(`edge_prob`) graph with partial correlations of random sign and magnitude
#'   proportional to `signal`, with the diagonal set just above the absolute value of the
#'   smallest eigenvalue so that the precision matrix is positive definite.
#'
#' @param n Number of samples.
#' @param p Number of variables.
#' @param structure `"ar1"`, `"gwishart"` or `"random"`.
#' @param rho Correlation parameter for `"ar1"`.
#' @param edge_prob Edge probability for `"gwishart"` and `"random"`.
#' @param signal Range of the off-diagonal magnitudes before the rescaling, for `"random"`.
#' @param node_names Optional variable names. By default `V1`, `V2`, ...
#' @param seed Optional seed. The caller's random stream is restored afterwards when a
#'   seed is given.
#'
#' @return A list with `data` (an `n` by `p` matrix), `graph` (true adjacency matrix),
#'   `partial_cor`, `omega`, `sigma` and `true_beta` (see [true_beta()]).
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 75, p = 20, structure = "random", edge_prob = 0.05, seed = 1)
#' dim(sim$data)
#' sum(sim$graph[upper.tri(sim$graph)])
#' @export
simulate_ggm <- function(n, p, structure = c("ar1", "gwishart", "random"), rho = 0.7,
                         edge_prob = 0.1, signal = c(0.3, 0.6), node_names = NULL,
                         seed = NULL) {
  structure <- match.arg(structure)
  if (length(n) != 1L || length(p) != 1L || n < 2 || p < 2)
    stop("`n` and `p` must be whole numbers of at least 2.", call. = FALSE)
  if (!is.null(node_names) && length(node_names) != p)
    stop("`node_names` must have length p.", call. = FALSE)
  run <- function() {
    if (structure == "ar1") {
      sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
      omega <- solve(sigma)
      graph <- 1 * (abs(omega) > 1e-12)
    } else if (structure == "gwishart") {
      .require_suggested("BDgraph", "simulate_ggm(structure = \"gwishart\")")
      g <- .random_graph(p, edge_prob)
      w <- BDgraph::rgwish(n = 1, adj = g + diag(p), b = 3, D = NULL, threshold = 1e-8)
      s <- solve(w)
      sigma <- stats::cov2cor(s)
      omega <- solve(sigma) * (g + diag(p))
      graph <- g
    } else {
      g <- .random_graph(p, edge_prob)
      vals <- stats::runif(sum(g) / 2, signal[1L], signal[2L]) *
        sample(c(-1, 1), sum(g) / 2, replace = TRUE)
      omega <- matrix(0, p, p)
      omega[upper.tri(omega)][g[upper.tri(g)] == 1] <- vals
      omega <- omega + t(omega)
      diag(omega) <- abs(min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values)) + 0.2
      sigma <- stats::cov2cor(solve(omega))
      omega <- solve(sigma) * (g + diag(p))
      graph <- g
    }
    diag(graph) <- 0
    data <- .rmvnorm_chol(n, sigma)
    list(data = data, graph = graph, sigma = sigma, omega = omega)
  }
  res <- if (is.null(seed)) run() else withr::with_seed(.resolve_seed(seed), run())
  nm <- if (is.null(node_names)) paste0("V", seq_len(p)) else as.character(node_names)
  psi <- .psi_from_omega(res$omega)
  psi[res$graph == 0 & row(psi) != col(psi)] <- 0
  add_names <- function(m) {
    dimnames(m) <- list(nm, nm)
    m
  }
  data <- res$data
  colnames(data) <- nm
  list(data = data, graph = add_names(res$graph), partial_cor = add_names(psi),
       omega = add_names(res$omega), sigma = add_names(res$sigma),
       true_beta = add_names(true_beta(res$sigma)))
}

.random_graph <- function(p, edge_prob) {
  g <- matrix(0, p, p)
  g[upper.tri(g)] <- stats::rbinom(p * (p - 1) / 2, 1L, edge_prob)
  g + t(g)
}

.rmvnorm_chol <- function(n, sigma) {
  matrix(stats::rnorm(n * nrow(sigma)), n, nrow(sigma)) %*% chol(sigma)
}

#' True regression coefficients of a covariance matrix
#'
#' Computes \eqn{\Sigma_{aB}\Sigma_{BB}^{-1}} for every node.
#'
#' @param sigma Covariance matrix.
#' @return A square matrix with `beta[a, b]` the coefficient of \eqn{X_b} when \eqn{X_a} is
#'   the response.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' true_beta(0.5^abs(outer(1:4, 1:4, "-")))
#' @export
true_beta <- function(sigma) {
  p <- ncol(sigma)
  out <- matrix(0, p, p, dimnames = dimnames(sigma))
  for (a in seq_len(p)) out[a, -a] <- solve(sigma[-a, -a, drop = FALSE], sigma[-a, a])
  out
}
