# Horseshoe regression with the Makalic-Schmidt / Bhattacharya Gibbs updates.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

#' Horseshoe regression
#'
#' Gibbs sampler for the linear regression
#' \eqn{y = X\beta + \epsilon}, \eqn{\epsilon \sim N(0, \sigma^2 I)}, with
#' \eqn{\beta_j \sim N(0, \lambda_j^2 \tau^2 \sigma^2)} and
#' \eqn{\lambda_j, \tau \sim C^+(0, 1)}. This is the regression fitted at every node
#' by [subho()]. The intercept, when requested, is the first column of the design and
#' carries the same prior, as in Bhadury et al. (2023).
#'
#' The coefficient update draws from \eqn{N(A^{-1}X^Ty, \sigma^2 A^{-1})} with a
#' Cholesky factor of \eqn{A} when the number of columns is below `0.8 * nrow(X)`, and
#' otherwise uses the algorithm of Bhattacharya et al. (2016), whose cost is cubic in
#' the number of rows.
#'
#' @param y Numeric response.
#' @param X Numeric predictor matrix, without an intercept column.
#' @param intercept Add an intercept column to the design?
#' @param n_iter Total Gibbs iterations.
#' @param burn Iterations discarded at the start.
#' @param thin Keep every `thin`-th iteration after `burn`.
#' @param sigma_prior `"sparseIG"` for \eqn{\sigma^2 \sim IG(0.1, 0.1)} or `"Jeffreys"`
#'   for \eqn{p(\sigma^2) \propto 1/\sigma^2}.
#' @param tau `NULL` (default) for the half-Cauchy prior on the global scale, or a fixed
#'   positive value of \eqn{\tau}.
#' @param seed Optional integer seed. The caller's random number stream is restored
#'   afterwards when a seed is given.
#'
#' @return A list of class `hs_fit`:
#' \describe{
#'   \item{beta_mean, beta_median}{Posterior mean and median of the predictor
#'     coefficients.}
#'   \item{intercept_draws}{Draws of the intercept, or `NULL`.}
#'   \item{beta_draws}{Matrix of draws, predictors in rows and saved iterations in
#'     columns.}
#'   \item{sigma2_draws, tau2_draws}{Draws of \eqn{\sigma^2} and \eqn{\tau^2}.}
#'   \item{pseudo_p}{\eqn{2\min\{P(\beta_j > 0), P(\beta_j < 0)\}} for each predictor.}
#' }
#'
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Carvalho CM, Polson NG, Scott JG (2010). The horseshoe estimator for sparse signals.
#' Biometrika, 97(2), 465-480.
#'
#' Makalic E, Schmidt DF (2015). A simple sampler for the horseshoe estimator. IEEE
#' Signal Processing Letters, 23(1), 179-182.
#'
#' Bhattacharya A, Chakraborty A, Mallick BK (2016). Fast sampling with Gaussian scale
#' mixture priors in high-dimensional regression. Biometrika.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 10), 60, 10)
#' y <- 2 * X[, 1] - X[, 2] + rnorm(60)
#' fit <- horseshoe_regression(y, X, n_iter = 400, burn = 100, thin = 2)
#' round(fit$beta_median, 2)
#' @export
horseshoe_regression <- function(y, X, intercept = TRUE, n_iter = 2000, burn = 200,
                                 thin = 10, sigma_prior = c("sparseIG", "Jeffreys"),
                                 tau = NULL, seed = NULL) {
  sigma_prior <- match.arg(sigma_prior)
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!is.numeric(X) || anyNA(X) || anyNA(y) || length(y) != nrow(X))
    stop("`y` and `X` must be numeric, without missing values, with length(y) == nrow(X).",
         call. = FALSE)
  .check_mcmc(n_iter, burn, thin)
  if (!is.null(tau) && (length(tau) != 1L || !is.numeric(tau) || tau <= 0))
    stop("`tau` must be NULL or a single positive number.", call. = FALSE)
  design <- if (intercept) cbind(1, X) else X
  run <- function() {
    .hs_draws(design, y, n_iter, burn, thin, sigma_prior, tau)
  }
  res <- if (is.null(seed)) run() else withr::with_seed(.resolve_seed(seed), run())
  draws <- res$beta
  if (intercept) {
    intercept_draws <- draws[1L, ]
    draws <- draws[-1L, , drop = FALSE]
  } else {
    intercept_draws <- NULL
  }
  rownames(draws) <- colnames(X)
  out <- list(
    beta_mean = stats::setNames(rowMeans(draws), colnames(X)),
    beta_median = stats::setNames(apply(draws, 1L, stats::median), colnames(X)),
    intercept_draws = intercept_draws,
    beta_draws = draws,
    sigma2_draws = as.numeric(res$sigma2),
    tau2_draws = as.numeric(res$tau2),
    pseudo_p = apply(draws, 1L, function(b) 2 * min(mean(b > 0), mean(b < 0)))
  )
  structure(out, class = "hs_fit")
}

.hs_draws <- function(design, y, n_iter, burn, thin, sigma_prior, tau) {
  a_sigma <- if (sigma_prior == "sparseIG") 0.1 else 0
  b_sigma <- if (sigma_prior == "sparseIG") 0.1 else 0
  hs_gibbs_cpp(design, y, as.integer(n_iter), as.integer(burn), as.integer(thin),
               a_sigma, b_sigma, 1, if (is.null(tau)) 1 else tau^2, !is.null(tau))
}

#' Coefficients kept by a credible interval
#'
#' Sets a posterior mean coefficient to zero when its central credible interval
#' contains zero, one of the selection rules discussed in Bhadury et al. (2023).
#'
#' @param fit A result of [horseshoe_regression()].
#' @param level Credible level, for example `0.5` or `0.95`.
#' @return The numeric vector of thresholded posterior means.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 5), 60, 5)
#' y <- 2 * X[, 1] + rnorm(60)
#' fit <- horseshoe_regression(y, X, n_iter = 300, burn = 100, thin = 2)
#' credible_beta(fit, level = 0.95)
#' @export
credible_beta <- function(fit, level = 0.95) {
  if (!inherits(fit, "hs_fit")) stop("`fit` must come from horseshoe_regression().", call. = FALSE)
  if (length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be between 0 and 1.", call. = FALSE)
  alpha <- 1 - level
  lo <- apply(fit$beta_draws, 1L, stats::quantile, probs = alpha / 2)
  hi <- apply(fit$beta_draws, 1L, stats::quantile, probs = 1 - alpha / 2)
  fit$beta_mean * (lo * hi > 0)
}

#' @export
print.hs_fit <- function(x, ...) {
  cat(sprintf("Horseshoe regression: %d predictors, %d saved draws\n",
              nrow(x$beta_draws), ncol(x$beta_draws)))
  cat(sprintf("Posterior mean of sigma^2: %.3f\n", mean(x$sigma2_draws)))
  invisible(x)
}
