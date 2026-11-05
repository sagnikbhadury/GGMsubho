# Step-up selection over the ranked, nested models.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

#' Log marginal likelihood of the working model
#'
#' For the design `cbind(1, X)` with the conjugate working prior
#' \eqn{\beta \mid \sigma^2 \sim N(0, \sigma^2 v I)},
#' \eqn{\sigma^2 \sim IG(c_0, d_0)}, returns
#' \deqn{\log p(y) = c_0\log d_0 + \log\Gamma(c_0 + n/2) - \log\Gamma(c_0) -
#'   (n/2)\log(2\pi) - \frac12\log|V^*{}^{-1}| - \frac{m}{2}\log v - c^*\log d^*,}
#' where \eqn{V^*{}^{-1} = v^{-1}I + X^TX}, \eqn{m} is the number of columns,
#' \eqn{c^* = c_0 + n/2} and \eqn{d^* = d_0 + (y^Ty - y^TXV^*X^Ty)/2}. With `v = 1` this
#' is the marginal likelihood in Bhadury et al. (2023). It is returned on the log scale
#' because the likelihood itself falls below the smallest representable number once `n`
#' reaches several hundred.
#'
#' @param y Numeric response.
#' @param X Predictor matrix without an intercept, or `NULL` for the intercept-only
#'   model.
#' @param c0,d0 Shape and rate of the inverse-gamma working prior.
#' @param v Prior variance multiplier of the working prior.
#' @return A single number.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 2), 40, 2)
#' y <- X[, 1] + rnorm(40)
#' log_marginal_working(y, X)
#' log_marginal_working(y, NULL)
#' @export
log_marginal_working <- function(y, X = NULL, c0 = 1, d0 = 1, v = 1) {
  y <- as.numeric(y)
  design <- if (is.null(X)) matrix(1, length(y), 1L) else cbind(1, as.matrix(X))
  .check_working(c0, d0, v)
  m <- ncol(design)
  r <- step_up_path_cpp(crossprod(design), as.numeric(crossprod(design, y)), sum(y^2),
                        length(y), seq_len(m) - 1L, c0, d0, v, numeric(m))
  r$log_marginal[m]
}

.check_working <- function(c0, d0, v) {
  if (!all(is.finite(c(c0, d0, v))) || c0 <= 0 || d0 <= 0 || v <= 0)
    stop("`c0`, `d0` and `v` must be positive.", call. = FALSE)
  invisible(TRUE)
}

# log prior of a model with k of the P candidate predictors
.log_model_prior <- function(k, P, model_prior, a_model, b_model) {
  if (model_prior == "beta-binomial") {
    return(lgamma(a_model + b_model) - lgamma(a_model) - lgamma(b_model) +
             lgamma(a_model + k) + lgamma(P + b_model - k) -
             lgamma(a_model + b_model + P))
  }
  w <- if (P > 1L) 1 / P else 0.5
  k * log(w) + (P - k) * log1p(-w)
}

#' Step-up model selection along a ranking of the predictors
#'
#' Ranks the predictors by `abs(beta_hat)` and compares the `K + 1` nested models that
#' contain the `k` highest ranked predictors, `k = 0, ..., K`, always with an intercept.
#' Each model is scored by its working marginal likelihood
#' ([log_marginal_working()]) times a model prior, and the highest score wins. This is
#' the selection step of SUBHO (Bhadury et al. 2023).
#'
#' @param y Numeric response.
#' @param X Predictor matrix without an intercept.
#' @param beta_hat Numeric vector of length `ncol(X)`, typically the posterior median or
#'   mean of [horseshoe_regression()].
#' @param K Largest number of predictors considered. `NULL` means `ncol(X)`.
#' @param model_prior `"beta-binomial"` (the marginal of the Beta-Bernoulli prior) or
#'   `"binomial"` (independent inclusion with probability `1 / ncol(X)`).
#' @param a_model,b_model Parameters of the Beta prior on the inclusion probability.
#' @param c0,d0,v Working prior; see [log_marginal_working()].
#'
#' @return A list with `selected` (logical vector), `k`, `coef` (posterior mean under the
#'   selected working model, zero elsewhere), `log_posterior` and `log_marginal` (one
#'   entry for each `k`), and `order` (predictor indices from highest to lowest rank).
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(80 * 8), 80, 8)
#' y <- 1.5 * X[, 1] - X[, 2] + rnorm(80)
#' fit <- horseshoe_regression(y, X, n_iter = 400, burn = 100, thin = 2)
#' step_up_select(y, X, fit$beta_median)$k
#' @export
step_up_select <- function(y, X, beta_hat, K = NULL,
                           model_prior = c("beta-binomial", "binomial"),
                           a_model = 1, b_model = 1, c0 = 1, d0 = 1, v = 1) {
  model_prior <- match.arg(model_prior)
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (length(beta_hat) != ncol(X)) stop("`beta_hat` must have length ncol(X).", call. = FALSE)
  .check_working(c0, d0, v)
  design <- cbind(1, X)
  .step_up_core(crossprod(design), as.numeric(crossprod(design, y)), sum(y^2), length(y),
                as.numeric(beta_hat), seq_len(ncol(X)) + 1L, K, model_prior, a_model,
                b_model, c0, d0, v)
}

# `columns` holds the positions in the Gram matrix of the candidate predictors,
# `beta_hat` their estimates, in the same order.
.step_up_core <- function(gram, xty, yty, n, beta_hat, columns, K, model_prior, a_model,
                          b_model, c0, d0, v) {
  P <- length(columns)
  if (is.null(K)) K <- P
  if (length(K) != 1L || !is.numeric(K) || K < 0 || K != round(K))
    stop("`K` must be NULL or a non-negative whole number.", call. = FALSE)
  K <- min(K, P)
  ranks <- order(abs(beta_hat), decreasing = TRUE)
  path <- c(1L, columns[ranks[seq_len(K)]]) - 1L
  lp <- .log_model_prior(0:K, P, model_prior, a_model, b_model)
  r <- step_up_path_cpp(gram, xty, yty, n, path, c0, d0, v, lp)
  selected <- logical(P)
  coef <- numeric(P)
  if (r$k > 0L) {
    keep <- ranks[seq_len(r$k)]
    selected[keep] <- TRUE
    coef[keep] <- r$coef[-1L]
  }
  list(selected = selected, k = r$k, coef = coef,
       log_posterior = as.numeric(r$log_posterior),
       log_marginal = as.numeric(r$log_marginal), order = ranks)
}
