# Graphical-SUBHO: the node-wise driver and the post-processing.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

#' Estimate a sparse Gaussian graphical model with SUBHO
#'
#' Graphical-SUBHO (Step-Up Bayesian HOrseshoe; Bhadury, Mitra and Gaskins 2023) learns
#' a network by regressing every variable on all the others.
#'
#' For each node \eqn{a}:
#' 1. fit [horseshoe_regression()] of \eqn{X_a} on the other columns plus an intercept;
#' 2. rank the predictors by \eqn{|\hat\beta_a^{(b)}|};
#' 3. score the nested models holding the top \eqn{k} predictors, \eqn{k = 0, \ldots, K},
#'    by working marginal likelihood times the Beta-Bernoulli model prior
#'    ([step_up_select()]), and keep the best.
#'
#' The neighbourhoods are then combined by the AND and OR rules (Meinshausen and
#' Buhlmann 2006), and partial correlations come from [partial_correlation()].
#'
#' @section Choosing the arguments:
#' `n_iter = 2000`, `burn = 200`, `thin = 10`, the \eqn{IG(0.1, 0.1)} prior on the
#' residual variance, `c0 = d0 = v = 1` and `a_model = b_model = 1` are the settings of
#' the authors' code. `estimator = "median"` ranks by the posterior median; `"mean"` ranks
#' by the posterior mean. `K = NULL` evaluates the full path, whose cost per node is
#' cubic in `K`; a smaller `K` saves time when neighbourhoods are known to be small.
#'
#' @section Parallel computing:
#' Each node has its own seed, `seed + a`, so a fit with a given `seed` is identical for
#' every value of `cores`. Workers need the package to be installed. When `cores > 1`,
#' limit the threads of a multi-threaded BLAS so that they do not compete with the
#' workers.
#'
#' @param data Numeric matrix or data frame, samples in rows and variables in columns,
#'   with no missing values. Column names become node names.
#' @param estimator `"median"` or `"mean"`: the posterior summary of the horseshoe
#'   coefficients used for the ranking and for the partial correlations.
#' @param K Largest neighbourhood considered at each node. `NULL` means `ncol(data) - 1`.
#' @param standardize Centre and scale the columns first?
#' @param n_iter,burn,thin Gibbs iterations, burn-in and thinning at each node.
#' @param sigma_prior,tau Passed to [horseshoe_regression()].
#' @param model_prior,a_model,b_model,c0,d0,v Passed to [step_up_select()].
#' @param cores Number of parallel workers.
#' @param seed Integer seed, or `NULL` to draw one from the current random stream.
#' @param verbose Report progress (serial runs only)?
#'
#' @return An object of class `ggm_fit`, with
#' \describe{
#'   \item{graph_and, graph_or}{Symmetric 0/1 adjacency matrices.}
#'   \item{psi}{Partial correlations from the shrunken coefficients `beta`; non-zero
#'     outside the graph as well. Multiply by a graph to restrict it.}
#'   \item{psi_sparse}{Partial correlations from `beta_sparse`.}
#'   \item{gamma}{Selection matrix, with `gamma[a, b] = 1` when \eqn{X_b} is selected for
#'     response \eqn{X_a}; not symmetric.}
#'   \item{beta}{Coefficients from the horseshoe fit: `beta[a, b]` is the coefficient of
#'     \eqn{X_b} when \eqn{X_a} is the response.}
#'   \item{beta_sparse}{Coefficients of the selected working models.}
#'   \item{neighbourhood_size}{Number of predictors selected at each node.}
#'   \item{residual_variance}{Mean squared residual at each node under `beta`.}
#'   \item{sigma2}{Posterior mean of \eqn{\sigma_a^2} at each node.}
#'   \item{settings, nodes, n, p, seed, elapsed, call}{The run's details.}
#' }
#'
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Meinshausen N, Buhlmann P (2006). High-dimensional graphs and variable selection with
#' the lasso. The Annals of Statistics, 34(3), 1436-1462.
#'
#' Peng J, Wang P, Zhou N, Zhu J (2009). Partial correlation estimation by joint sparse
#' regression models. Journal of the American Statistical Association, 104(486), 735-746.
#'
#' @seealso [plot_network()], [graph_metrics()], [simulate_ggm()]
#'
#' @examples
#' sim <- simulate_ggm(n = 50, p = 8, structure = "ar1", seed = 1)
#' fit <- subho(sim$data, n_iter = 300, burn = 100, thin = 2, seed = 1)
#' fit
#' graph_metrics(fit$graph_and, sim$graph)
#' @export
subho <- function(data, estimator = c("median", "mean"), K = NULL, standardize = TRUE,
                  n_iter = 2000, burn = 200, thin = 10,
                  sigma_prior = c("sparseIG", "Jeffreys"), tau = NULL,
                  model_prior = c("beta-binomial", "binomial"), a_model = 1, b_model = 1,
                  c0 = 1, d0 = 1, v = 1, cores = 1, seed = NULL, verbose = FALSE) {
  cl <- match.call()
  estimator <- match.arg(estimator)
  sigma_prior <- match.arg(sigma_prior)
  model_prior <- match.arg(model_prior)
  .check_mcmc(n_iter, burn, thin)
  .check_working(c0, d0, v)
  cores <- .check_cores(cores)
  seed <- .resolve_seed(seed)
  Z <- .check_data(data, standardize)
  nodes <- colnames(Z)
  p <- ncol(Z)
  design_all <- cbind(1, Z)

  payload <- list(
    Z = Z, gram = crossprod(design_all), estimator = estimator, K = K,
    n_iter = n_iter, burn = burn, thin = thin, sigma_prior = sigma_prior, tau = tau,
    model_prior = model_prior, a_model = a_model, b_model = b_model,
    c0 = c0, d0 = d0, v = v
  )
  t0 <- proc.time()[["elapsed"]]
  res <- .run_nodes(p, .subho_node, payload, seed, cores, verbose)
  elapsed <- proc.time()[["elapsed"]] - t0

  beta <- beta_sparse <- gamma <- .empty_matrix(nodes)
  for (a in seq_len(p)) {
    beta[a, -a] <- res[[a]]$beta
    beta_sparse[a, -a] <- res[[a]]$beta_sparse
    gamma[a, -a] <- res[[a]]$selected
  }
  .new_fit(
    "SUBHO", nodes, nrow(Z),
    list(
      graph_and = symmetrize_graph(gamma, "and"),
      graph_or = symmetrize_graph(gamma, "or"),
      psi = partial_correlation(beta),
      psi_sparse = partial_correlation(beta_sparse),
      gamma = gamma, beta = beta, beta_sparse = beta_sparse,
      neighbourhood_size = stats::setNames(vapply(res, `[[`, integer(1), "k"), nodes),
      residual_variance = stats::setNames(vapply(res, `[[`, numeric(1), "residual"), nodes),
      sigma2 = stats::setNames(vapply(res, `[[`, numeric(1), "sigma2"), nodes),
      settings = list(estimator = estimator, K = K, standardize = standardize,
                      n_iter = n_iter, burn = burn, thin = thin,
                      sigma_prior = sigma_prior, tau = tau, model_prior = model_prior,
                      a_model = a_model, b_model = b_model, c0 = c0, d0 = d0, v = v)
    ),
    cl, elapsed, seed
  )
}

# One node. Exported to workers through the namespace.
.subho_node <- function(a, payload) {
  Z <- payload$Z
  n <- nrow(Z)
  y <- Z[, a]
  design <- cbind(1, Z[, -a, drop = FALSE])
  hs <- .hs_draws(design, y, payload$n_iter, payload$burn, payload$thin,
                  payload$sigma_prior, payload$tau)
  draws <- hs$beta[-1L, , drop = FALSE]
  beta <- if (payload$estimator == "median") apply(draws, 1L, stats::median) else rowMeans(draws)
  columns <- seq_len(ncol(Z) + 1L)[-c(1L, a + 1L)]
  xty <- as.numeric(crossprod(cbind(1, Z), y))
  s <- .step_up_core(payload$gram, xty, sum(y^2), n, beta, columns, payload$K,
                     payload$model_prior, payload$a_model, payload$b_model,
                     payload$c0, payload$d0, payload$v)
  list(beta = beta, beta_sparse = s$coef, selected = as.numeric(s$selected), k = s$k,
       residual = sum((y - mean(hs$beta[1L, ]) - design[, -1L, drop = FALSE] %*% beta)^2) / n,
       sigma2 = mean(hs$sigma2))
}

#' Combine neighbourhoods into an undirected graph
#'
#' The AND rule keeps edge \eqn{(a, b)} when each node selects the other; the OR rule
#' keeps it when either does (Meinshausen and Buhlmann 2006).
#'
#' @param gamma Square matrix; a non-zero entry means selected.
#' @param rule `"and"` or `"or"`.
#' @return A symmetric 0/1 matrix with a zero diagonal.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' g <- rbind(c(0, 1, 1), c(1, 0, 0), c(0, 0, 0))
#' symmetrize_graph(g, "and")
#' symmetrize_graph(g, "or")
#' @export
symmetrize_graph <- function(gamma, rule = c("and", "or")) {
  rule <- match.arg(rule)
  s <- gamma != 0
  out <- if (rule == "and") s & t(s) else s | t(s)
  storage.mode(out) <- "numeric"
  diag(out) <- 0
  dimnames(out) <- dimnames(gamma)
  out
}

#' Partial correlations from node-wise coefficients
#'
#' Computes
#' \eqn{\hat\psi_{ab} = \mathrm{sign}(\hat\beta_a^{(b)})
#' [\{\hat\beta_a^{(b)}\hat\beta_b^{(a)}\}^{1/2}]_{\le 1}} when the two coefficients have
#' the same sign, and zero otherwise (Peng et al. 2009; Bhadury et al. 2023).
#'
#' @param beta Square coefficient matrix with `beta[a, b]` the coefficient of \eqn{X_b}
#'   for response \eqn{X_a}.
#' @return A symmetric matrix with a unit diagonal. It need not be positive definite.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Peng J, Wang P, Zhou N, Zhu J (2009). Partial correlation estimation by joint sparse
#' regression models. Journal of the American Statistical Association, 104(486), 735-746.
#' @examples
#' b <- rbind(c(0, 0.5, -0.1), c(0.3, 0, 0.2), c(0.1, 0.4, 0))
#' partial_correlation(b)
#' @export
partial_correlation <- function(beta) {
  if (!is.matrix(beta) || nrow(beta) != ncol(beta))
    stop("`beta` must be a square matrix.", call. = FALSE)
  product <- beta * t(beta)
  product[product < 0] <- 0
  psi <- pmin(sqrt(product), 1) * sign(beta)
  dimnames(psi) <- dimnames(beta)
  diag(psi) <- 1
  psi
}

#' Plug-in precision matrix from a SUBHO fit
#'
#' Rescales the partial correlations with the node-wise residual variances:
#' \eqn{\hat\omega_{aa} = 1/\hat\sigma_a^2} and
#' \eqn{\hat\omega_{ab} = -\hat\psi_{ab}(\hat\omega_{aa}\hat\omega_{bb})^{1/2}}
#' (Peng et al. 2009), on the scale of the data used in the fit. Positive definiteness is
#' not guaranteed.
#'
#' @param fit A result of [subho()].
#' @param rule `"and"`, `"or"` or `"none"`: set entries outside this graph to zero, or
#'   keep all of them.
#' @return A symmetric matrix.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Peng J, Wang P, Zhou N, Zhu J (2009). Partial correlation estimation by joint sparse
#' regression models. Journal of the American Statistical Association, 104(486), 735-746.
#' @examples
#' sim <- simulate_ggm(n = 50, p = 6, structure = "ar1", seed = 1)
#' fit <- subho(sim$data, n_iter = 300, burn = 100, thin = 2, seed = 1)
#' round(precision_estimate(fit), 2)
#' @export
precision_estimate <- function(fit, rule = c("and", "or", "none")) {
  rule <- match.arg(rule)
  if (!inherits(fit, "ggm_fit") || is.null(fit$residual_variance))
    stop("`fit` must come from subho().", call. = FALSE)
  w <- 1 / fit$residual_variance
  psi <- fit$psi
  if (rule != "none") psi <- psi * ggm_graph(fit, rule)
  omega <- -psi * sqrt(outer(w, w))
  diag(omega) <- w
  dimnames(omega) <- dimnames(fit$psi)
  omega
}

#' Adjacency matrix of a fitted network
#'
#' @param fit A `ggm_fit` object.
#' @param rule `"and"` or `"or"`, for neighbourhood methods. Ignored by methods that give
#'   one graph.
#' @return A symmetric 0/1 matrix.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 6, structure = "ar1", seed = 1)
#' fit <- subho(sim$data, n_iter = 300, burn = 100, thin = 2, seed = 1)
#' ggm_graph(fit, "or")
#' @export
ggm_graph <- function(fit, rule = c("and", "or")) {
  rule <- match.arg(rule)
  if (!inherits(fit, "ggm_fit")) stop("`fit` must be a ggm_fit object.", call. = FALSE)
  if (!is.null(fit$graph)) return(fit$graph)
  if (rule == "and") fit$graph_and else fit$graph_or
}

#' @export
print.ggm_fit <- function(x, ...) {
  cat(sprintf("%s network: %d nodes, %d samples\n", x$method, x$p, x$n))
  ut <- upper.tri(x$psi)
  if (!is.null(x$graph)) {
    cat(sprintf("Edges: %d\n", sum(x$graph[ut] != 0)))
  } else {
    cat(sprintf("Edges: %d (AND rule), %d (OR rule)\n",
                sum(x$graph_and[ut] != 0), sum(x$graph_or[ut] != 0)))
  }
  if (!is.null(x$elapsed)) cat(sprintf("Elapsed: %.1f seconds\n", x$elapsed))
  invisible(x)
}

#' Summarise a fitted network
#'
#' @param object A `ggm_fit` object.
#' @param rule `"and"` or `"or"`, for neighbourhood methods.
#' @param ... Not used.
#' @return A one-row data frame with the number of nodes, samples and edges, the density,
#'   the counts of positive and negative edges, the median absolute partial correlation
#'   on edges, the isolated nodes and the largest degree.
#' @export
summary.ggm_fit <- function(object, rule = c("and", "or"), ...) {
  rule <- match.arg(rule)
  g <- ggm_graph(object, rule)
  ut <- upper.tri(g)
  e <- g[ut] != 0
  psi <- object$psi[ut]
  data.frame(
    method = object$method, nodes = object$p, samples = object$n, edges = sum(e),
    density = mean(e), positive = sum(psi[e] > 0), negative = sum(psi[e] < 0),
    median_abs_psi = if (any(e)) stats::median(abs(psi[e])) else NA_real_,
    isolated_nodes = sum(rowSums(g != 0) == 0), max_degree = max(rowSums(g != 0))
  )
}
