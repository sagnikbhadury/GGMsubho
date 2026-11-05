# Comparison methods of the SUBHO article.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

#' Graphical horseshoe
#'
#' Gibbs sampler of Li, Craig and Bhadra (2019), with a horseshoe prior on the
#' off-diagonal entries of the precision matrix. An edge is kept when the central
#' `level` credible interval of the entry excludes zero.
#'
#' @inheritParams subho
#' @param burn Burn-in iterations.
#' @param n_save Iterations saved after burn-in.
#' @param level Credible level for the graph.
#' @return A `ggm_fit` object with `graph`, `psi` and `omega` (posterior mean).
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Li Y, Craig BA, Bhadra A (2019). The graphical horseshoe estimator for inverse
#' covariance matrices. Journal of Computational and Graphical Statistics, 28(3), 747-757.
#' @examples
#' sim <- simulate_ggm(n = 50, p = 6, structure = "ar1", seed = 1)
#' fit <- ghs_ggm(sim$data, burn = 50, n_save = 100, seed = 1)
#' fit
#' @export
ghs_ggm <- function(data, level = 0.5, standardize = TRUE, burn = 500, n_save = 2000,
                    seed = NULL) {
  cl <- match.call()
  seed <- .resolve_seed(seed)
  X <- .check_data(data, standardize, min_p = 3L)
  nodes <- colnames(X)
  t0 <- proc.time()[["elapsed"]]
  draws <- withr::with_seed(seed, .ghs_gibbs(X, burn, n_save))
  elapsed <- proc.time()[["elapsed"]] - t0
  omega <- apply(draws, c(1L, 2L), mean)
  alpha <- 1 - level
  lo <- apply(draws, c(1L, 2L), stats::quantile, probs = alpha / 2)
  hi <- apply(draws, c(1L, 2L), stats::quantile, probs = 1 - alpha / 2)
  graph <- 1 * (lo * hi > 0)
  diag(graph) <- 0
  dimnames(graph) <- dimnames(omega) <- list(nodes, nodes)
  .new_fit("GHS", nodes, nrow(X),
           list(graph = graph, psi = .psi_from_omega(omega), omega = omega), cl, elapsed, seed)
}

.ghs_gibbs <- function(X, burn, n_save) {
  n <- nrow(X)
  S <- crossprod(X)
  p <- nrow(S)
  out <- array(0, c(p, p, n_save))
  Omega <- Sigma <- diag(p)
  Lambda2 <- Nu <- matrix(1, p, p)
  tau2 <- 1
  xi <- 1
  rinv <- function(k, shape, rate) 1 / stats::rgamma(k, shape, rate = rate)
  for (iter in seq_len(burn + n_save)) {
    for (i in seq_len(p)) {
      ind <- seq_len(p)[-i]
      Sigma11 <- Sigma[ind, ind, drop = FALSE]
      sigma12 <- Sigma[ind, i]
      s21 <- S[ind, i]
      s22 <- S[i, i]
      gamma <- stats::rgamma(1L, n / 2 + 1, rate = s22 / 2)
      inv_omega11 <- Sigma11 - tcrossprod(sigma12) / Sigma[i, i]
      inv_c <- s22 * inv_omega11 + diag(1 / (Lambda2[ind, i] * tau2), p - 1L)
      R <- chol(inv_c)
      mu <- -chol2inv(R) %*% s21
      beta <- as.vector(mu + backsolve(R, stats::rnorm(p - 1L)))
      lambda2 <- rinv(p - 1L, 1, 1 / Nu[ind, i] + beta^2 / (2 * tau2))
      nu <- rinv(p - 1L, 1, 1 + 1 / lambda2)
      Omega[i, ind] <- Omega[ind, i] <- beta
      temp <- as.vector(inv_omega11 %*% beta)
      Omega[i, i] <- gamma + sum(beta * temp)
      Sigma[ind, ind] <- inv_omega11 + tcrossprod(temp) / gamma
      Sigma[i, ind] <- Sigma[ind, i] <- -temp / gamma
      Sigma[i, i] <- 1 / gamma
      Lambda2[i, ind] <- Lambda2[ind, i] <- lambda2
      Nu[i, ind] <- Nu[ind, i] <- nu
    }
    w <- Omega[lower.tri(Omega)]
    l <- Lambda2[lower.tri(Lambda2)]
    tau2 <- rinv(1L, (length(w) + 1) / 2, 1 / xi + sum(w^2 / (2 * l)))
    xi <- rinv(1L, 1, 1 + 1 / tau2)
    if (iter > burn) out[, , iter - burn] <- Omega
    if (iter %% 50L == 0L) Sys.sleep(0)
  }
  out
}

#' Inverse Wishart estimate of a network
#'
#' With the conjugate prior \eqn{\Sigma \sim IW(m, I)}, \eqn{m = p + 2}, the precision
#' matrix has posterior mean \eqn{(n + m)(I + X^TX)^{-1}}. An edge is kept when the
#' central `level` credible interval of the entry, computed from independent posterior
#' draws, excludes zero. This estimator is not sparse and is included for comparison.
#'
#' @inheritParams subho
#' @param level Credible level for the graph.
#' @param n_draws Number of posterior draws.
#' @return A `ggm_fit` object with `graph`, `psi` and `omega`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 6, structure = "ar1", seed = 1)
#' iw_ggm(sim$data, n_draws = 200, seed = 1)
#' @export
iw_ggm <- function(data, level = 0.5, n_draws = 2000, standardize = TRUE, seed = NULL) {
  cl <- match.call()
  seed <- .resolve_seed(seed)
  X <- .check_data(data, standardize)
  nodes <- colnames(X)
  p <- ncol(X)
  t0 <- proc.time()[["elapsed"]]
  df <- p + 2 + nrow(X)
  inv_sn <- solve(diag(p) + crossprod(X))
  draws <- withr::with_seed(seed, stats::rWishart(n_draws, df, inv_sn))
  alpha <- 1 - level
  lo <- apply(draws, c(1L, 2L), stats::quantile, probs = alpha / 2)
  hi <- apply(draws, c(1L, 2L), stats::quantile, probs = 1 - alpha / 2)
  graph <- 1 * (lo * hi > 0)
  diag(graph) <- 0
  omega <- df * inv_sn
  dimnames(graph) <- dimnames(omega) <- list(nodes, nodes)
  .new_fit("Inverse Wishart", nodes, nrow(X),
           list(graph = graph, psi = .psi_from_omega(omega), omega = omega),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

#' Neighbourhood selection with the lasso
#'
#' The approach of Meinshausen and Buhlmann (2006), with the penalty chosen by
#' cross-validation at each node ([glmnet::cv.glmnet()]).
#'
#' @inheritParams subho
#' @param s Penalty passed to the `coef` method of `cv.glmnet`: `"lambda.1se"` or
#'   `"lambda.min"`.
#' @return A `ggm_fit` object with `graph_and`, `graph_or`, `psi`, `gamma` and `beta`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Meinshausen N, Buhlmann P (2006). High-dimensional graphs and variable selection with
#' the lasso. The Annals of Statistics, 34(3), 1436-1462.
#' @examples
#' sim <- simulate_ggm(n = 60, p = 6, structure = "ar1", seed = 1)
#' lasso_ggm(sim$data, seed = 1)
#' @export
lasso_ggm <- function(data, s = c("lambda.1se", "lambda.min"), standardize = TRUE, cores = 1,
                      seed = NULL) {
  cl <- match.call()
  s <- match.arg(s)
  cores <- .check_cores(cores)
  seed <- .resolve_seed(seed)
  Z <- .check_data(data, standardize, min_p = 3L)
  nodes <- colnames(Z)
  t0 <- proc.time()[["elapsed"]]
  res <- .run_nodes(ncol(Z), .lasso_node, list(Z = Z, s = s), seed, cores)
  beta <- .empty_matrix(nodes)
  for (a in seq_along(nodes)) beta[a, -a] <- res[[a]]
  .new_fit("LASSO", nodes, nrow(Z),
           list(graph_and = symmetrize_graph(beta, "and"), graph_or = symmetrize_graph(beta, "or"),
                psi = partial_correlation(beta), gamma = 1 * (beta != 0), beta = beta),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

.lasso_node <- function(a, payload) {
  fit <- glmnet::cv.glmnet(payload$Z[, -a, drop = FALSE], payload$Z[, a], family = "gaussian",
                           alpha = 1)
  as.numeric(stats::coef(fit, s = payload$s))[-1L]
}

#' Graphical lasso with cross-validation
#'
#' Fits the graphical lasso (Friedman, Hastie and Tibshirani 2008) through
#' `CVglasso::CVglasso()`. The graph is the support of the estimated precision matrix.
#' Requires the suggested package 'CVglasso'.
#'
#' @inheritParams subho
#' @param nlam,lam.min.ratio,K_folds,cv_cores Number of penalties, smallest to largest
#'   penalty ratio, cross-validation folds and cores.
#' @return A `ggm_fit` object with `graph`, `psi` and `omega`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse covariance estimation with
#' the graphical lasso. Biostatistics, 9(3), 432-441.
#' @examplesIf requireNamespace("CVglasso", quietly = TRUE)
#' sim <- simulate_ggm(n = 60, p = 6, structure = "ar1", seed = 1)
#' glasso_ggm(sim$data, seed = 1)
#' @export
glasso_ggm <- function(data, nlam = 10, lam.min.ratio = 0.01, K_folds = 5, cv_cores = 1,
                       standardize = TRUE, seed = NULL) {
  .require_suggested("CVglasso", "glasso_ggm()")
  cl <- match.call()
  seed <- .resolve_seed(seed)
  X <- .check_data(data, standardize)
  nodes <- colnames(X)
  t0 <- proc.time()[["elapsed"]]
  res <- withr::with_seed(seed, CVglasso::CVglasso(
    X = X, nlam = nlam, lam.min.ratio = lam.min.ratio, diagonal = FALSE, path = FALSE,
    tol = 1e-4, maxit = 10000, K = K_folds, crit.cv = "loglik", start = "warm",
    cores = cv_cores, trace = "none"))
  omega <- res$Omega
  dimnames(omega) <- list(nodes, nodes)
  graph <- 1 * (omega != 0 & t(omega) != 0)
  diag(graph) <- 0
  .new_fit("GLASSO", nodes, nrow(X),
           list(graph = graph, psi = .psi_from_omega(omega), omega = omega, lambda = res$Tuning[2L]),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

#' G-Wishart network with BDgraph
#'
#' Fits the birth-death MCMC of `BDgraph::bdgraph()` with an edge prior of 0.5 and three
#' prior degrees of freedom. The graph is the median graph, the edges with posterior
#' inclusion probability of at least 0.5. Requires the suggested package 'BDgraph'.
#'
#' @inheritParams subho
#' @param iter,burnin Iterations and burn-in passed to `bdgraph()`.
#' @return A `ggm_fit` object with `graph`, `psi`, `omega` and `edge_probability`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Mohammadi R, Wit EC (2019). BDgraph: an R package for Bayesian structure learning in
#' graphical models. Journal of Statistical Software, 89(3), 1-30.
#' @examplesIf requireNamespace("BDgraph", quietly = TRUE)
#' \donttest{
#' sim <- simulate_ggm(n = 60, p = 6, structure = "ar1", seed = 1)
#' bdgraph_ggm(sim$data, iter = 1000, burnin = 500, seed = 1)
#' }
#' @export
bdgraph_ggm <- function(data, iter = 20000, burnin = 8000, standardize = TRUE, seed = NULL) {
  .require_suggested("BDgraph", "bdgraph_ggm()")
  cl <- match.call()
  seed <- .resolve_seed(seed)
  X <- .check_data(data, standardize)
  nodes <- colnames(X)
  t0 <- proc.time()[["elapsed"]]
  res <- withr::with_seed(seed, BDgraph::bdgraph(
    X, method = "ggm", algorithm = "bdmcmc", iter = iter, burnin = burnin, g.prior = 0.5,
    df.prior = 3, g.start = "empty", save = TRUE, cores = 1, verbose = FALSE))
  prob <- as.matrix(BDgraph::plinks(res))
  prob[lower.tri(prob)] <- t(prob)[lower.tri(prob)]
  omega <- res$K_hat
  dimnames(prob) <- dimnames(omega) <- list(nodes, nodes)
  graph <- 1 * (prob >= 0.5)
  diag(graph) <- 0
  .new_fit("BDgraph", nodes, nrow(X),
           list(graph = graph, psi = .psi_from_omega(omega), omega = omega,
                edge_probability = prob),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

#' Neighbourhood selection with a hyper-g prior
#'
#' Runs `BAS::bas.lm()` at each node with `prior = "hyper-g"` and `alpha = 3`, sampling at
#' most \eqn{2^{19}} models. The neighbourhood is the median probability model. Requires
#' the suggested package 'BAS'.
#'
#' @inheritParams subho
#' @param n_models Largest number of models sampled at a node. `NULL` means
#'   `2^min(p, 19)`, as in the authors' code.
#' @return A `ggm_fit` object with `graph_and`, `graph_or`, `psi`, `gamma` and `beta`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Liang F, Paulo R, Molina G, Clyde MA, Berger JO (2008). Mixtures of g priors for
#' Bayesian variable selection. Journal of the American Statistical Association,
#' 103(481), 410-423.
#' @examplesIf requireNamespace("BAS", quietly = TRUE)
#' \donttest{
#' sim <- simulate_ggm(n = 60, p = 6, structure = "ar1", seed = 1)
#' bas_ggm(sim$data, seed = 1)
#' }
#' @export
bas_ggm <- function(data, n_models = NULL, standardize = TRUE, cores = 1, seed = NULL) {
  .require_suggested("BAS", "bas_ggm()")
  cl <- match.call()
  cores <- .check_cores(cores)
  seed <- .resolve_seed(seed)
  Z <- .check_data(data, standardize, min_p = 3L)
  nodes <- colnames(Z)
  if (is.null(n_models)) n_models <- 2^min(ncol(Z), 19)
  t0 <- proc.time()[["elapsed"]]
  res <- .run_nodes(ncol(Z), .bas_node, list(Z = Z, n_models = n_models), seed, cores)
  beta <- .empty_matrix(nodes)
  for (a in seq_along(nodes)) beta[a, -a] <- res[[a]]
  .new_fit("Hyper-g (BAS)", nodes, nrow(Z),
           list(graph_and = symmetrize_graph(beta, "and"), graph_or = symmetrize_graph(beta, "or"),
                psi = partial_correlation(beta), gamma = 1 * (beta != 0), beta = beta),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

.bas_node <- function(a, payload) {
  Z <- payload$Z
  pred <- paste0("x", seq_len(ncol(Z) - 1L))
  df <- as.data.frame(Z[, -a, drop = FALSE])
  names(df) <- pred
  df$y <- Z[, a]
  fit <- do.call(BAS::bas.lm, list(formula = stats::as.formula("y ~ ."), data = df, alpha = 3,
                                   n.models = payload$n_models, prior = "hyper-g"))
  cf <- stats::coef(fit, estimator = "MPM")
  out <- stats::setNames(numeric(length(pred)), pred)
  nm <- cf$namesx
  keep <- nm %in% pred
  out[nm[keep]] <- cf$postmean[keep]
  unname(out)
}

#' Neighbourhood selection with projection predictive inference
#'
#' The approach of Williams et al. (2018): a horseshoe reference model is fitted with
#' `rstanarm::stan_glm()`, and `projpred::cv_varsel()` selects a submodel by forward
#' search with approximate leave-one-out cross-validation. This is very slow beyond a few
#' dozen variables. Requires the suggested packages 'rstanarm' and 'projpred'.
#'
#' @inheritParams subho
#' @param iter,chains,adapt_delta Settings of `stan_glm()`.
#' @param alpha Passed to `projpred::suggest_size()`.
#' @param fallback_size Submodel size used when `suggest_size()` returns `NA`.
#' @return A `ggm_fit` object with `graph_and`, `graph_or`, `psi`, `gamma` and `beta`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Williams DR, Piironen J, Vehtari A, Rast P (2018). Bayesian estimation of Gaussian
#' graphical models with predictive covariance selection. arXiv:1801.05725.
#'
#' Piironen J, Paasiniemi M, Vehtari A (2020). Projective inference in high-dimensional
#' problems: prediction and feature selection. Electronic Journal of Statistics, 14(1),
#' 2155-2197.
#' @examplesIf interactive() && requireNamespace("projpred", quietly = TRUE) && requireNamespace("rstanarm", quietly = TRUE)
#' sim <- simulate_ggm(n = 60, p = 4, structure = "ar1", seed = 1)
#' projpred_ggm(sim$data, iter = 500, seed = 1)
#' @export
projpred_ggm <- function(data, iter = 1500, chains = 2, adapt_delta = 0.95, alpha = 0.1,
                         fallback_size = 19, standardize = TRUE, cores = 1, seed = NULL) {
  .require_suggested("rstanarm", "projpred_ggm()")
  .require_suggested("projpred", "projpred_ggm()")
  cl <- match.call()
  cores <- .check_cores(cores)
  seed <- .resolve_seed(seed)
  Z <- .check_data(data, standardize, min_p = 3L)
  nodes <- colnames(Z)
  t0 <- proc.time()[["elapsed"]]
  payload <- list(Z = Z, iter = iter, chains = chains, adapt_delta = adapt_delta, alpha = alpha,
                  fallback_size = fallback_size)
  res <- .run_nodes(ncol(Z), .projpred_node, payload, seed, cores)
  beta <- .empty_matrix(nodes)
  for (a in seq_along(nodes)) beta[a, -a] <- res[[a]]
  .new_fit("Projpred-HS", nodes, nrow(Z),
           list(graph_and = symmetrize_graph(beta, "and"), graph_or = symmetrize_graph(beta, "or"),
                psi = partial_correlation(beta), gamma = 1 * (beta != 0), beta = beta),
           cl, proc.time()[["elapsed"]] - t0, seed)
}

.projpred_node <- function(a, payload) {
  Z <- payload$Z
  pred <- paste0("x", seq_len(ncol(Z) - 1L))
  df <- as.data.frame(Z[, -a, drop = FALSE])
  names(df) <- pred
  df$y <- Z[, a]
  fit <- suppressWarnings(rstanarm::stan_glm(
    stats::as.formula("y ~ ."), data = df, family = stats::gaussian(),
    prior = rstanarm::hs(1, 1), adapt_delta = payload$adapt_delta, chains = payload$chains,
    iter = payload$iter, refresh = 0))
  vs <- suppressWarnings(projpred::cv_varsel(fit, method = "forward", cv_method = "LOO",
                                             verbose = FALSE))
  s <- suppressWarnings(projpred::suggest_size(vs, alpha = payload$alpha))
  if (is.na(s)) s <- min(payload$fallback_size, length(pred))
  out <- stats::setNames(numeric(length(pred)), pred)
  if (s > 0) {
    terms <- projpred::ranking(vs)$fulldata[seq_len(s)]
    proj <- suppressWarnings(projpred::project(vs, nterms = s, refit_prj = TRUE, verbose = FALSE))
    draws <- as.matrix(proj)
    out[terms] <- colMeans(draws[, terms, drop = FALSE])
  }
  unname(out)
}
