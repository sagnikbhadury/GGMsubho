# Neighbourhood selection with BASAD priors, a comparison method in the SUBHO article.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713;
# prior: Narisetty and He (2014), The Annals of Statistics 42(2), 789-817.

#' Starting inclusion probability for BASAD
#'
#' Finds \eqn{q} with \eqn{P(\mathrm{Binomial}(p, q) > k) = 0.1}, for
#' \eqn{k = \max(10, \log n)}, by damped Newton steps on the squared difference.
#'
#' @param n,p Numbers of samples and predictors.
#' @param eps Convergence tolerance.
#' @param start Starting value.
#' @param max_iter Largest number of steps.
#' @return The value of \eqn{q}. When `p <= k` the condition holds for every \eqn{q} and
#'   `start` is returned.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Narisetty NN, He X (2014). Bayesian variable selection with shrinking and diffusing
#' priors. The Annals of Statistics, 42(2), 789-817.
#' @examples
#' basad_start_q(150, 150)
#' @export
basad_start_q <- function(n, p, eps = 1e-5, start = 0.05, max_iter = 1e5) {
  k <- max(10, log(n))
  if (p <= k) return(start)
  g <- function(q) (stats::pbinom(k, p, q) - 0.9)^2
  g_prime <- function(q) {
    i <- 0:floor(k)
    2 * sum(stats::dbinom(i, p, q) * (i - p * q) / (q * (1 - q))) *
      (sum(stats::dbinom(i, p, q)) - 0.9)
  }
  q <- start
  for (count in seq_len(max_iter)) {
    q_next <- q - 0.1 * g(q) / g_prime(q)
    if (!is.finite(q_next)) return(q)
    if (q_next < 0) q_next <- count / 1000
    if (q_next > 1) q_next <- 1 - count / 1000
    if (abs(q_next - q) < eps) return(q)
    q <- q_next
  }
  warning("basad_start_q() did not converge.", call. = FALSE)
  q
}

# Gibbs sampler for one regression; the design already contains an intercept column.
.basad_regression <- function(X, y, n_iter, burn, thin, start_q, fdr) {
  n <- nrow(X)
  p <- ncol(X)
  n_save <- (n_iter - burn) %/% thin
  sigma2 <- 1
  tau0 <- sigma2 / (10 * n)
  tau1 <- max(p^2.1 / (100 * n), log(n)) * sigma2
  z <- stats::rbinom(p, 1L, start_q)
  beta_out <- z_out <- zp_out <- matrix(0, p, n_save)
  primal <- p < 0.8 * n
  if (primal) {
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)
  }
  for (i in seq_len(n_iter)) {
    d <- z * tau1 + (1 - z) * tau0
    if (primal) {
      A <- XtX
      diag(A) <- diag(A) + 1 / d
      U <- chol(A)
      beta <- backsolve(U, backsolve(U, Xty, transpose = TRUE) +
                          sqrt(sigma2) * stats::rnorm(p))
    } else {
      s <- sqrt(sigma2)
      u <- sqrt(sigma2 * d) * stats::rnorm(p)
      v <- X %*% u / s + stats::rnorm(n)
      XD <- sweep(X, 2L, sqrt(d), "*")
      M <- tcrossprod(XD)
      diag(M) <- diag(M) + 1
      w <- solve(M, y / s - v)
      beta <- u + (sigma2 * d) * crossprod(X, w) / s
    }
    q <- stats::rbeta(1L, 0.5 + sum(z), 0.5 + p - sum(z))
    s1 <- q * stats::dnorm(beta, 0, sqrt(sigma2 * tau1))
    s0 <- (1 - q) * stats::dnorm(beta, 0, sqrt(sigma2 * tau0))
    zp <- as.numeric(s1 / (s1 + s0))
    zp[!is.finite(zp)] <- 1
    z <- stats::rbinom(p, 1L, zp)
    d <- z * tau1 + (1 - z) * tau0
    resid <- y - X %*% beta
    sigma2 <- 1 / stats::rgamma(1L, (n + p + 0.2) / 2,
                                rate = (sum(resid^2) + sum(beta^2 / d) + 0.2) / 2)
    if (i > burn && (i - burn) %% thin == 0L) {
      k <- (i - burn) %/% thin
      if (k <= n_save) {
        beta_out[, k] <- beta
        z_out[, k] <- z
        zp_out[, k] <- zp
      }
    }
  }
  prob <- rowMeans(zp_out)
  sorted <- sort(prob, decreasing = TRUE)
  cross <- which(cumsum(1 - sorted) / seq_along(sorted) > fdr)
  z_fdr <- if (length(cross)) as.numeric(prob > sorted[min(cross)]) else rep(1, p)
  list(beta = rowMeans(beta_out), z_median = as.numeric(apply(z_out, 1L, stats::median) > 0),
       z_fdr = z_fdr)
}

#' Neighbourhood selection with BASAD priors
#'
#' Spike-and-slab regression at each node with the shrinking and diffusing prior of
#' Narisetty and He (2014). The spike and slab variances are
#' \eqn{\tau_0^2 = 1/(10n)} and \eqn{\tau_1^2 = \max(p^{2.1}/(100n), \log n)}, and the
#' inclusion probability has a Beta(0.5, 0.5) prior updated within the sampler, starting
#' at [basad_start_q()]. Selection uses the median probability model or a Bayesian FDR
#' rule on the posterior inclusion probabilities.
#'
#' @inheritParams subho
#' @param selection `"median"` for the median probability model or `"fdr"` for the
#'   Bayesian FDR rule.
#' @param fdr Target of the Bayesian FDR rule.
#' @return A `ggm_fit` object with `graph_and`, `graph_or`, `psi`, `gamma` and `beta`.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#'
#' Narisetty NN, He X (2014). Bayesian variable selection with shrinking and diffusing
#' priors. The Annals of Statistics, 42(2), 789-817.
#' @examples
#' sim <- simulate_ggm(n = 50, p = 5, structure = "ar1", seed = 1)
#' fit <- basad_ggm(sim$data, n_iter = 300, burn = 100, thin = 2, seed = 1)
#' fit
#' @export
basad_ggm <- function(data, selection = c("median", "fdr"), fdr = 0.1, standardize = TRUE,
                      n_iter = 7000, burn = 3000, thin = 10, cores = 1, seed = NULL,
                      verbose = FALSE) {
  cl <- match.call()
  selection <- match.arg(selection)
  .check_mcmc(n_iter, burn, thin)
  cores <- .check_cores(cores)
  seed <- .resolve_seed(seed)
  Z <- .check_data(data, standardize)
  nodes <- colnames(Z)
  p <- ncol(Z)
  payload <- list(Z = Z, n_iter = n_iter, burn = burn, thin = thin, fdr = fdr,
                  start_q = basad_start_q(nrow(Z), p))
  t0 <- proc.time()[["elapsed"]]
  res <- .run_nodes(p, .basad_node, payload, seed, cores, verbose)
  elapsed <- proc.time()[["elapsed"]] - t0
  beta <- gamma <- .empty_matrix(nodes)
  field <- if (selection == "median") "z_median" else "z_fdr"
  for (a in seq_len(p)) {
    beta[a, -a] <- res[[a]]$beta
    gamma[a, -a] <- res[[a]][[field]]
  }
  .new_fit("BASAD", nodes, nrow(Z),
           list(graph_and = symmetrize_graph(gamma, "and"), graph_or = symmetrize_graph(gamma, "or"),
                psi = partial_correlation(beta), gamma = gamma, beta = beta),
           cl, elapsed, seed)
}

.basad_node <- function(a, payload) {
  Z <- payload$Z
  r <- .basad_regression(cbind(1, Z[, -a, drop = FALSE]), Z[, a], payload$n_iter,
                         payload$burn, payload$thin, payload$start_q, payload$fdr)
  list(beta = r$beta[-1L], z_median = r$z_median[-1L], z_fdr = r$z_fdr[-1L])
}
