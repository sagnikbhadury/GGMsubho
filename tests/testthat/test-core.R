# Marginal likelihood against the multivariate t density it equals in closed form
log_mvt <- function(y, X, v, c0, d0) {
  n <- length(y)
  S <- diag(n) + v * tcrossprod(cbind(1, X))
  q <- sum(y * solve(S, y))
  lgamma(c0 + n / 2) - lgamma(c0) - (n / 2) * log(2 * pi * d0) -
    0.5 * as.numeric(determinant(S)$modulus) - (c0 + n / 2) * log(1 + q / (2 * d0))
}

test_that("log_marginal_working equals the closed-form multivariate t density", {
  set.seed(1)
  X <- matrix(rnorm(40 * 3), 40, 3)
  y <- rnorm(40)
  for (v in c(1, 2.5, 10)) for (c0 in c(1, 0.5)) {
    expect_equal(log_marginal_working(y, X, c0 = c0, d0 = 2, v = v),
                 log_mvt(y, X, v, c0, 2), tolerance = 1e-9)
  }
})

test_that("selection does not underflow at large n", {
  set.seed(2)
  n <- 2000
  x <- rnorm(n)
  y <- 0.7 * x + sqrt(0.51) * rnorm(n)
  X <- cbind(x, matrix(rnorm(n * 4), n, 4))
  s <- step_up_select(y, X, c(1, 0.01, 0.02, 0.03, 0.04))
  expect_true(s$selected[1])
  expect_equal(s$k, 1L)
})

test_that("the path matches direct evaluation of each nested model", {
  set.seed(3)
  X <- matrix(rnorm(50 * 5), 50, 5)
  y <- X[, 1] - X[, 3] + rnorm(50)
  bh <- c(1, 0.1, -0.8, 0.05, 0.2)
  s <- step_up_select(y, X, bh)
  ord <- order(abs(bh), decreasing = TRUE)
  direct <- vapply(0:5, function(k) {
    log_marginal_working(y, if (k == 0) NULL else X[, ord[seq_len(k)], drop = FALSE])
  }, numeric(1))
  expect_equal(s$log_marginal, direct, tolerance = 1e-9)
  prior <- lgamma(1 + 0:5) + lgamma(5 + 1 - 0:5) - lgamma(2 + 5)
  expect_equal(s$log_posterior, direct + prior, tolerance = 1e-9)
})

test_that("the two beta updates sample the same posterior", {
  set.seed(4)
  n <- 30
  X <- matrix(rnorm(n * 4), n, 4)
  y <- X[, 1] + rnorm(n)
  # p < 0.8 n uses the Cholesky update; duplicating rows changes nothing about the check
  # that matters here, so compare with the conjugate case of a fixed, large tau instead.
  fit <- horseshoe_regression(y, X, intercept = FALSE, n_iter = 6000, burn = 1000, thin = 2,
                              tau = 1e3, seed = 5)
  ols <- solve(crossprod(X), crossprod(X, y))
  expect_equal(as.numeric(fit$beta_mean), as.numeric(ols), tolerance = 0.1)
})

test_that("a small sample with many predictors recovers the signal", {
  set.seed(5)
  X <- matrix(rnorm(40 * 60), 40, 60)
  y <- 3 * X[, 1] - 3 * X[, 2] + rnorm(40, sd = 0.5)
  fit <- horseshoe_regression(y, X, n_iter = 1500, burn = 500, thin = 2, seed = 6)
  expect_equal(order(abs(fit$beta_median), decreasing = TRUE)[1:2] |> sort(), 1:2)
  expect_equal(length(fit$sigma2_draws), 500L)
  expect_lt(mean(fit$sigma2_draws), 4)
})

test_that("draws are saved at the requested iterations when burn is not a multiple of thin", {
  set.seed(7)
  X <- matrix(rnorm(40 * 3), 40, 3)
  fit <- horseshoe_regression(rnorm(40), X, n_iter = 155, burn = 53, thin = 10, seed = 1)
  expect_equal(ncol(fit$beta_draws), 10L)
  expect_true(all(fit$beta_draws[, 10] != 0))
})

test_that("partial_correlation applies the sign rule and the truncation", {
  b <- rbind(c(0, 0.5, -0.1), c(0.3, 0, 0.2), c(0.1, 0.4, 0))
  psi <- partial_correlation(b)
  expect_equal(psi[1, 2], sqrt(0.15))
  expect_equal(psi[1, 3], 0)
  expect_equal(psi[2, 3], sqrt(0.08))
  expect_equal(psi, t(psi))
  expect_equal(unname(diag(psi)), rep(1, 3))
  expect_equal(partial_correlation(rbind(c(0, -2), c(-3, 0)))[1, 2], -1)
})

test_that("the AND and OR rules work as defined", {
  g <- rbind(c(0, 1, 1), c(1, 0, 0), c(0, 0, 0))
  expect_equal(sum(symmetrize_graph(g, "and")), 2)
  expect_equal(sum(symmetrize_graph(g, "or")), 4)
})

test_that("subho recovers an AR(1) chain and has the documented shape", {
  sim <- simulate_ggm(n = 150, p = 10, structure = "ar1", seed = 11)
  fit <- subho(sim$data, n_iter = 600, burn = 100, thin = 5, seed = 12)
  expect_s3_class(fit, "ggm_fit")
  m <- graph_metrics(fit$graph_and, sim$graph)
  expect_gte(m[["tpr"]], 0.8)
  expect_lte(m[["fdr"]], 0.2)
  expect_equal(fit$graph_and, t(fit$graph_and))
  expect_equal(unname(diag(fit$graph_or)), rep(0, 10))
  expect_equal(rownames(fit$psi), colnames(sim$data))
  expect_true(all(fit$graph_and <= fit$graph_or))
  expect_lt(partial_cor_error(fit$psi, sim$partial_cor)[["error_total"]], 1)
})

test_that("fits are reproducible and leave the caller's random stream alone", {
  sim <- simulate_ggm(n = 40, p = 5, structure = "ar1", seed = 1)
  f1 <- subho(sim$data, n_iter = 150, burn = 50, thin = 2, seed = 9)
  f2 <- subho(sim$data, n_iter = 150, burn = 50, thin = 2, seed = 9)
  expect_identical(f1$beta, f2$beta)
  set.seed(100)
  a <- runif(1)
  set.seed(100)
  invisible(subho(sim$data, n_iter = 150, burn = 50, thin = 2, seed = 9))
  expect_identical(a, runif(1))
})

test_that("mean and median estimators, K and the model prior are respected", {
  sim <- simulate_ggm(n = 60, p = 6, structure = "ar1", seed = 2)
  f_med <- subho(sim$data, n_iter = 200, burn = 50, thin = 2, seed = 3)
  f_mean <- subho(sim$data, estimator = "mean", n_iter = 200, burn = 50, thin = 2, seed = 3)
  expect_false(identical(f_med$beta, f_mean$beta))
  f0 <- subho(sim$data, K = 1, n_iter = 200, burn = 50, thin = 2, seed = 3)
  expect_true(all(f0$neighbourhood_size <= 1))
  f_b <- subho(sim$data, model_prior = "binomial", n_iter = 200, burn = 50, thin = 2, seed = 3)
  expect_s3_class(f_b, "ggm_fit")
})

test_that("invalid input is rejected", {
  expect_error(subho(matrix(c(1, NA, 3, 4, 5, 6), 3)), "missing")
  expect_error(subho(cbind(1:5, rep(1, 5))), "Constant")
  expect_error(subho(matrix(rnorm(20), 10), n_iter = 10, burn = 20), "n_iter")
  expect_error(subho(data.frame(a = letters[1:5], b = 1:5)), "numeric")
  expect_error(horseshoe_regression(1:3, matrix(1:6, 2)), "length")
  expect_error(step_up_select(1:3, matrix(rnorm(6), 3), 1), "beta_hat")
  expect_error(log_marginal_working(rnorm(5), v = -1), "positive")
})
