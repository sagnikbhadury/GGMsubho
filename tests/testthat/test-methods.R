sim <- simulate_ggm(n = 100, p = 6, structure = "ar1", seed = 21)

test_that("simulate_ggm gives a consistent truth", {
  expect_equal(unname(diag(sim$sigma)), rep(1, 6))
  expect_equal(sum(sim$graph) / 2, 5)
  expect_equal(sim$partial_cor[1, 2], 0.7 / sqrt(1.49), tolerance = 1e-8)
  # beta and omega agree: beta[a, b] = -omega[a, b] / omega[a, a]
  o <- solve(sim$sigma)
  expect_equal(sim$true_beta[2, 3], -o[2, 3] / o[2, 2], tolerance = 1e-8)
  r <- simulate_ggm(n = 20, p = 30, structure = "random", edge_prob = 0.1, seed = 1)
  expect_true(all(eigen(r$sigma, only.values = TRUE)$values > 0))
  expect_equal(which(r$partial_cor[upper.tri(r$partial_cor)] != 0),
               which(r$graph[upper.tri(r$graph)] != 0))
})

test_that("graph_metrics counts correctly and handles an empty graph", {
  truth <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
  est <- rbind(c(0, 1, 1), c(1, 0, 0), c(1, 0, 0))
  m <- graph_metrics(est, truth)
  expect_equal(m[c("tp", "fp", "fn", "tn")], c(tp = 1, fp = 1, fn = 1, tn = 0))
  expect_equal(m[["fdr"]], 0.5)
  empty <- graph_metrics(matrix(0, 3, 3), truth)
  expect_equal(empty[["fdr"]], 0)
  expect_equal(empty[["mcc"]], 0)
})

test_that("the other metrics agree with hand calculation", {
  expect_equal(partial_cor_error(c(0.1, 0.5, 0), c(0, 0.4, 0))[["error_total"]], 0.02)
  t1 <- top_edge_metrics(c(0.9, 0.05, 0.4), c(1, 0, 0), top = 1)
  expect_equal(t1[["tp"]], 1)
  expect_equal(t1[["threshold"]], 0.9)
  s <- correlation_summary(sim$partial_cor, sim$sigma)
  expect_equal(sum(s$partial$count), 15)
  expect_equal(s$marginal$count[s$marginal$bin == "0"], 10)
})

test_that("the comparison methods run and find strong edges", {
  fits <- list(
    basad = basad_ggm(sim$data, n_iter = 500, burn = 200, thin = 2, seed = 1),
    basad_f = basad_ggm(sim$data, selection = "fdr", n_iter = 500, burn = 200, thin = 2, seed = 1),
    ghs = ghs_ggm(sim$data, burn = 100, n_save = 200, seed = 1),
    iw = iw_ggm(sim$data, n_draws = 300, seed = 1),
    lasso = lasso_ggm(sim$data, seed = 1)
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_s3_class(f, "ggm_fit")
    expect_gte(graph_metrics(ggm_graph(f, "or"), sim$graph)[["tpr"]], 0.6)
    expect_equal(sign(f$psi[1, 2]), 1)
  }
  expect_output(print(fits$ghs), "Edges")
})

test_that("BASAD handles more predictors than samples", {
  s <- simulate_ggm(n = 15, p = 25, structure = "ar1", seed = 3)
  f <- basad_ggm(s$data, n_iter = 60, burn = 20, thin = 2, seed = 1)
  expect_equal(dim(f$gamma), c(25L, 25L))
  expect_true(basad_start_q(150, 150) > 0 && basad_start_q(150, 150) < 1)
  expect_equal(pbinom(10, 150, basad_start_q(150, 150)), 0.9, tolerance = 0.01)
})

test_that("network tools return consistent objects", {
  net <- network_data(sim$partial_cor)
  expect_equal(nrow(net$edges), 5L)
  expect_true(all(net$edges$partial_correlation > 0))
  expect_equal(igraph::ecount(net$graph), 5L)
  # a precision matrix is converted, with the sign reversed
  o <- solve(sim$sigma)
  o[abs(o) < 1e-10] <- 0
  npr <- network_data(o, type = "precision")
  expect_equal(npr$edges$partial_correlation, net$edges$partial_correlation, tolerance = 1e-6)
  expect_s3_class(plot_network(sim$partial_cor), "ggplot")
  expect_s3_class(plot_network(sim$partial_cor, layout = "stress", curvature = 0), "ggplot")
  cl <- cluster_nodes(sim$data, k = 2)
  expect_equal(length(cl), 6L)
  expect_s3_class(plot_network(sim$partial_cor, clusters = cl), "ggplot")
  expect_s3_class(plot_network(matrix(0, 3, 3)), "ggplot")
  expect_s3_class(plot_matrix(sim$partial_cor), "ggplot")
  expect_error(network_data(matrix(1:4, 2)), "symmetric")
  built <- ggplot2::ggplot_build(plot_network(sim$partial_cor, clusters = cl))
  expect_true(length(built$data) >= 2L)
})

test_that("connectivity ranking breaks ties by strength and overlap counts", {
  f <- subho(sim$data, n_iter = 200, burn = 50, thin = 2, seed = 8)
  top <- top_connected(f)
  expect_equal(nrow(top), 6L)
  expect_true(all(diff(top$degree) <= 0))
  ov <- top_overlap(list(a = f, b = f), top = 3)
  expect_equal(unname(ov), matrix(3L, 2, 2))
  expect_s3_class(plot_top_overlap(ov), "ggplot")
  expect_true(summary(f)$edges == sum(f$graph_and) / 2)
  expect_equal(dim(precision_estimate(f)), c(6L, 6L))
  expect_true(isSymmetric(precision_estimate(f)))
})

test_that("results do not depend on the number of cores", {
  skip_on_cran()
  skip_if_not_installed("GGMsubho")
  f1 <- subho(sim$data, n_iter = 200, burn = 50, thin = 2, seed = 5, cores = 1)
  f2 <- subho(sim$data, n_iter = 200, burn = 50, thin = 2, seed = 5, cores = 2)
  expect_identical(f1$beta, f2$beta)
  expect_identical(f1$gamma, f2$gamma)
})

test_that("optional packages produce a fit when they are installed", {
  skip_on_cran()
  if (requireNamespace("CVglasso", quietly = TRUE)) {
    expect_s3_class(glasso_ggm(sim$data, seed = 1), "ggm_fit")
  }
  if (requireNamespace("BAS", quietly = TRUE)) {
    f <- bas_ggm(sim$data, seed = 1)
    expect_gte(graph_metrics(f$graph_or, sim$graph)[["tpr"]], 0.6)
  }
  if (requireNamespace("BDgraph", quietly = TRUE)) {
    f <- bdgraph_ggm(sim$data, iter = 800, burnin = 300, seed = 1)
    expect_s3_class(f, "ggm_fit")
    g <- simulate_ggm(n = 30, p = 8, structure = "gwishart", edge_prob = 0.2, seed = 1)
    expect_equal(unname(diag(g$sigma)), rep(1, 8))
  }
  if (requireNamespace("pROC", quietly = TRUE)) {
    expect_equal(edge_auc(c(0.9, 0.1, 0.4, 0.2), c(1, 0, 1, 0)), 1)
  }
})
