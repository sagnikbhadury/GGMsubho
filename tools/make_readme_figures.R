# Builds every figure and number in README.md. All data are simulated.
#   Rscript tools/make_readme_figures.R <library with GGMsubho> <output dir for figures>
args <- commandArgs(trailingOnly = TRUE)
library(GGMsubho, lib.loc = args[1])
library(ggplot2)
out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
CORES <- 6
save_fig <- function(p, file, w, h) ggsave(file.path(out, file), p, width = w, height = h, dpi = 150, bg = "white")

sim <- simulate_ggm(n = 120, p = 40, structure = "random", edge_prob = 0.06,
                    node_names = sprintf("Gene%02d", 1:40), seed = 2026)
cat("true edges:", sum(sim$graph[upper.tri(sim$graph)]), "\n")
tp <- sim$partial_cor[upper.tri(sim$partial_cor)]; tp <- tp[tp != 0]
cat("true |psi| range:", round(range(abs(tp)), 3), "\n")

fit <- subho(sim$data, cores = CORES, seed = 1)
print(fit)
print(summary(fit, "and"))
print(round(graph_metrics(fit$graph_and, sim$graph), 3))
print(round(graph_metrics(fit$graph_or, sim$graph), 3))
print(head(network_data(fit)$edges, 5), digits = 3)
print(top_connected(fit, top = 5), digits = 3)
print(round(partial_cor_error(fit$psi * fit$graph_and, sim$partial_cor), 3))

save_fig(plot_network(fit, subtitle = "AND rule"), "network_subho.png", 7, 7.6)
save_fig(plot_network(sim$partial_cor, title = "True network"), "network_truth.png", 7, 7.6)
cl <- cluster_nodes(sim$data, k = 4)
save_fig(plot_network(fit, clusters = cl, subtitle = "Edges coloured by cluster of nodes"), "network_clusters.png", 7, 7.6)
save_fig(plot_network(fit, layout = "stress", curvature = 0, subtitle = "Force-directed layout"), "network_stress.png", 7, 7.6)

t0 <- proc.time()[3]
fits <- list(SUBHO = fit,
             `Graphical horseshoe` = ghs_ggm(sim$data, seed = 1),
             LASSO = lasso_ggm(sim$data, cores = CORES, seed = 1))
cat("competitors secs:", round(proc.time()[3] - t0, 1), "\n")
tab <- t(sapply(fits, function(f) {
  g <- ggm_graph(f, "and")
  m <- graph_metrics(g, sim$graph)
  c(m[c("edges", "fdr", "tpr", "mcc")], seconds = f$elapsed)
}))
print(round(tab, 3))
panel <- c(list(Truth = sim$partial_cor), fits)
save_fig(plot_network_panels(panel, ncol = 2), "network_panels.png", 9, 9)
save_fig(plot_top_overlap(top_overlap(fits, top = 10), "Shared nodes among the 10 most connected"), "top_overlap.png", 6, 5)
save_fig(plot_matrix(fit, title = "SUBHO, AND rule") , "matrix.png", 6, 5)

big <- simulate_ggm(n = 150, p = 300, structure = "random", edge_prob = 0.01, seed = 7)
fb <- subho(big$data, cores = CORES, seed = 2)
cat("large network secs:", round(fb$elapsed, 1), "on", CORES, "cores\n")
print(round(graph_metrics(fb$graph_and, big$graph), 3))
save_fig(plot_network(fb, show_labels = FALSE, node_size_range = c(0.3, 2.5),
                      edge_width_range = c(0.1, 0.6),
                      subtitle = "n = 150, p = 300, simulated"), "network_large.png", 7, 7.6)
file.size(list.files(out, full.names = TRUE))
