# README figure and numbers for the TNBC panel. Fits the first 60 genes as an illustration.
#   Rscript tools/make_readme_tnbc.R <library with GGMsubho> <output dir for figures>
args <- commandArgs(trailingOnly = TRUE)
library(GGMsubho, lib.loc = args[1])
library(ggplot2)
x <- tnbc_metabric[, 1:60]
fit <- subho(x, cores = 6, seed = 1)
print(fit)
print(summary(fit, "and"))
print(top_connected(fit, top = 5), digits = 3)
ggsave(file.path(args[2], "network_tnbc.png"),
       plot_network(fit, title = "METABRIC TNBC", subtitle = "SUBHO, AND rule, first 60 genes of the panel"),
       width = 7, height = 7.6, dpi = 150, bg = "white")
