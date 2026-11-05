# Builds data/tnbc_metabric.rda from the processed 1000-gene panel.
#   Rscript tools/make_tnbc_data.R <path to Low1000pval_TNBCpositive_data.rds> <package dir>
# The sample identifiers are dropped; the expression values are unchanged.
args <- commandArgs(trailingOnly = TRUE)
x <- as.matrix(readRDS(args[1]))
stopifnot(identical(dim(x), c(299L, 1000L)), !anyNA(x), !anyDuplicated(colnames(x)))
tnbc_metabric <- x
rownames(tnbc_metabric) <- NULL
dir.create(file.path(args[2], "data"), showWarnings = FALSE)
save(tnbc_metabric, file = file.path(args[2], "data", "tnbc_metabric.rda"),
     compress = "xz", compression_level = 9, version = 2)
cat("KB:", round(file.size(file.path(args[2], "data", "tnbc_metabric.rda")) / 1024), "\n")
