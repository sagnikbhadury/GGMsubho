# Internal helpers shared by the estimators.
# Method reference for the package: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

.check_data <- function(data, standardize = TRUE, min_p = 2L) {
  if (is.data.frame(data)) {
    if (!all(vapply(data, is.numeric, logical(1))))
      stop("All columns of `data` must be numeric.", call. = FALSE)
    data <- as.matrix(data)
  }
  if (!is.matrix(data) || !is.numeric(data))
    stop("`data` must be a numeric matrix or data frame with samples in rows.", call. = FALSE)
  if (anyNA(data) || any(!is.finite(data)))
    stop("`data` contains missing or non-finite values; handle them before fitting.", call. = FALSE)
  if (ncol(data) < min_p)
    stop(sprintf("`data` must have at least %d columns.", min_p), call. = FALSE)
  if (nrow(data) < 3L)
    stop("`data` must have at least 3 rows.", call. = FALSE)
  sds <- apply(data, 2, stats::sd)
  if (any(sds == 0))
    stop("Constant column(s) in `data`: ",
         paste(which(sds == 0), collapse = ", "), call. = FALSE)
  if (is.null(colnames(data))) colnames(data) <- paste0("V", seq_len(ncol(data)))
  if (anyDuplicated(colnames(data)))
    stop("Column names of `data` must be unique.", call. = FALSE)
  nodes <- colnames(data)
  if (standardize) data <- scale(data)
  attr(data, "scaled:center") <- NULL
  attr(data, "scaled:scale") <- NULL
  rownames(data) <- NULL
  colnames(data) <- nodes
  data
}

.check_mcmc <- function(n_iter, burn, thin) {
  ok <- function(x) length(x) == 1L && is.numeric(x) && is.finite(x) && x == round(x)
  if (!ok(n_iter) || !ok(burn) || !ok(thin) || burn < 0 || thin < 1 || n_iter <= burn)
    stop("`n_iter`, `burn` and `thin` must be whole numbers with 0 <= burn < n_iter and thin >= 1.",
         call. = FALSE)
  if ((n_iter - burn) %/% thin < 2L)
    stop("`n_iter`, `burn` and `thin` leave fewer than 2 saved draws.", call. = FALSE)
  invisible(TRUE)
}

.resolve_seed <- function(seed) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max %/% 2L, 1L))
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed))
    stop("`seed` must be a single number or NULL.", call. = FALSE)
  as.integer(seed)
}

.check_cores <- function(cores) {
  if (length(cores) != 1L || !is.numeric(cores) || cores < 1 || cores != round(cores))
    stop("`cores` must be a positive whole number.", call. = FALSE)
  as.integer(cores)
}

# Runs one seeded task. The seed depends only on the node, so a fit is identical for
# any number of cores, and the caller's random number stream is left untouched.
.seeded_task <- function(a, node_fun, payload, seed) {
  withr::with_seed(seed + a, node_fun(a, payload),
                   .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                   .rng_sample_kind = "Rejection")
}

.run_nodes <- function(p, fun, payload, seed, cores = 1L, verbose = FALSE) {
  if (cores <= 1L) {
    return(lapply(seq_len(p), function(a) {
      if (verbose && (a %% max(1L, p %/% 10L) == 0L))
        message(sprintf("node %d of %d", a, p))
      .seeded_task(a, fun, payload, seed)
    }))
  }
  cl <- parallel::makeCluster(min(cores, p))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # Workers run the installed package's code; one BLAS thread is assumed per worker.
  lib_paths <- unique(c(dirname(getNamespaceInfo("GGMsubho", "path")), .libPaths()))
  # Sent as a call: .libPaths() keeps its state in its closure, which would be copied.
  parallel::clusterCall(cl, eval, call(".libPaths", lib_paths))
  parallel::clusterCall(cl, loadNamespace, "GGMsubho")
  parallel::parLapply(cl, seq_len(p), .seeded_task, node_fun = fun, payload = payload, seed = seed)
}

.empty_matrix <- function(nodes) {
  matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
}

.require_suggested <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("%s needs the '%s' package. Install it with install.packages(\"%s\").",
                 what, pkg, pkg), call. = FALSE)
  invisible(TRUE)
}

.new_fit <- function(method, nodes, n, fields, call, elapsed, seed = NULL) {
  structure(c(list(method = method, nodes = nodes, n = n, p = length(nodes)), fields,
              list(call = call, elapsed = elapsed, seed = seed)),
            class = "ggm_fit")
}

# Psi from the precision matrix: psi_ab = -omega_ab / sqrt(omega_aa omega_bb)
.psi_from_omega <- function(omega) {
  d <- sqrt(diag(omega))
  psi <- -omega / outer(d, d)
  diag(psi) <- 1
  dimnames(psi) <- dimnames(omega)
  psi
}
