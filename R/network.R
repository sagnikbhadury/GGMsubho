# Network construction and figures.
# Method reference: Bhadury, Mitra and Gaskins (2023), arXiv:2308.02713.

# Partial correlation matrix and graph from a fit or a matrix.
.network_matrix <- function(x, rule, type) {
  if (inherits(x, "ggm_fit")) {
    return(list(psi = x$psi * ggm_graph(x, rule), method = x$method))
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x))
    stop("`x` must be a ggm_fit object or a square numeric matrix.", call. = FALSE)
  if (!isTRUE(all.equal(unname(x), unname(t(x)), tolerance = 1e-8)))
    stop("The matrix `x` must be symmetric.", call. = FALSE)
  if (is.null(colnames(x))) dimnames(x) <- list(paste0("V", seq_len(ncol(x))), paste0("V", seq_len(ncol(x))))
  if (is.null(rownames(x))) rownames(x) <- colnames(x)
  psi <- if (type == "precision") .psi_from_omega(x) else x
  list(psi = psi, method = NULL)
}

#' Edges and nodes of a fitted network
#'
#' Turns a fit or a partial correlation matrix into the tables and the `igraph` object
#' behind [plot_network()].
#'
#' @param x A `ggm_fit` object, or a symmetric matrix of partial correlations or a
#'   precision matrix (see `type`).
#' @param rule `"and"` or `"or"`, for neighbourhood methods.
#' @param threshold Edges with a smaller absolute partial correlation are dropped. An edge
#'   whose two coefficients disagree in sign has a partial correlation of exactly zero, so
#'   the default removes it.
#' @param type For a matrix `x`: `"partial_correlation"` or `"precision"`. A precision
#'   matrix is converted with \eqn{\psi_{ab} = -\omega_{ab}/(\omega_{aa}\omega_{bb})^{1/2}}.
#' @param keep_isolated Keep nodes with no edges?
#'
#' @return A list with `edges` (`from`, `to`, `partial_correlation`), `nodes` (`name`,
#'   `degree`, `strength`, the sum of absolute partial correlations) and `graph`, an
#'   `igraph` object.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 6, structure = "ar1", seed = 1)
#' net <- network_data(sim$partial_cor)
#' net$edges
#' @export
network_data <- function(x, rule = c("and", "or"), threshold = 1e-8,
                         type = c("partial_correlation", "precision"), keep_isolated = FALSE) {
  rule <- match.arg(rule)
  type <- match.arg(type)
  psi <- .network_matrix(x, rule, type)$psi
  psi[row(psi) == col(psi)] <- 0
  keep <- upper.tri(psi) & abs(psi) >= threshold
  idx <- which(keep, arr.ind = TRUE)
  nm <- colnames(psi)
  edges <- data.frame(from = nm[idx[, 1L]], to = nm[idx[, 2L]],
                      partial_correlation = psi[idx], stringsAsFactors = FALSE)
  active <- abs(psi) >= threshold
  nodes <- data.frame(name = nm, degree = rowSums(active), strength = rowSums(abs(psi) * active),
                      stringsAsFactors = FALSE, row.names = NULL)
  if (!keep_isolated) nodes <- nodes[nodes$degree > 0, , drop = FALSE]
  list(edges = edges, nodes = nodes,
       graph = igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes))
}

#' Draw a network
#'
#' Draws a network with nodes on a circle, by default with curved edges coloured by the
#' sign of the partial correlation and scaled by its magnitude, and with nodes sized by the
#' sum of absolute partial correlations. Isolated nodes are left out unless
#' `keep_isolated = TRUE`. For hundreds of nodes, use `show_labels = FALSE`.
#'
#' Edges can instead be coloured by groups of nodes, as in the gene network figures of
#' Bhadury et al. (2023), by giving `clusters` (for example from [cluster_nodes()]): edges
#' inside a cluster take its colour, and edges between clusters share one colour.
#'
#' @inheritParams network_data
#' @param clusters Optional vector of node clusters, either named by node or in node order.
#' @param layout Layout passed to `ggraph::ggraph()`; `"circle"` by default. Others include
#'   `"stress"` and `"fr"`.
#' @param title,subtitle Text of the figure. By default the title is the method's name.
#' @param show_labels Draw node labels?
#' @param label_size Size of the labels.
#' @param node_colour Colour of the nodes.
#' @param edge_alpha Transparency of the edges.
#' @param edge_width_range,node_size_range Ranges of edge widths and node sizes.
#' @param curvature Curvature of the edges; 0 draws straight lines.
#' @param base_family Font family; `""` uses the default of the graphics device.
#'
#' @return A `ggplot` object, so it can be changed with `+` and saved with
#'   `ggplot2::ggsave()`. The data behind it are available from [network_data()].
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 12, structure = "random", edge_prob = 0.15, seed = 1)
#' plot_network(sim$partial_cor, title = "True network")
#' @export
plot_network <- function(x, rule = c("and", "or"), threshold = 1e-8,
                         type = c("partial_correlation", "precision"), clusters = NULL,
                         keep_isolated = FALSE, layout = "circle", title = NULL, subtitle = NULL,
                         show_labels = TRUE, label_size = 2.6, node_colour = "blue",
                         edge_alpha = 0.75, edge_width_range = c(0.2, 1.4),
                         node_size_range = c(0.8, 5), curvature = 0.2, base_family = "") {
  rule <- match.arg(rule)
  type <- match.arg(type)
  net <- network_data(x, rule, threshold, type, keep_isolated)
  if (is.null(title) && inherits(x, "ggm_fit")) title <- paste(x$method, "network")
  if (nrow(net$edges) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = title, subtitle = paste(c(subtitle, "no edges"), collapse = ": ")) +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                            plot.subtitle = ggplot2::element_text(hjust = 0.5)))
  }
  g <- net$graph
  igraph::E(g)$magnitude <- abs(net$edges$partial_correlation)
  igraph::E(g)$sign_label <- factor(ifelse(net$edges$partial_correlation > 0, "Positive", "Negative"),
                                    levels = c("Positive", "Negative"))
  has_clusters <- !is.null(clusters)
  if (has_clusters) {
    igraph::E(g)$edge_group <- .edge_groups(clusters, net, .network_nodes(x))
  }
  is_circle <- identical(layout, "circle")
  p <- ggraph::ggraph(g, layout = layout)
  if (has_clusters) {
    aes_edge <- ggplot2::aes(colour = edge_group, width = magnitude)
  } else {
    aes_edge <- ggplot2::aes(colour = sign_label, width = magnitude)
  }
  if (curvature != 0) {
    p <- p + ggraph::geom_edge_arc(aes_edge, strength = curvature, alpha = edge_alpha)
  } else {
    p <- p + ggraph::geom_edge_link(aes_edge, alpha = edge_alpha)
  }
  p <- p +
    ggraph::scale_edge_width(name = "Magnitude", range = edge_width_range) +
    ggraph::geom_node_point(ggplot2::aes(size = strength), colour = node_colour) +
    ggplot2::scale_size(name = "Connectivity (sum of |partial correlation|)",
                        range = node_size_range)
  p <- p + if (has_clusters) {
    ggraph::scale_edge_colour_manual(name = "Edge", values = .cluster_palette(levels(igraph::E(g)$edge_group)))
  } else {
    ggraph::scale_edge_colour_manual(name = "Sign", values = c(Positive = "#D55E00", Negative = "#0072B2"),
                                     drop = FALSE)
  }
  if (show_labels) {
    if (is_circle) {
      p <- p + ggraph::geom_node_text(
        ggplot2::aes(label = name, x = x * 1.15, y = y * 1.15,
                     angle = ifelse(atan(-(x / y)) * (180 / pi) < 0,
                                    90 + atan(-(x / y)) * (180 / pi),
                                    270 + atan(-x / y) * (180 / pi)),
                     hjust = ifelse(x > 0, 0, 1)),
        size = label_size, fontface = "bold", family = base_family)
    } else {
      p <- p + ggraph::geom_node_text(ggplot2::aes(label = name), size = label_size, repel = TRUE,
                                      fontface = "bold", family = base_family)
    }
  }
  if (is_circle) p <- p + ggplot2::expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))
  p + ggplot2::coord_fixed() + ggplot2::theme_void(base_family = base_family) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 8, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.position = "bottom", legend.box = "vertical") +
    ggplot2::labs(title = title, subtitle = subtitle)
}

.network_nodes <- function(x) if (inherits(x, "ggm_fit")) x$nodes else colnames(x)

.edge_groups <- function(clusters, net, all_nodes) {
  if (is.null(names(clusters))) {
    if (is.null(all_nodes) || length(clusters) != length(all_nodes))
      stop("An unnamed `clusters` must have one entry per node.", call. = FALSE)
    names(clusters) <- all_nodes
  }
  if (!all(c(net$edges$from, net$edges$to) %in% names(clusters)))
    stop("`clusters` does not cover every node in the network.", call. = FALSE)
  a <- as.character(clusters[net$edges$from])
  b <- as.character(clusters[net$edges$to])
  lev <- paste("Cluster", sort(unique(as.vector(clusters))))
  factor(ifelse(a == b, paste("Cluster", a), "Between clusters"), levels = c(lev, "Between clusters"))
}

# The cluster colours of the published gene network figures, extended when needed.
.cluster_palette <- function(levels) {
  base <- c("red", "green", "blue", "purple", "orange", "brown", "magenta", "gold", "darkgreen", "grey40")
  within <- setdiff(levels, "Between clusters")
  cols <- rep_len(base, length(within))
  stats::setNames(c(cols, "aquamarine"), c(within, "Between clusters"))
}

#' Draw several networks side by side
#'
#' Draws each fit with [plot_network()] in a grid, without labels by default, as used to
#' compare methods in Bhadury et al. (2023). Requires the suggested package 'patchwork'.
#'
#' @param fits A named list of `ggm_fit` objects or matrices. The names become titles.
#' @param ncol Number of columns.
#' @param show_labels Draw node labels?
#' @param legend Keep legends?
#' @param ... Other arguments passed to [plot_network()].
#' @return A `patchwork` object.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examplesIf requireNamespace("patchwork", quietly = TRUE)
#' sim <- simulate_ggm(n = 50, p = 12, structure = "random", edge_prob = 0.15, seed = 1)
#' plot_network_panels(list(Truth = sim$partial_cor, Again = sim$partial_cor))
#' @export
plot_network_panels <- function(fits, ncol = NULL, show_labels = FALSE, legend = FALSE, ...) {
  .require_suggested("patchwork", "plot_network_panels()")
  if (!is.list(fits) || inherits(fits, "ggm_fit") || is.null(names(fits)))
    stop("`fits` must be a named list.", call. = FALSE)
  plots <- lapply(names(fits), function(nm) {
    p <- plot_network(fits[[nm]], title = nm, show_labels = show_labels, ...)
    if (legend) p else p + ggplot2::theme(legend.position = "none")
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

#' Most connected nodes
#'
#' Ranks the nodes by their number of connections, \eqn{G_a = \sum_b \hat G_{ab}}, and
#' breaks ties by the sum of absolute partial correlations,
#' \eqn{P_a = \sum_b |\hat\psi_{ab}|}.
#'
#' @param fit A `ggm_fit` object.
#' @param rule `"and"` or `"or"`, for neighbourhood methods.
#' @param top Number of nodes returned; `NULL` returns all.
#' @return A data frame with `node`, `degree` and `strength`, best first.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 8, structure = "ar1", seed = 1)
#' fit <- subho(sim$data, n_iter = 300, burn = 100, thin = 2, seed = 1)
#' top_connected(fit, top = 3)
#' @export
top_connected <- function(fit, rule = c("and", "or"), top = NULL) {
  rule <- match.arg(rule)
  g <- ggm_graph(fit, rule)
  psi <- fit$psi
  psi[row(psi) == col(psi)] <- 0
  out <- data.frame(node = fit$nodes, degree = as.numeric(rowSums(g != 0)),
                    strength = as.numeric(rowSums(abs(psi))), stringsAsFactors = FALSE)
  out <- out[order(-out$degree, -out$strength), , drop = FALSE]
  rownames(out) <- NULL
  if (is.null(top)) out else utils::head(out, top)
}

#' Overlap of the most connected nodes across methods
#'
#' Counts, for every pair of fits, the nodes shared by their [top_connected()] lists.
#'
#' @param fits A named list of `ggm_fit` objects.
#' @param top Number of most connected nodes per fit.
#' @param rule `"and"` or `"or"`, for neighbourhood methods.
#' @return A symmetric integer matrix.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 8, structure = "ar1", seed = 1)
#' fits <- list(LASSO = lasso_ggm(sim$data, seed = 1),
#'              IW = iw_ggm(sim$data, n_draws = 200, seed = 1))
#' top_overlap(fits, top = 3)
#' @export
top_overlap <- function(fits, top = 50, rule = c("and", "or")) {
  rule <- match.arg(rule)
  sets <- lapply(fits, function(f) top_connected(f, rule, top)$node)
  m <- length(sets)
  out <- matrix(0L, m, m, dimnames = list(names(fits), names(fits)))
  for (i in seq_len(m)) for (j in seq_len(m)) out[i, j] <- length(intersect(sets[[i]], sets[[j]]))
  out
}

#' Draw the overlap of the most connected nodes
#'
#' @param overlap Matrix from [top_overlap()].
#' @param title Figure title.
#' @return A `ggplot` object.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' m <- matrix(c(10L, 4L, 4L, 10L), 2, dimnames = list(c("A", "B"), c("A", "B")))
#' plot_top_overlap(m)
#' @export
plot_top_overlap <- function(overlap, title = NULL) {
  nm <- rownames(overlap)
  df <- data.frame(method_x = factor(rep(nm, times = ncol(overlap)), levels = nm),
                   method_y = factor(rep(colnames(overlap), each = nrow(overlap)),
                                     levels = rev(colnames(overlap))),
                   overlap = as.vector(overlap))
  ggplot2::ggplot(df, ggplot2::aes(method_x, method_y, fill = overlap)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = overlap), colour = "white", size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#9ecae1", high = "#08306b", name = "Shared nodes") +
    ggplot2::coord_fixed() + ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Draw an adjacency or partial correlation matrix
#'
#' @param x A `ggm_fit` object or a square matrix.
#' @param rule `"and"` or `"or"`, for neighbourhood methods.
#' @param title Figure title.
#' @return A `ggplot` object.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 12, structure = "random", edge_prob = 0.15, seed = 1)
#' plot_matrix(sim$partial_cor, title = "True partial correlations")
#' @export
plot_matrix <- function(x, rule = c("and", "or"), title = NULL) {
  rule <- match.arg(rule)
  m <- if (inherits(x, "ggm_fit")) x$psi * ggm_graph(x, rule) else x
  m[row(m) == col(m)] <- NA
  nm <- if (is.null(colnames(m))) as.character(seq_len(ncol(m))) else colnames(m)
  df <- data.frame(row_node = factor(rep(nm, times = length(nm)), levels = rev(nm)),
                   col_node = factor(rep(nm, each = length(nm)), levels = nm),
                   value = as.vector(m))
  ggplot2::ggplot(df, ggplot2::aes(col_node, row_node, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                                  limits = c(-1, 1), na.value = "grey85",
                                  name = "Partial\ncorrelation") +
    ggplot2::coord_fixed() + ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = if (length(nm) > 30) ggplot2::element_blank() else ggplot2::element_text(size = 6),
                   axis.text.x = if (length(nm) > 30) ggplot2::element_blank() else
                     ggplot2::element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
}

#' Cluster the nodes of a data set
#'
#' Hierarchical clustering of the variables on Euclidean distance with complete linkage,
#' cut into `k` groups. Useful for colouring edges in [plot_network()].
#'
#' @param data Numeric matrix with samples in rows and variables in columns.
#' @param k Number of clusters.
#' @param method Agglomeration method passed to [stats::hclust()].
#' @return An integer vector named by variable.
#' @references
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian
#' Graphical Model Estimation. arXiv:2308.02713.
#' \url{https://arxiv.org/abs/2308.02713}
#' @examples
#' sim <- simulate_ggm(n = 50, p = 12, structure = "ar1", seed = 1)
#' cluster_nodes(sim$data, k = 3)
#' @export
cluster_nodes <- function(data, k = 4, method = "complete") {
  Z <- .check_data(data, standardize = FALSE)
  stats::cutree(stats::hclust(stats::dist(t(Z), method = "euclidean"), method = method), k = k)
}
