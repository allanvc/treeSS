# Plot methods for scan results -----------------------------------------
#
# plot.treess() draws the detected cluster(s) of any spatial scan result
# (circular_scan, treespatial_scan, sequential_scan). It has two modes:
#
#   * map = NULL  -> base graphics: region centroids as points, cluster
#                    members filled by cluster. No extra dependencies.
#   * map = <sf>  -> ggplot2 + sf choropleth, one facet per cluster,
#                    joined to the scan's region_id through `map_id`
#                    (and, if the map's key differs from region_id,
#                    through `key`). Returns the ggplot object.
#
# Both modes are thin wrappers around get_cluster_regions(); users who
# prefer another toolkit (leaflet, tmap, ...) can call that function
# directly and join its output to their own layer.
#
# plot.tree_scan() has no geography; it draws a dot chart of the top
# tree cuts by log-likelihood ratio.

utils::globalVariables(".data")  # tidy-eval pronoun used inside aes()

#' Plot the Clusters Detected by a Scan
#'
#' Draws the most likely cluster and, optionally, the distinct secondary
#' clusters (single-pass scans) or every iteration (sequential scans) of
#' a spatial scan result. Without a polygon layer the method uses base
#' graphics and the region centroids stored in the result; with an
#' \pkg{sf} polygon layer it builds a faceted \pkg{ggplot2} choropleth.
#'
#' The method is a convenience wrapper around
#' \code{\link{get_cluster_regions}}: it performs the region-to-polygon
#' join, replicates the map once per cluster panel, and applies a fixed
#' palette. Users who prefer another mapping toolkit (e.g.
#' \pkg{leaflet} or \pkg{tmap}) can call \code{get_cluster_regions()}
#' directly and join its output to their own layer.
#'
#' @param x A \code{"circular_scan"}, \code{"treespatial_scan"}, or
#'   \code{"sequential_scan"} object (all inherit from \code{"treess"}).
#' @param map Optional \pkg{sf} data frame with one polygon per region.
#'   If \code{NULL} (default), the region centroids stored in the result
#'   are plotted with base graphics.
#' @param map_id Name of the column of \code{map} that identifies the
#'   region. Required when \code{map} is supplied.
#' @param key Optional two-column \code{data.frame} translating the
#'   scan's region identifiers into the values of \code{map[[map_id]]}:
#'   the first column must hold \code{region_id} values as used in the
#'   scan, the second the matching value in \code{map[[map_id]]}. Omit
#'   when \code{map[[map_id]]} already contains the scan's
#'   \code{region_id} values.
#' @param n_clusters Integer. Number of clusters to draw for single-pass
#'   scans: \code{1} (default) draws the most likely cluster only;
#'   larger values add the distinct secondary clusters retained by
#'   \code{\link{filter_clusters}}, one panel each. Ignored for
#'   \code{"sequential_scan"} objects, where every iteration is drawn.
#' @param palette Character vector of fill colours, recycled over
#'   clusters. The default is a four-colour qualitative palette.
#' @param na_fill Fill colour for regions outside the cluster of a panel.
#' @param wrap Integer. Width at which panel labels are wrapped.
#' @param ... For the base-graphics mode, further arguments passed to
#'   \code{\link[graphics]{plot}} (e.g. \code{cex}, \code{main}).
#'   Ignored in the \pkg{ggplot2} mode.
#'
#' @return In the base-graphics mode, the cluster-region table from
#'   \code{get_cluster_regions()} is returned invisibly. In the
#'   \pkg{ggplot2} mode a \code{ggplot} object is returned, which can be
#'   modified further with \code{+}.
#'
#' @seealso \code{\link{get_cluster_regions}},
#'   \code{\link{filter_clusters}}, \code{\link{sequential_scan}}
#'
#' @examples
#' data(london_collisions); data(london_tree)
#' res <- treespatial_scan(
#'   london_collisions,
#'   cases = cases, population = population, region_id = region_id,
#'   x = x, y = y, node_id = node_id, tree = london_tree,
#'   nsim = 99, seed = 42
#' )
#'
#' # Base graphics on centroids: most likely cluster
#' plot(res)
#'
#' # Two distinct clusters, one panel each
#' plot(res, n_clusters = 2)
#'
#' # Choropleth with the shipped borough polygons (requires sf, ggplot2).
#' # The scan uses integer region_id; the map is keyed by borough NAME,
#' # so a two-column key translates one into the other.
#' if (requireNamespace("sf", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(london_boroughs_map)
#'   key <- unique(london_collisions[, c("region_id", "borough")])
#'   plot(res, map = london_boroughs_map, map_id = "NAME", key = key)
#' }
#'
#' @export
plot.treess <- function(x, map = NULL, map_id = NULL, key = NULL,
                        n_clusters = 1L, palette = NULL,
                        na_fill = "gray95", wrap = 30L, ...) {

  if (is.null(x$regions)) {
    stop("This scan result has no region table; nothing to map. ",
         "For 'tree_scan' objects use plot() on the tree cuts instead.",
         call. = FALSE)
  }

  cr <- if (inherits(x, "sequential_scan")) {
    get_cluster_regions(x, overlap = TRUE)
  } else {
    get_cluster_regions(x, n_clusters = as.integer(n_clusters),
                        overlap = TRUE)
  }
  if (nrow(cr) == 0L) {
    message("No clusters to plot.")
    return(invisible(cr))
  }

  panels <- unique(cr$panel)
  cr$panel <- factor(cr$panel, levels = panels)
  n_pan <- length(panels)

  if (is.null(palette)) {
    palette <- c("#C44E52", "#4C72B0", "#55A868", "#8172B2")
  }
  pal <- rep_len(palette, n_pan)
  names(pal) <- as.character(seq_len(n_pan))

  wrap_label <- function(s) {
    vapply(strsplit(s, "\n", fixed = TRUE), function(parts) {
      paste(vapply(parts, function(p)
        paste(strwrap(p, width = wrap), collapse = "\n"),
        character(1L)), collapse = "\n")
    }, character(1L))
  }

  if (is.null(map)) {
    return(invisible(.plot_treess_base(cr, panels, pal, na_fill,
                                       wrap_label, ...)))
  }

  # --- ggplot2 + sf mode -------------------------------------------------
  if (!requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'sf' and 'ggplot2' are required to plot on a polygon ",
         "layer; install them or call plot() without 'map'.",
         call. = FALSE)
  }
  if (!inherits(map, "sf")) {
    stop("'map' must be an sf data frame.", call. = FALSE)
  }
  if (is.null(map_id) || !map_id %in% names(map)) {
    stop("'map_id' must name a column of 'map'.", call. = FALSE)
  }

  # Translate scan region_id into the map's key, if they differ.
  if (!is.null(key)) {
    key <- as.data.frame(key)
    if (ncol(key) < 2L) {
      stop("'key' must have two columns: region_id and the matching ",
           "value of map[[map_id]].", call. = FALSE)
    }
    cr$.map_key <- key[[2L]][match(cr$region_id, key[[1L]])]
  } else {
    cr$.map_key <- cr$region_id
  }

  members <- cr[!is.na(cr$cluster), ]
  unmatched <- setdiff(unique(members$.map_key), map[[map_id]])
  if (length(unmatched) > 0L) {
    warning(length(unmatched), " cluster region(s) have no polygon in ",
            "'map' (e.g. ", paste(utils::head(unmatched, 3L),
                                  collapse = ", "), ").",
            call. = FALSE)
  }

  pieces <- lapply(panels, function(p) {
    d <- members[members$panel == p, ]
    m <- map
    m$panel   <- p
    m$cluster <- d$cluster[match(m[[map_id]], d$.map_key)]
    m
  })
  mm <- do.call(rbind, pieces)
  mm$panel   <- factor(mm$panel, levels = panels)
  mm$cluster <- factor(mm$cluster, levels = names(pal))

  ggplot2::ggplot(mm) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[["cluster"]]),
                     color = "gray40", linewidth = 0.12, alpha = 0.80,
                     show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, na.value = na_fill,
                               na.translate = FALSE) +
    ggplot2::facet_wrap(~ panel, nrow = 1L,
                        labeller = ggplot2::labeller(panel = wrap_label)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(color = "gray55", size = 8)
    )
}


# Base-graphics fallback: centroids, one panel per cluster.
.plot_treess_base <- function(cr, panels, pal, na_fill, wrap_label, ...) {
  n_pan <- length(panels)
  titles <- wrap_label(as.character(panels))
  n_lines <- max(lengths(strsplit(titles, "\n", fixed = TRUE)))
  op <- graphics::par(mfrow = c(1L, n_pan),
                      mar = c(2.5, 2.5, 1.5 + 1.1 * n_lines, 0.8),
                      mgp = c(1.5, 0.4, 0))
  on.exit(graphics::par(op), add = TRUE)

  for (i in seq_len(n_pan)) {
    d <- cr[cr$panel == panels[i], ]
    in_cl <- !is.na(d$cluster)
    graphics::plot(d$x, d$y, type = "n", asp = 1, xlab = "", ylab = "",
                   main = titles[i], cex.main = 0.9,
                   cex.axis = 0.7, ...)
    graphics::points(d$x[!in_cl], d$y[!in_cl], pch = 21,
                     bg = na_fill, col = "gray50", cex = 1)
    graphics::points(d$x[in_cl], d$y[in_cl], pch = 21,
                     bg = pal[i], col = "gray30", cex = 1.4)
  }
  cr
}


#' Plot the Top Tree Cuts of a Tree-Only Scan
#'
#' A \code{"tree_scan"} result has no geography. Its \code{plot} method
#' draws a dot chart of the highest-scoring tree cuts by log-likelihood
#' ratio, marking those significant at the scan's \code{alpha}.
#'
#' @param x A \code{"tree_scan"} object.
#' @param n Integer. Number of top cuts to display (default 10).
#' @param ... Further arguments passed to \code{\link[graphics]{dotchart}}.
#'
#' @return The table of plotted cuts, invisibly.
#'
#' @seealso \code{\link{tree_scan}}, \code{\link{plot.treess}}
#'
#' @examples
#' tree <- data.frame(node_id   = c(1, 2, 3, 4, 5, 6, 7, 8),
#'                    parent_id = c(NA, 1, 1, 2, 2, 3, 3, 3))
#' leaf_data <- data.frame(node_id = c(4, 5, 6, 7, 8),
#'                         cases   = c(50, 5, 3, 2, 4))
#' res <- tree_scan(leaf_data, cases = cases, node_id = node_id,
#'                  tree = tree, nsim = 99, seed = 1)
#' plot(res, n = 8)
#'
#' @export
plot.tree_scan <- function(x, n = 10L, ...) {
  cuts <- utils::head(x$all_cuts[order(-x$all_cuts$llr), ], n)
  sig  <- cuts$pvalue < x$alpha
  cols <- ifelse(sig, "#C44E52", "gray40")
  graphics::dotchart(rev(cuts$llr),
                     labels = rev(as.character(cuts$node_id)),
                     pch = 19, color = rev(cols),
                     xlab = "Log-likelihood ratio",
                     main = sprintf("Top %d tree cuts (red: p < %s)",
                                    nrow(cuts), format(x$alpha)),
                     ...)
  invisible(cuts)
}
