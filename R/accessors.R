# Accessor generics for scan results ------------------------------------
#
# All scan functions in treeSS return objects that inherit from the
# common parent class "treess":
#
#   circular_scan()    -> c("circular_scan",    "treess")
#   tree_scan()        -> c("tree_scan",        "treess")
#   treespatial_scan() -> c("treespatial_scan", "treess")
#   sequential_scan()  -> c("sequential_scan",  "treess")
#
# The accessors below dispatch on "treess" whenever the underlying
# element is stored the same way across subclasses, and on the specific
# subclass when the representation differs (e.g. sequential_scan stores
# one cluster per iteration instead of a single most likely cluster).
# User code should rely on these accessors instead of reaching into the
# internal list structure (e.g. prefer most_likely_cluster(res) over
# res$most_likely_cluster).

#' Extract the Most Likely Cluster from a Scan Result
#'
#' Returns the most likely cluster detected by a scan, without requiring
#' the user to access the internal list structure of the result object.
#'
#' @param x A scan result object: \code{"circular_scan"},
#'   \code{"tree_scan"}, \code{"treespatial_scan"}, or
#'   \code{"sequential_scan"} (all of which inherit from
#'   \code{"treess"}).
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A named \code{list} describing the most likely cluster. The
#'   components depend on the scan type: for \code{circular_scan},
#'   region IDs, cases, expected count, population, relative risk, and
#'   log-likelihood ratio; for \code{tree_scan}, additionally the node
#'   ID and leaf IDs; for \code{treespatial_scan}, both the tree node
#'   (with leaf IDs) and the spatial zone (with region IDs). For
#'   \code{sequential_scan}, the cluster detected in the first
#'   iteration is returned (subsequent iterations are secondary
#'   clusters by construction; see
#'   \code{\link{secondary_clusters}}).
#'
#' @seealso \code{\link{secondary_clusters}}, \code{\link{pvalue}},
#'   \code{\link{get_cluster_regions}}
#'
#' @examples
#' ex <- generate_example_data(n_regions = 16, seed = 1)
#' dat <- data.frame(cases = ex$cases, population = ex$population,
#'                   region_id = ex$region_id, x = ex$x, y = ex$y,
#'                   node_id = ex$node_id)
#' res <- treespatial_scan(dat, cases = cases,
#'                         population = population,
#'                         region_id = region_id, x = x, y = y,
#'                         node_id = node_id, tree = ex$tree,
#'                         nsim = 19, seed = 42)
#' most_likely_cluster(res)
#'
#' @export
most_likely_cluster <- function(x, ...) {
  UseMethod("most_likely_cluster")
}

#' @rdname most_likely_cluster
#' @export
most_likely_cluster.treess <- function(x, ...) {
  x$most_likely_cluster
}

#' @rdname most_likely_cluster
#' @export
most_likely_cluster.sequential_scan <- function(x, ...) {
  cl <- x$clusters
  if (is.null(cl) || nrow(cl) == 0L) {
    return(NULL)
  }
  row <- cl[1L, , drop = FALSE]
  out <- as.list(row)
  # Unwrap list-columns (region_ids / leaf_ids) so the returned value
  # has the same shape as for the other scan types.
  for (nm in c("region_ids", "leaf_ids")) {
    if (nm %in% names(out)) out[[nm]] <- out[[nm]][[1L]]
  }
  out
}

#' Extract Secondary Clusters from a Scan Result
#'
#' Returns the table of secondary cluster candidates stored in a scan
#' result. Note that these are unfiltered candidates: overlapping
#' candidates are not removed and, for \code{circular_scan} and
#' \code{treespatial_scan}, per-candidate p-values are not attached.
#' Use \code{\link{filter_clusters}} to obtain distinct significant
#' clusters.
#'
#' @param x A scan result object (see
#'   \code{\link{most_likely_cluster}}).
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return A \code{data.frame}. For \code{circular_scan} and
#'   \code{treespatial_scan}, the top secondary candidates ordered by
#'   decreasing log-likelihood ratio. For \code{tree_scan}, the table
#'   of significant cuts at the alpha used in the original call. For
#'   \code{sequential_scan}, the per-iteration cluster table beyond the
#'   first iteration (the first iteration is the most likely cluster;
#'   see \code{\link{most_likely_cluster}}).
#'
#' @seealso \code{\link{filter_clusters}},
#'   \code{\link{most_likely_cluster}}
#'
#' @examples
#' ex <- generate_example_data(n_regions = 16, seed = 1)
#' dat <- data.frame(cases = ex$cases, population = ex$population,
#'                   region_id = ex$region_id, x = ex$x, y = ex$y,
#'                   node_id = ex$node_id)
#' res <- treespatial_scan(dat, cases = cases,
#'                         population = population,
#'                         region_id = region_id, x = x, y = y,
#'                         node_id = node_id, tree = ex$tree,
#'                         nsim = 19, seed = 42)
#' head(secondary_clusters(res))
#'
#' @export
secondary_clusters <- function(x, ...) {
  UseMethod("secondary_clusters")
}

#' @rdname secondary_clusters
#' @export
secondary_clusters.treess <- function(x, ...) {
  x$secondary_clusters
}

#' @rdname secondary_clusters
#' @export
secondary_clusters.tree_scan <- function(x, ...) {
  x$significant_cuts
}

#' @rdname secondary_clusters
#' @export
secondary_clusters.sequential_scan <- function(x, ...) {
  cl <- x$clusters
  if (is.null(cl) || nrow(cl) <= 1L) {
    return(cl[0L, , drop = FALSE])
  }
  cl[-1L, , drop = FALSE]
}

#' Extract the Monte Carlo P-Value from a Scan Result
#'
#' Returns the Monte Carlo p-value(s) of a scan result without
#' requiring the user to access the internal list structure.
#'
#' @param x A scan result object (see
#'   \code{\link{most_likely_cluster}}).
#' @param ... Further arguments passed to methods (currently unused).
#'
#' @return For \code{circular_scan}, \code{tree_scan}, and
#'   \code{treespatial_scan}, a single numeric value: the Monte Carlo
#'   p-value of the most likely cluster. For \code{sequential_scan},
#'   a named numeric vector with one p-value per iteration.
#'
#' @seealso \code{\link{most_likely_cluster}}
#'
#' @examples
#' ex <- generate_example_data(n_regions = 16, seed = 1)
#' dat <- data.frame(cases = ex$cases, population = ex$population,
#'                   region_id = ex$region_id, x = ex$x, y = ex$y,
#'                   node_id = ex$node_id)
#' res <- treespatial_scan(dat, cases = cases,
#'                         population = population,
#'                         region_id = region_id, x = x, y = y,
#'                         node_id = node_id, tree = ex$tree,
#'                         nsim = 19, seed = 42)
#' pvalue(res)
#'
#' @export
pvalue <- function(x, ...) {
  UseMethod("pvalue")
}

#' @rdname pvalue
#' @export
pvalue.treess <- function(x, ...) {
  x$pvalue
}

#' @rdname pvalue
#' @export
pvalue.sequential_scan <- function(x, ...) {
  cl <- x$clusters
  if (is.null(cl) || nrow(cl) == 0L) {
    return(numeric(0))
  }
  stats::setNames(cl$pvalue, paste0("iteration_", cl$iteration))
}
