#' Filter Overlapping Secondary Clusters
#'
#' Removes overlapping secondary clusters from the output of
#' \code{\link{treespatial_scan}}, \code{\link{circular_scan}}, or
#' \code{\link{tree_scan}}.
#'
#' @param result An object of class \code{"treespatial_scan"},
#'   \code{"circular_scan"}, or \code{"tree_scan"}.
#' @param alpha Numeric. Significance level for filtering. Default is
#'   \code{NULL}, which uses the alpha from the original analysis.
#'
#' @return A \code{data.frame} of distinct significant clusters ordered by
#'   descending LLR.
#'
#' @details
#' For \code{treespatial_scan} objects, two clusters overlap if they have
#' BOTH tree overlap (one node is ancestor/descendant of the other) AND
#' spatial overlap (same center or identical set of regions). This follows
#' Cancado et al. (2025).
#'
#' For \code{circular_scan} objects, the criterion follows Kulldorff (1997):
#' a secondary cluster is distinct if its center is not contained in any
#' previously retained cluster's zone, and no previously retained cluster's
#' center is contained in the candidate zone.
#'
#' For \code{tree_scan} objects, a secondary cluster is distinct if its node
#' is NOT an ancestor or descendant of any previously retained node. This
#' follows Kulldorff et al. (2003).
#'
#' When a dominant signal saturates the candidate pool (making distinct
#' secondary clusters hard to find), consider using
#' \code{\link{iterative_scan}} instead, which re-runs the scan after
#' removing each detected cluster.
#'
#' @references
#' Kulldorff, M. (1997). A spatial scan statistic. \emph{Communications in
#' Statistics - Theory and Methods}, 26(6), 1481-1496.
#'
#' Kulldorff, M., Fang, Z., & Walsh, S. J. (2003). A tree-based scan
#' statistic for database disease surveillance. \emph{Biometrics}, 59(2),
#' 323-331.
#'
#' Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
#' (2025). A tree-spatial scan statistic. \emph{Environmental and Ecological
#' Statistics}, 32, 953-978. \doi{10.1007/s10651-025-00670-w}
#'
#' @seealso \code{\link{treespatial_scan}}, \code{\link{circular_scan}},
#'   \code{\link{tree_scan}}, \code{\link{iterative_scan}}
#'
#' @export
filter_clusters <- function(result, alpha = NULL) {

  is_treespatial <- inherits(result, "treespatial_scan")
  is_circular <- inherits(result, "circular_scan")
  is_tree <- inherits(result, "tree_scan")

  if (!is_treespatial && !is_circular && !is_tree) {
    stop("'result' must be of class 'treespatial_scan', 'circular_scan', ",
         "or 'tree_scan'.", call. = FALSE)
  }

  if (is.null(alpha)) alpha <- result$alpha

  # --- tree_scan: use all_cuts (p-values already computed) ---
  if (is_tree) {
    sig <- result$all_cuts[result$all_cuts$pvalue < alpha, ]
    if (nrow(sig) == 0) {
      message("No significant clusters found at alpha = ", alpha, ".")
      return(data.frame())
    }
    sig <- sig[order(-sig$llr), ]
    return(.filter_tree_only(sig, result$tree))
  }

  # --- circular_scan / treespatial_scan: use secondary_clusters ---
  all_pairs <- result$secondary_clusters
  if (is.null(all_pairs) || nrow(all_pairs) == 0) {
    message("No secondary clusters to filter.")
    return(data.frame())
  }

  # Compute p-values (vectorized)
  sim_llr <- result$simulated_llr
  nsim <- result$nsim
  sim_sorted <- sort(sim_llr)
  n_geq <- nsim - findInterval(all_pairs$llr, sim_sorted, left.open = TRUE)
  all_pairs$pvalue <- (n_geq + 1) / (nsim + 1)

  sig_pairs <- all_pairs[all_pairs$pvalue < alpha, ]
  if (nrow(sig_pairs) == 0) {
    message("No significant clusters found at alpha = ", alpha, ".")
    return(data.frame())
  }
  sig_pairs <- sig_pairs[order(-sig_pairs$llr), ]

  # Pre-compute zone regions
  zones <- build_zones(result$regions,
                       max_pop = result$total_population * 0.5)
  pair_regions <- vector("list", nrow(sig_pairs))
  for (i in seq_len(nrow(sig_pairs))) {
    zi <- .find_zone(zones, sig_pairs$center[i], sig_pairs$n_regions[i])
    if (!is.null(zi)) {
      pair_regions[[i]] <- zones[[zi]]$region_idx
    } else {
      pair_regions[[i]] <- integer(0)
    }
  }

  if (is_circular) {
    .filter_circular(sig_pairs, pair_regions)
  } else {
    .filter_treespatial(sig_pairs, pair_regions, result$tree)
  }
}


# ---- Tree-only criterion (Kulldorff et al. 2003) ----
# Distinct if NOT ancestor/descendant of any retained node.
.filter_tree_only <- function(sig_cuts, tree) {

  all_descendants <- list()
  for (nd in tree$node_id) {
    all_descendants[[nd]] <- .get_descendants(tree, nd)
  }

  retained <- list(sig_cuts[1, ])

  if (nrow(sig_cuts) > 1) {
    for (i in 2:nrow(sig_cuts)) {
      cand_node <- sig_cuts$node_id[i]
      cand_desc <- all_descendants[[cand_node]]

      is_distinct <- TRUE
      for (j in seq_along(retained)) {
        kept_node <- retained[[j]]$node_id
        kept_desc <- all_descendants[[kept_node]]

        if (cand_node %in% kept_desc || kept_node %in% cand_desc) {
          is_distinct <- FALSE
          break
        }
      }

      if (is_distinct) {
        retained[[length(retained) + 1]] <- sig_cuts[i, ]
      }
    }
  }

  do.call(rbind, retained)
}


# ---- Kulldorff (1997) criterion for circular scan ----
.filter_circular <- function(sig_pairs, pair_regions) {

  retained <- list(sig_pairs[1, ])
  retained_regions <- list(pair_regions[[1]])
  retained_centers <- sig_pairs$center[1]

  if (nrow(sig_pairs) > 1) {
    for (i in 2:nrow(sig_pairs)) {
      cand_center <- sig_pairs$center[i]
      cand_regs <- pair_regions[[i]]

      is_distinct <- TRUE
      for (j in seq_along(retained)) {
        if (cand_center %in% retained_regions[[j]] ||
            retained_centers[j] %in% cand_regs) {
          is_distinct <- FALSE
          break
        }
      }

      if (is_distinct) {
        retained[[length(retained) + 1]] <- sig_pairs[i, ]
        retained_regions[[length(retained_regions) + 1]] <- cand_regs
        retained_centers <- c(retained_centers, cand_center)
      }
    }
  }

  do.call(rbind, retained)
}


# ---- Tree-spatial overlap criterion (Cancado et al. 2025, Sec 5.1.1) ----
# A candidate (z, g) is REJECTED as redundant if:
#   * the candidate's branch g overlaps with the retained branch (one is
#     ancestor/descendant of the other in the tree), AND
#   * the candidate's spatial zone intersects the retained zone (any
#     region in common).
# Equivalently, a candidate is RETAINED as distinct if at least one of:
#   * branches are disjoint (no ancestor/descendant relation), OR
#   * spatial zones are disjoint (no region in common).
# This matches the three scenarios listed in Sec. 5.1.1 of Cancado et al.
.filter_treespatial <- function(sig_pairs, pair_regions, tree) {

  all_descendants <- list()
  for (nd in tree$node_id) {
    all_descendants[[nd]] <- .get_descendants(tree, nd)
  }

  retained <- list(sig_pairs[1, ])
  retained_regions <- list(pair_regions[[1]])

  if (nrow(sig_pairs) > 1) {
    for (i in 2:nrow(sig_pairs)) {
      candidate <- sig_pairs[i, ]
      cand_node <- candidate$node_id
      cand_desc <- all_descendants[[cand_node]]
      cand_regs <- pair_regions[[i]]

      is_distinct <- TRUE
      for (j in seq_along(retained)) {
        kept <- retained[[j]]
        kept_node <- kept$node_id
        kept_desc <- all_descendants[[kept_node]]
        kept_regs <- retained_regions[[j]]

        # Branch overlap: same node OR ancestor/descendant
        node_overlap <- (cand_node == kept_node) ||
                         (cand_node %in% kept_desc) ||
                         (kept_node %in% cand_desc)
        # Spatial overlap: any region in common
        spatial_overlap <- length(intersect(cand_regs, kept_regs)) > 0

        if (node_overlap && spatial_overlap) {
          is_distinct <- FALSE
          break
        }
      }

      if (is_distinct) {
        retained[[length(retained) + 1]] <- candidate
        retained_regions[[length(retained_regions) + 1]] <- cand_regs
      }
    }
  }

  do.call(rbind, retained)
}
