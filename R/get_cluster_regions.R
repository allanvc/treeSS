#' Extract Cluster Membership for Each Region
#'
#' Returns a data.frame indicating which cluster (if any) each region
#' belongs to. This is the primary output for visualization — use it with
#' any mapping package (ggplot2, leaflet, tmap, etc.).
#'
#' This is a generic with methods for objects returned by
#' \code{\link{treespatial_scan}}, \code{\link{circular_scan}}, and
#' \code{\link{iterative_scan}}.
#'
#' @param result Either a \code{"treespatial_scan"}/\code{"circular_scan"}
#'   object (single-pass scan) or an \code{"iterative_scan"} object.
#' @param n_clusters Integer. (Single-pass methods only.) Number of
#'   clusters to extract. \code{1} returns only the most likely cluster.
#'   Values greater than 1 use \code{\link{filter_clusters}} internally
#'   to identify distinct secondary clusters. Default is \code{1}.
#'   Ignored for iterative-scan inputs (every iteration is returned).
#' @param overlap Logical. If \code{TRUE} (default), returns a facet-ready
#'   data.frame: all regions are replicated for each cluster panel, with a
#'   \code{panel} column for faceting and \code{cluster} marked only for
#'   member regions. If \code{FALSE}, returns one row per region assigned
#'   to its highest-ranked cluster (for single-panel maps).
#' @param ... Further arguments passed to methods.
#'
#' @return A \code{data.frame} with columns from \code{result$regions} plus:
#'   \describe{
#'     \item{cluster}{Integer cluster number (1 = most likely / first
#'       iteration, 2 = first secondary / second iteration, etc.), or
#'       \code{NA} if the region is not in the cluster.}
#'     \item{node_id}{The tree node of the cluster, or \code{NA}.}
#'     \item{llr}{The log-likelihood ratio of the cluster, or \code{NA}.}
#'     \item{pvalue}{The p-value of the cluster, or \code{NA}.}
#'     \item{pvalue_adjusted, significant}{(Iterative method only) The
#'       Holm-Bonferroni adjusted p-value and corresponding significance
#'       flag for the iteration.}
#'     \item{panel}{(Only when \code{overlap = TRUE}) A label for
#'       \code{facet_wrap}, e.g. \code{"#1 P209 (LR=39.6)"} for single-
#'       pass scans, or \code{"Iter 1: P209 (LR=39.6, p_adj=0.005)"} for
#'       iterative scans.}
#'   }
#'
#' @seealso \code{\link{treespatial_scan}}, \code{\link{circular_scan}},
#'   \code{\link{iterative_scan}}, \code{\link{filter_clusters}}
#'
#' @export
get_cluster_regions <- function(result, n_clusters = 1L, overlap = TRUE,
                                 ...) {
  UseMethod("get_cluster_regions")
}


#' @rdname get_cluster_regions
#' @export
get_cluster_regions.default <- function(result, n_clusters = 1L,
                                          overlap = TRUE, ...) {

  if (!inherits(result, "treespatial_scan") &&
      !inherits(result, "circular_scan")) {
    stop("'result' must be a treespatial_scan, circular_scan, or ",
         "iterative_scan object.", call. = FALSE)
  }
  .get_cluster_regions_singlepass(result, n_clusters, overlap)
}


#' @rdname get_cluster_regions
#' @export
get_cluster_regions.iterative_scan <- function(result, n_clusters = 1L,
                                                 overlap = TRUE, ...) {

  iter <- result
  if (iter$n_iter == 0) {
    return(data.frame())
  }
  if (is.null(iter$regions)) {
    stop("Iterative scan has no $regions table (tree-only scan). ",
         "Mapping by region is not applicable.", call. = FALSE)
  }

  full_regions <- iter$regions
  parts <- vector("list", iter$n_iter)

  for (k in seq_len(iter$n_iter)) {
    cl_regs <- iter$clusters$region_ids[[k]]
    in_cl   <- full_regions$region_id %in% cl_regs

    cr_k <- cbind(
      full_regions,
      data.frame(
        cluster         = ifelse(in_cl, k, NA_integer_),
        node_id         = ifelse(in_cl,
                                  as.character(iter$clusters$node_id[k]),
                                  NA_character_),
        llr             = ifelse(in_cl, iter$clusters$llr[k], NA_real_),
        pvalue          = ifelse(in_cl, iter$clusters$pvalue[k], NA_real_),
        pvalue_adjusted = ifelse(in_cl, iter$clusters$pvalue_adjusted[k],
                                  NA_real_),
        significant     = ifelse(in_cl, iter$clusters$significant[k], NA),
        stringsAsFactors = FALSE
      )
    )

    if (overlap) {
      sig_marker <- if (isTRUE(iter$clusters$significant[k])) "" else " (n.s.)"
      cr_k$panel <- paste0("Iter ", k, ": ", iter$clusters$node_id[k],
                           " (LR=", round(iter$clusters$llr[k], 1),
                           ", p_adj=",
                           format.pval(iter$clusters$pvalue_adjusted[k],
                                        digits = 2),
                           sig_marker, ")")
    }

    parts[[k]] <- cr_k
  }

  if (overlap) {
    do.call(rbind, parts)
  } else {
    # Non-overlap: one row per region, using the LOWEST iteration in which
    # it appears (mirrors single-pass non-overlap semantics).
    full <- do.call(rbind, parts)
    full <- full[order(full$cluster, na.last = TRUE), ]
    full[!duplicated(full$region_id), ]
  }
}


# Internal: original single-pass logic moved out of the generic
.get_cluster_regions_singlepass <- function(result, n_clusters, overlap) {

  regions <- result$regions
  n <- nrow(regions)

  # --- Collect all cluster region assignments ---
  cluster_list <- list()

  # Cluster 1: most likely
  mlc <- result$most_likely_cluster
  mlc_idx <- match(mlc$region_ids, regions$region_id)
  mlc_idx <- mlc_idx[!is.na(mlc_idx)]
  mlc_node <- if (!is.null(mlc$node_id)) as.character(mlc$node_id) else "Zone"
  cluster_list[[1]] <- list(
    idx     = mlc_idx,
    cluster = 1L,
    node_id = mlc_node,
    llr     = mlc$llr,
    pvalue  = result$pvalue
  )

  # Secondary clusters
  if (n_clusters > 1) {
    dc <- filter_clusters(result)
    if (!is.null(dc) && nrow(dc) > 1) {
      zones <- build_zones(regions, max_pop = result$total_population * 0.5)
      sec_rows <- dc[-1, , drop = FALSE]
      n_sec <- min(nrow(sec_rows), n_clusters - 1)
      has_node <- "node_id" %in% names(sec_rows)

      for (k in seq_len(n_sec)) {
        zi <- .find_zone(zones, sec_rows$center[k], sec_rows$n_regions[k])
        if (!is.null(zi)) {
          sec_node <- if (has_node) as.character(sec_rows$node_id[k]) else "Zone"
          cluster_list[[k + 1]] <- list(
            idx     = zones[[zi]]$region_idx,
            cluster = as.integer(k + 1),
            node_id = sec_node,
            llr     = sec_rows$llr[k],
            pvalue  = sec_rows$pvalue[k]
          )
        }
      }
    }
  }

  # --- Build output ---
  if (overlap) {
    # Facet-ready format: ALL regions replicated for each cluster panel.
    # In each panel, cluster members are marked, non-members have NA.
    # Ready for: merge(map, cr) then ggplot + facet_wrap(~ panel)
    rows <- list()

    for (cl in cluster_list) {
      chunk <- regions
      chunk$cluster <- NA_integer_
      chunk$node_id <- NA_character_
      chunk$llr     <- NA_real_
      chunk$pvalue  <- NA_real_
      chunk$panel   <- paste0("#", cl$cluster, " ", cl$node_id,
                              " (LR=", round(cl$llr, 1), ")")

      # Mark cluster members
      chunk$cluster[cl$idx] <- cl$cluster
      chunk$node_id[cl$idx] <- cl$node_id
      chunk$llr[cl$idx]     <- cl$llr
      chunk$pvalue[cl$idx]  <- cl$pvalue

      rows <- c(rows, list(chunk))
    }

    out <- do.call(rbind, rows)
    rownames(out) <- NULL

  } else {
    # Single-map format: one row per region, highest-ranked cluster wins.
    out <- regions
    out$cluster <- NA_integer_
    out$node_id <- NA_character_
    out$llr     <- NA_real_
    out$pvalue  <- NA_real_

    for (cl in cluster_list) {
      free <- is.na(out$cluster[cl$idx])
      out$cluster[cl$idx[free]] <- cl$cluster
      out$node_id[cl$idx[free]] <- cl$node_id
      out$llr[cl$idx[free]]     <- cl$llr
      out$pvalue[cl$idx[free]]  <- cl$pvalue
    }
  }

  out
}

#' @keywords internal
.find_zone <- function(zones, center, n_regions) {
  for (i in seq_along(zones)) {
    if (zones[[i]]$center == center &&
        length(zones[[i]]$region_idx) == n_regions) {
      return(i)
    }
  }
  NULL
}
