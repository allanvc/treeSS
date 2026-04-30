#' Iterative Scan
#'
#' Runs the scan statistic iteratively, removing the cases of each detected
#' cluster before the next iteration. Supports tree-spatial, circular, and
#' tree-only scans, dispatched by which arguments are supplied.
#'
#' \strong{Note on methodology.} This iterative procedure is \emph{not} part
#' of the original tree-spatial scan statistic of Cancado et al. (2025). It is
#' an extension inspired by the conditional scan approach of Zhang et al.
#' (2010), which sequentially removes cases attributed to detected clusters
#' before re-running the scan to find additional, distinct anomalies. To
#' reproduce the secondary-cluster procedure described in Section 5.1.1 of
#' Cancado et al. (2025), use \code{\link{filter_clusters}} or
#' \code{\link{get_cluster_regions}} with \code{n_clusters > 1} on a single
#' scan result instead.
#'
#' \strong{Multiple-testing correction.} Because each iteration is a separate
#' hypothesis test on data that has been modified by the previous iteration,
#' the raw p-values overstate significance. This function collects raw p-values
#' from all iterations performed and applies the Holm-Bonferroni correction
#' (\code{\link[stats]{p.adjust}} with \code{method = "holm"}) at the end. The
#' returned \code{clusters} data.frame includes both the raw p-value
#' (\code{pvalue}) and the adjusted p-value (\code{pvalue_adjusted}), plus a
#' logical \code{significant} column indicating whether the cluster is
#' significant after correction at level \code{alpha}.
#'
#' The loop stops early only when the residual signal is exhausted (the most
#' likely cluster has \code{LR = 0} or zero cases), not based on raw p-values,
#' since the multiple-testing correction depends on the full set of tests.
#'
#' @param cases Numeric vector. For tree-spatial scan: one entry per
#'   (region, leaf) row. For circular: one entry per region. For tree-only:
#'   one entry per leaf.
#' @param population Numeric vector parallel to \code{cases}.
#' @param region_id,x,y Vectors parallel to \code{cases} (omit for
#'   tree-only scan).
#' @param node_id Vector parallel to \code{cases} (omit for circular and
#'   tree-only scans).
#' @param tree Tree as a 2-column data.frame (\code{node_id, parent_id}),
#'   or use \code{tree_node_id}/\code{tree_parent_id}. Omit for circular
#'   scan.
#' @param tree_node_id,tree_parent_id Optional. Parallel vectors as an
#'   alternative to \code{tree}.
#' @param max_iter Integer. Maximum number of iterations.
#' @param alpha Numeric. Significance threshold applied to Holm-Bonferroni
#'   adjusted p-values.
#' @param nsim Integer. MC simulations per iteration.
#' @param max_pop_pct Numeric. Passed to inner scans.
#' @param model Character. \code{"poisson"} or \code{"binomial"}.
#' @param seed Integer or \code{NULL}.
#' @param verbose Logical.
#' @param n_cores Integer. OpenMP threads.
#'
#' @return An object of class \code{"iterative_scan"} with components:
#'   \describe{
#'     \item{clusters}{A data.frame with one row per iteration, containing
#'       the raw p-value (\code{pvalue}), the Holm-Bonferroni adjusted
#'       p-value (\code{pvalue_adjusted}), a logical \code{significant}
#'       column, plus the cluster's node, regions, cases, expected, RR, and
#'       LLR.}
#'     \item{iterations}{A list with the full scan result of each iteration.}
#'     \item{regions, tree, alpha, n_iter, scan_type}{Bookkeeping fields.}
#'   }
#'
#' @references
#' Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L. H.
#' (2025). A tree-spatial scan statistic. \emph{Environmental and Ecological
#' Statistics}, 32, 953-978.
#'
#' Kulldorff, M. (1997). A spatial scan statistic. \emph{Communications in
#' Statistics - Theory and Methods}, 26(6), 1481-1496.
#'
#' Zhang, Z., Assuncao, R., & Kulldorff, M. (2010). Spatial scan statistics
#' adjusted for multiple clusters. \emph{Journal of Probability and
#' Statistics}, 2010, 642379.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65-70.
#'
#' @seealso \code{\link{treespatial_scan}}, \code{\link{circular_scan}},
#'   \code{\link{tree_scan}}, \code{\link{filter_clusters}},
#'   \code{\link{get_cluster_regions}}
#'
#' @export
#' @examples
#' \donttest{
#' data(london_collisions); data(london_tree)
#' result <- iterative_scan(
#'   cases       = london_collisions$cases,
#'   population  = london_collisions$population,
#'   region_id   = london_collisions$region_id,
#'   x           = london_collisions$x,
#'   y           = london_collisions$y,
#'   node_id     = london_collisions$node_id,
#'   tree        = london_tree,
#'   max_iter = 3, nsim = 99, seed = 42
#' )
#' print(result)
#' }
iterative_scan <- function(cases = NULL, population = NULL,
                           region_id = NULL, x = NULL, y = NULL,
                           node_id = NULL,
                           tree = NULL,
                           tree_node_id = NULL, tree_parent_id = NULL,
                           max_iter = 5L, alpha = 0.05,
                           nsim = 999L, max_pop_pct = 0.5,
                           model = c("poisson", "binomial"),
                           seed = NULL, verbose = TRUE,
                           n_cores = 1L) {

  model <- match.arg(model)

  has_tree <- !is.null(tree) || (!is.null(tree_node_id) &&
                                  !is.null(tree_parent_id))
  has_node <- !is.null(node_id)
  has_geo  <- !is.null(region_id) && !is.null(x) && !is.null(y)

  if (has_tree && has_node && has_geo) {
    scan_type <- "treespatial"
    tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)
  } else if (!has_tree && has_geo) {
    scan_type <- "circular"
  } else if (has_tree && !has_node && !has_geo) {
    scan_type <- "tree"
    tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)
  } else {
    stop("Inputs are not consistent. Provide:\n",
         "  tree-spatial: cases, population, region_id, x, y, node_id, tree\n",
         "  circular:     cases, population, region_id, x, y\n",
         "  tree-only:    cases, population, tree (no region_id/x/y/node_id)",
         call. = FALSE)
  }

  # State carried across iterations
  if (scan_type == "treespatial") {
    prepared <- .build_inputs(cases, population, region_id, x, y,
                              node_id, tree)
    regions       <- prepared$regions
    cases_current <- prepared$cases_matrix    # matrix
  } else if (scan_type == "circular") {
    regions       <- .build_regions_circular(cases, population,
                                              region_id, x, y)
    cases_current <- as.numeric(cases)        # vector
  } else {
    # tree-only
    if (is.null(population)) {
      stop("Tree-only scan requires 'population' (vector of leaf-level ",
           "denominators).", call. = FALSE)
    }
    cases_current <- as.numeric(cases)
    regions <- NULL
  }

  clusters_rows <- list()
  iterations <- list()

  for (k in seq_len(max_iter)) {
    iter_seed <- if (is.null(seed)) NULL else seed + k - 1L

    if (verbose) message(sprintf("Iteration %d/%d ...", k, max_iter))

    if (scan_type == "treespatial") {
      # Convert cases matrix back to vectors for treespatial_scan
      vec <- .matrix_to_vectors(cases_current, regions, tree)
      res <- treespatial_scan(
        cases       = vec$cases,
        population  = vec$population,
        region_id   = vec$region_id,
        x           = vec$x,
        y           = vec$y,
        node_id     = vec$node_id,
        tree        = tree,
        max_pop_pct = max_pop_pct,
        nsim        = nsim, alpha = alpha,
        model       = model,
        seed        = iter_seed, n_cores = n_cores
      )
    } else if (scan_type == "circular") {
      res <- circular_scan(
        cases       = cases_current,
        population  = regions$population,
        region_id   = regions$region_id,
        x           = regions$x,
        y           = regions$y,
        max_pop_pct = max_pop_pct,
        nsim        = nsim, alpha = alpha,
        model       = model,
        seed        = iter_seed, n_cores = n_cores
      )
    } else {
      res <- tree_scan(tree, cases_current, population = population,
                       nsim = nsim, alpha = alpha,
                       model = model,
                       seed = iter_seed, n_cores = n_cores)
    }

    mlc <- res$most_likely_cluster
    pvalue <- res$pvalue

    # Stop only if the residual signal is gone (no cluster found, LLR=0).
    # Otherwise, continue collecting iterations - significance is decided
    # AT THE END via Holm-Bonferroni multiple-testing correction.
    if (is.null(mlc$llr) || mlc$llr <= 0 ||
        is.null(mlc$cases) || mlc$cases == 0) {
      if (verbose) {
        message(sprintf("  Iteration %d: no residual signal (LR=0). Stopping.",
                        k))
      }
      break
    }

    # --- Record this cluster ---
    if (scan_type == "tree") {
      row <- data.frame(
        iteration  = k,
        node_id    = as.character(mlc$node_id),
        cases      = mlc$cases,
        expected   = mlc$expected,
        population = mlc$population,
        llr        = mlc$llr,
        pvalue     = pvalue,
        stringsAsFactors = FALSE
      )
      row$leaf_ids <- I(list(mlc$leaf_ids))
    } else {
      row <- data.frame(
        iteration  = k,
        node_id    = if (has_tree) as.character(mlc$node_id) else NA_character_,
        n_regions  = length(mlc$region_ids),
        cases      = mlc$cases,
        expected   = mlc$expected,
        population = mlc$population,
        rr         = mlc$rr,
        llr        = mlc$llr,
        pvalue     = pvalue,
        stringsAsFactors = FALSE
      )
      row$region_ids <- I(list(mlc$region_ids))
    }

    clusters_rows[[k]] <- row
    iterations[[k]] <- res

    if (verbose) {
      if (scan_type == "tree") {
        message(sprintf("  Cluster %d: node = %s, LR = %.2f, raw p = %.3f",
                        k, mlc$node_id, mlc$llr, pvalue))
      } else if (scan_type == "treespatial") {
        message(sprintf("  Cluster %d: node = %s, %d regions, LR = %.2f, raw p = %.3f",
                        k, mlc$node_id, length(mlc$region_ids), mlc$llr, pvalue))
      } else {
        message(sprintf("  Cluster %d: %d regions, LR = %.2f, raw p = %.3f",
                        k, length(mlc$region_ids), mlc$llr, pvalue))
      }
    }

    # --- Remove cases for next iteration ---
    cases_current <- .zero_cluster_cases(cases_current, mlc, tree, regions,
                                         scan_type)
  }

  # --- Apply Holm-Bonferroni correction across the m collected iterations ---
  if (length(clusters_rows) > 0) {
    raw_p   <- vapply(clusters_rows, function(r) r$pvalue, numeric(1))
    adj_p   <- stats::p.adjust(raw_p, method = "holm")
    sig_vec <- adj_p < alpha
    for (k in seq_along(clusters_rows)) {
      clusters_rows[[k]]$pvalue_adjusted <- adj_p[k]
      clusters_rows[[k]]$significant     <- sig_vec[k]
    }
    if (verbose) {
      message(sprintf("\nHolm-Bonferroni adjustment over %d iteration(s):",
                      length(clusters_rows)))
      for (k in seq_along(clusters_rows)) {
        message(sprintf("  Iter %d: raw p = %.3f, adjusted p = %.3f -> %s",
                        k, raw_p[k], adj_p[k],
                        if (sig_vec[k]) "SIGNIFICANT" else "not significant"))
      }
    }
  }

  # For tree-only mode, regions wasn't created in this scope
  out_regions <- if (scan_type == "tree") NULL else regions

  out <- list(
    clusters   = if (length(clusters_rows) > 0) do.call(rbind, clusters_rows)
                 else data.frame(),
    iterations = iterations,
    regions    = out_regions,
    tree       = tree,
    alpha      = alpha,
    n_iter     = length(clusters_rows),
    scan_type  = scan_type
  )
  class(out) <- "iterative_scan"
  out
}


# ---- Helper: zero out cluster cases for next iteration ----
#' @keywords internal
.zero_cluster_cases <- function(cases, mlc, tree, regions, scan_type) {

  if (scan_type == "circular") {
    region_idx <- match(mlc$region_ids, regions$region_id)
    region_idx <- region_idx[!is.na(region_idx)]
    cases[region_idx] <- 0
    return(cases)
  }

  if (scan_type == "tree") {
    # Tree-only: zero leaves under the detected node
    leaves <- .get_leaves(tree)
    desc <- .get_descendants(tree, as.character(mlc$node_id))
    desc_leaves <- intersect(desc, leaves)
    leaf_idx <- match(desc_leaves, leaves)
    leaf_idx <- leaf_idx[!is.na(leaf_idx)]
    cases[leaf_idx] <- 0
    return(cases)
  }

  # Tree-spatial: zero leaf descendants x regions
  node_id <- as.character(mlc$node_id)
  descendants <- .get_descendants(tree, node_id)

  parents_set <- unique(tree$parent_id[!is.na(tree$parent_id)])
  all_leaves <- setdiff(tree$node_id, parents_set)
  desc_leaves <- intersect(descendants, all_leaves)

  leaf_idx <- match(desc_leaves, rownames(cases))
  leaf_idx <- leaf_idx[!is.na(leaf_idx)]

  region_idx <- match(mlc$region_ids, regions$region_id)
  region_idx <- region_idx[!is.na(region_idx)]

  if (length(leaf_idx) > 0 && length(region_idx) > 0) {
    cases[leaf_idx, region_idx] <- 0L
  }
  cases
}


#' Print Method for iterative_scan Objects
#'
#' @param x An object of class \code{"iterative_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.iterative_scan <- function(x, ...) {
  cat("Iterative Scan\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Scan type:", x$scan_type, "\n")
  cat("Iterations performed:", x$n_iter, "\n")
  cat("Alpha:", x$alpha,
      "  (applied to Holm-Bonferroni adjusted p-values)\n\n")

  if (nrow(x$clusters) > 0) {
    n_sig <- sum(x$clusters$significant, na.rm = TRUE)
    cat("Significant clusters after Holm-Bonferroni: ", n_sig,
        " of ", nrow(x$clusters), "\n\n", sep = "")

    cat("All iterations:\n")
    display <- x$clusters
    display$region_ids <- NULL
    display$leaf_ids <- NULL
    print(display, row.names = FALSE)
  } else {
    cat("No significant clusters detected.\n")
  }

  invisible(x)
}
