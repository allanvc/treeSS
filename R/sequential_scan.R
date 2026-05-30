#' Sequential Scan for Secondary Clusters
#'
#' Implements the sequential adjustment of Zhang, Assuncao and Kulldorff
#' (2010) for detecting secondary clusters. After the most likely cluster
#' (MLC) is detected and found significant, the regions composing it
#' (optionally together with a buffer of nearest-neighbour regions) are
#' \emph{removed} from the dataset -- treated as a \dQuote{lake} with no
#' population and no cases -- and the scan is re-run on the reduced data
#' with a fresh Monte Carlo simulation. The procedure is iterated until
#' the MLC of the current reduced data is no longer significant, or
#' \code{max_iter} is reached.
#'
#' This differs from \code{\link{filter_clusters}}, which scans the
#' single-pass candidate pool for non-overlapping clusters (the original
#' Cancado et al. 2025 criterion). The sequential approach computes a
#' fresh, unconservative p-value for each secondary cluster, conditional
#' on the prior detections.
#'
#' Dispatch on the supplied arguments selects one of three scan types:
#' \itemize{
#'   \item tree-spatial: all of \code{cases}, \code{population},
#'         \code{region_id}, \code{x}, \code{y}, \code{node_id},
#'         \code{tree} are required;
#'   \item circular: \code{cases}, \code{population}, \code{region_id},
#'         \code{x}, \code{y} (no \code{node_id} / \code{tree});
#'   \item tree-only: \code{cases}, \code{population}, \code{tree}
#'         (no spatial inputs). \code{buffer_size} is ignored for the
#'         tree-only case.
#' }
#'
#' @param cases Numeric vector. For tree-spatial scan: one entry per
#'   (region, leaf) row. For circular: one entry per region. For
#'   tree-only: one entry per leaf.
#' @param population Numeric vector parallel to \code{cases}.
#' @param region_id,x,y Vectors parallel to \code{cases} (omit for
#'   tree-only scan).
#' @param node_id Vector parallel to \code{cases} (omit for circular and
#'   tree-only scans).
#' @param tree Tree as a 2-column data.frame
#'   (\code{node_id, parent_id}), or use
#'   \code{tree_node_id}/\code{tree_parent_id}. Omit for circular scan.
#' @param tree_node_id,tree_parent_id Optional. Parallel vectors as an
#'   alternative to \code{tree}.
#' @param max_iter Integer. Maximum number of iterations (including the
#'   primary MLC). Default 5.
#' @param alpha Numeric. Significance level. Iteration stops when the
#'   MLC p-value of the current reduced data exceeds \code{alpha}.
#'   Default 0.05.
#' @param nsim Integer. Monte Carlo replications per iteration.
#'   Default 999.
#' @param max_pop_pct Numeric. Maximum zone-size constraint, passed to
#'   the inner scans. Default 0.5.
#' @param buffer_size Integer. Number of nearest-neighbour regions to
#'   remove together with each detected cluster, computed by Euclidean
#'   distance from each cluster region to the remaining ones. Default
#'   0. Zhang et al. (2010) report that the type I error and power are
#'   essentially insensitive to the buffer size and that
#'   \code{buffer_size = 0} is generally adequate. Ignored for
#'   tree-only scans.
#' @param model Character. \code{"poisson"} (default) or
#'   \code{"binomial"}.
#' @param seed Integer or \code{NULL}. Seed for the first iteration;
#'   subsequent iterations advance the RNG state to avoid using the
#'   same seed across iterations.
#' @param verbose Logical. Print progress at each iteration.
#' @param n_cores Integer. OpenMP threads passed to inner scans.
#'
#' @return An object of class \code{"sequential_scan"} with components:
#'   \describe{
#'     \item{clusters}{A data.frame with one row per iteration performed,
#'       with columns \code{iteration}, \code{node_id} (only for
#'       tree-spatial / tree-only scans), \code{n_regions} (only for
#'       scans with a spatial component), \code{cases}, \code{expected},
#'       \code{population}, \code{rr} (only for scans with a spatial
#'       component), \code{llr}, \code{pvalue}, and a logical
#'       \code{significant} flag (\code{pvalue <= alpha}). The list-
#'       columns \code{region_ids} and/or \code{leaf_ids} carry the
#'       cluster composition.}
#'     \item{iterations}{A list with the full scan result of each
#'       iteration.}
#'     \item{regions, tree, alpha, nsim, buffer_size, n_iter,
#'           scan_type}{Bookkeeping fields. \code{n_iter} is the number
#'       of iterations actually performed.}
#'   }
#'
#' @references
#' Zhang, Z., Assuncao, R., & Kulldorff, M. (2010). Spatial scan
#' statistics adjusted for multiple clusters. \emph{Journal of
#' Probability and Statistics}, 2010, 642379.
#'
#' Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal,
#' L. H. (2025). A tree-spatial scan statistic. \emph{Environmental
#' and Ecological Statistics}, 32, 953-978.
#' \doi{10.1007/s10651-025-00670-w}
#'
#' Kulldorff, M. (1997). A spatial scan statistic.
#' \emph{Communications in Statistics - Theory and Methods}, 26(6),
#' 1481-1496.
#'
#' @seealso \code{\link{treespatial_scan}}, \code{\link{circular_scan}},
#'   \code{\link{tree_scan}}, \code{\link{filter_clusters}},
#'   \code{\link{get_cluster_regions}}
#'
#' @export
#' @examples
#' \donttest{
#' data(london_collisions); data(london_tree)
#' result <- sequential_scan(
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
sequential_scan <- function(cases = NULL, population = NULL,
                            region_id = NULL, x = NULL, y = NULL,
                            node_id = NULL,
                            tree = NULL,
                            tree_node_id = NULL, tree_parent_id = NULL,
                            max_iter = 5L, alpha = 0.05,
                            nsim = 999L, max_pop_pct = 0.5,
                            buffer_size = 0L,
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

  stopifnot(is.numeric(max_iter),    length(max_iter)    == 1L,
            max_iter >= 1L)
  stopifnot(is.numeric(alpha),       length(alpha)       == 1L,
            alpha > 0,    alpha < 1)
  stopifnot(is.numeric(nsim),        length(nsim)        == 1L,
            nsim >= 1L)
  stopifnot(is.numeric(max_pop_pct), length(max_pop_pct) == 1L,
            max_pop_pct > 0, max_pop_pct <= 1)
  stopifnot(is.numeric(buffer_size), length(buffer_size) == 1L,
            buffer_size >= 0L)

  max_iter    <- as.integer(max_iter)
  buffer_size <- as.integer(buffer_size)

  # Working copies of the per-row vectors that we will shrink at each
  # iteration by dropping rows that belong to the union of detected
  # cluster regions (plus buffer). For the tree-only scan we drop leaves
  # instead.
  cur_cases      <- as.numeric(cases)
  cur_population <- as.numeric(population)

  if (scan_type != "tree") {
    cur_region_id <- region_id
    cur_x         <- as.numeric(x)
    cur_y         <- as.numeric(y)
    if (scan_type == "treespatial") {
      cur_node_id <- node_id
    }
  }

  removed_region_ids <- character(0)
  removed_leaf_ids   <- character(0)
  clusters_rows <- list()
  iterations    <- list()

  for (k in seq_len(max_iter)) {
    iter_seed <- if (is.null(seed)) NULL else seed + k - 1L

    if (verbose) message(sprintf("Iteration %d/%d ...", k, max_iter))

    # Early stop: no data left
    if (length(cur_cases) == 0L || sum(cur_cases) <= 0 ||
        sum(cur_population) <= 0) {
      if (verbose) message(sprintf(
        "  Iteration %d: reduced dataset exhausted. Stopping.", k))
      break
    }

    # Run the appropriate scan on the current reduced data
    if (scan_type == "treespatial") {
      res <- treespatial_scan(
        cases       = cur_cases,
        population  = cur_population,
        region_id   = cur_region_id,
        x           = cur_x,
        y           = cur_y,
        node_id     = cur_node_id,
        tree        = tree,
        max_pop_pct = max_pop_pct,
        nsim        = nsim, alpha = alpha,
        model       = model,
        seed        = iter_seed, n_cores = n_cores
      )
    } else if (scan_type == "circular") {
      res <- circular_scan(
        cases       = cur_cases,
        population  = cur_population,
        region_id   = cur_region_id,
        x           = cur_x,
        y           = cur_y,
        max_pop_pct = max_pop_pct,
        nsim        = nsim, alpha = alpha,
        model       = model,
        seed        = iter_seed, n_cores = n_cores
      )
    } else {
      res <- tree_scan(
        tree       = tree,
        cases      = cur_cases,
        population = cur_population,
        nsim       = nsim, alpha = alpha,
        model      = model,
        seed       = iter_seed, n_cores = n_cores
      )
    }

    mlc <- res$most_likely_cluster
    pvalue <- res$pvalue
    is_sig <- !is.na(pvalue) && pvalue <= alpha

    # Stop on dead signal (LLR=0 or no cases): nothing left to find
    if (is.null(mlc$llr) || mlc$llr <= 0 ||
        is.null(mlc$cases) || mlc$cases == 0) {
      if (verbose) message(sprintf(
        "  Iteration %d: no residual signal (LR = 0). Stopping.", k))
      break
    }

    # Record this cluster
    if (scan_type == "tree") {
      row <- data.frame(
        iteration  = k,
        node_id    = as.character(mlc$node_id),
        cases      = mlc$cases,
        expected   = mlc$expected,
        population = mlc$population,
        llr        = mlc$llr,
        pvalue     = pvalue,
        significant = is_sig,
        stringsAsFactors = FALSE
      )
      row$leaf_ids <- I(list(mlc$leaf_ids))
    } else {
      row <- data.frame(
        iteration  = k,
        node_id    = if (scan_type == "treespatial")
                       as.character(mlc$node_id) else NA_character_,
        n_regions  = length(mlc$region_ids),
        cases      = mlc$cases,
        expected   = mlc$expected,
        population = mlc$population,
        rr         = mlc$rr,
        llr        = mlc$llr,
        pvalue     = pvalue,
        significant = is_sig,
        stringsAsFactors = FALSE
      )
      row$region_ids <- I(list(mlc$region_ids))
    }

    clusters_rows[[k]] <- row
    iterations[[k]]    <- res

    if (verbose) {
      if (scan_type == "tree") {
        message(sprintf(
          "  Cluster %d: node = %s, LR = %.2f, p = %.3f%s",
          k, mlc$node_id, mlc$llr, pvalue,
          if (is_sig) " *" else ""))
      } else if (scan_type == "treespatial") {
        message(sprintf(
          "  Cluster %d: node = %s, %d region(s), LR = %.2f, p = %.3f%s",
          k, mlc$node_id, length(mlc$region_ids), mlc$llr, pvalue,
          if (is_sig) " *" else ""))
      } else {
        message(sprintf(
          "  Cluster %d: %d region(s), LR = %.2f, p = %.3f%s",
          k, length(mlc$region_ids), mlc$llr, pvalue,
          if (is_sig) " *" else ""))
      }
    }

    # Stop after recording if not significant -- the sequential procedure
    # of Zhang et al. (2010) tests in order, stopping at the first
    # non-significant MLC.
    if (!is_sig) {
      if (verbose) message(sprintf(
        "  Iteration %d: MLC is not significant at alpha = %g. Stopping.",
        k, alpha))
      break
    }

    # Stop if budget reached
    if (k == max_iter) break

    # --- Remove the MLC regions (+ optional buffer) from working data ---
    if (scan_type == "tree") {
      # Remove the leaves descending from the detected node from the
      # current data.
      leaves_all <- .get_leaves(tree)
      desc <- .get_descendants(tree, as.character(mlc$node_id))
      desc_leaves <- intersect(desc, leaves_all)
      # cur_cases / cur_population are vectors aligned to leaves_all
      keep <- !(leaves_all %in% desc_leaves)
      cur_cases      <- cur_cases[keep]
      cur_population <- cur_population[keep]
      # Now the next tree_scan needs a tree consistent with the kept
      # leaves: prune the tree to the kept leaves and their ancestors.
      removed_leaf_ids <- union(removed_leaf_ids, desc_leaves)
      tree <- .prune_tree_to_leaves(tree, leaves_all[keep])
    } else {
      regs_to_remove <- mlc$region_ids
      if (buffer_size > 0L) {
        regs_to_remove <- .add_buffer_regions(
          cluster_region_ids = regs_to_remove,
          all_region_id      = cur_region_id,
          all_x              = cur_x,
          all_y              = cur_y,
          buffer_size        = buffer_size
        )
      }
      removed_region_ids <- union(removed_region_ids, regs_to_remove)
      keep <- !(cur_region_id %in% regs_to_remove)
      cur_cases      <- cur_cases[keep]
      cur_population <- cur_population[keep]
      cur_region_id  <- cur_region_id[keep]
      cur_x          <- cur_x[keep]
      cur_y          <- cur_y[keep]
      if (scan_type == "treespatial") {
        cur_node_id <- cur_node_id[keep]
      }
    }
  }

  # Build the regions table that we attach to the result, for plotting
  # downstream. We use the ORIGINAL (uncut) inputs so that
  # get_cluster_regions() can colour every region.
  if (scan_type == "tree") {
    out_regions <- NULL
  } else if (scan_type == "circular") {
    # Circular inputs are already one row per region.
    out_regions <- .build_regions_circular(
      cases = cases, population = population,
      region_id = region_id, x = x, y = y
    )
  } else {
    # Tree-spatial: inputs are long (region x leaf). Build the unique
    # regions table directly, taking population/x/y from the first
    # occurrence of each region_id (matches .build_inputs).
    reg_ids <- unique(region_id)
    first_idx <- match(reg_ids, region_id)
    out_regions <- data.frame(
      region_id  = reg_ids,
      population = population[first_idx],
      x          = x[first_idx],
      y          = y[first_idx],
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    clusters    = if (length(clusters_rows) > 0L)
                    do.call(rbind, clusters_rows) else data.frame(),
    iterations  = iterations,
    regions     = out_regions,
    tree        = tree,
    alpha       = alpha,
    nsim        = nsim,
    buffer_size = buffer_size,
    n_iter      = length(clusters_rows),
    scan_type   = scan_type
  )
  class(out) <- "sequential_scan"
  out
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
#' Compute the union of cluster regions and their nearest-neighbour buffer.
#' Distance is Euclidean on (x, y); regions are picked by smallest minimum
#' distance to any cluster region.
.add_buffer_regions <- function(cluster_region_ids, all_region_id,
                                 all_x, all_y, buffer_size) {

  cl_mask  <- all_region_id %in% cluster_region_ids
  if (!any(cl_mask) || all(cl_mask) || buffer_size == 0L) {
    return(unique(cluster_region_ids))
  }

  out_idx <- which(!cl_mask)
  cl_x <- all_x[cl_mask]; cl_y <- all_y[cl_mask]

  # Min distance from each outside region to any cluster region.
  min_dist <- vapply(out_idx, function(i) {
    d <- sqrt((cl_x - all_x[i])^2 + (cl_y - all_y[i])^2)
    min(d)
  }, numeric(1L))

  k <- min(buffer_size, length(out_idx))
  picked <- out_idx[order(min_dist)[seq_len(k)]]
  unique(c(cluster_region_ids, all_region_id[picked]))
}


#' @keywords internal
#' Prune a tree to a subset of leaves: keep those leaves and all their
#' ancestors; drop any internal node that has no kept descendant.
.prune_tree_to_leaves <- function(tree, leaves_kept) {
  if (length(leaves_kept) == 0L) {
    return(tree[integer(0), , drop = FALSE])
  }
  # Ancestors of each kept leaf
  needed <- character(0)
  for (lf in leaves_kept) {
    nd <- as.character(lf)
    while (length(nd) == 1L && !is.na(nd) && nd != "") {
      needed <- c(needed, nd)
      parent <- tree$parent_id[match(nd, tree$node_id)]
      if (is.na(parent)) break
      nd <- as.character(parent)
      if (nd %in% needed) break
    }
  }
  needed <- unique(needed)
  tree[tree$node_id %in% needed, , drop = FALSE]
}
