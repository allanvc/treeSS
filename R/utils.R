#' @keywords internal
.validate_regions <- function(regions) {
  required <- c("region_id", "population", "x", "y")
  missing_cols <- setdiff(required, names(regions))
  if (length(missing_cols) > 0) {
    stop("'regions' is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Required: region_id, population, x, y",
         call. = FALSE)
  }
  if (any(regions$population < 0)) {
    stop("'population' must be non-negative.", call. = FALSE)
  }
  if (nrow(regions) < 2) {
    stop("At least 2 regions are required.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.validate_tree <- function(tree) {
  required <- c("node_id", "parent_id")
  missing_cols <- setdiff(required, names(tree))
  if (length(missing_cols) > 0) {
    stop("'tree' is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Required: node_id, parent_id",
         call. = FALSE)
  }
  # Check for a single root (parent_id is NA)
  roots <- is.na(tree$parent_id)
  if (sum(roots) == 0) {
    stop("Tree must have at least one root node (parent_id = NA).",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.validate_cases_matrix <- function(cases, tree, regions) {
  if (!is.matrix(cases) && !is.data.frame(cases)) {
    stop("'cases' must be a matrix or data.frame.", call. = FALSE)
  }
  cases <- as.matrix(cases)
  leaves <- .get_leaves(tree)
  if (nrow(cases) != length(leaves)) {
    stop("Number of rows in 'cases' (", nrow(cases), ") must match the ",
         "number of leaf nodes in the tree (", length(leaves), ").",
         call. = FALSE)
  }
  if (ncol(cases) != nrow(regions)) {
    stop("Number of columns in 'cases' (", ncol(cases), ") must match the ",
         "number of regions (", nrow(regions), ").",
         call. = FALSE)
  }
  if (any(cases < 0, na.rm = TRUE)) {
    stop("'cases' must contain non-negative values.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.get_leaves <- function(tree) {
  all_nodes <- tree$node_id
  parent_nodes <- unique(tree$parent_id[!is.na(tree$parent_id)])
  leaves <- setdiff(all_nodes, parent_nodes)
  leaves
}

#' @keywords internal
.get_branches <- function(tree) {
  # Returns a list where each element is a branch (node + all its descendants)
  all_nodes <- tree$node_id
  branches <- list()

  for (node in all_nodes) {
    descendants <- .get_descendants(tree, node)
    branches[[as.character(node)]] <- list(
      node = node,
      leaves = intersect(descendants, .get_leaves(tree))
    )
  }

  branches
}

#' @keywords internal
.get_descendants <- function(tree, node) {
  children <- tree$node_id[!is.na(tree$parent_id) & tree$parent_id == node]
  if (length(children) == 0) {
    return(node)
  }
  desc <- node
  for (child in children) {
    desc <- c(desc, .get_descendants(tree, child))
  }
  desc
}

#' @keywords internal
.poisson_llr <- function(cz, nz, C, N) {
  # Log-likelihood ratio for the Poisson model
  # cz: cases inside zone, nz: population inside zone
  # C: total cases, N: total population
  cz_bar <- C - cz
  nz_bar <- N - nz

  expected <- C * nz / N

  if (cz <= expected || cz == 0 || nz == 0 || nz_bar == 0) {
    return(0)
  }

  llr <- 0
  if (cz > 0) {
    llr <- llr + cz * log(cz / expected)
  }
  if (cz_bar > 0) {
    llr <- llr + cz_bar * log(cz_bar / (C - expected))
  }

  llr
}


#' @keywords internal
#' Bottom-up aggregation: given a leaf x region matrix, return a
#' (all_nodes) x region matrix where each non-leaf row is the sum of its
#' children rows.
.aggregate_leaves_to_all <- function(cases_matrix, tree) {
  leaves    <- .get_leaves(tree)
  all_nodes <- tree$node_id
  n_nodes   <- length(all_nodes)
  n_regions <- ncol(cases_matrix)

  full_cases <- matrix(0, nrow = n_nodes, ncol = n_regions)
  rownames(full_cases) <- as.character(all_nodes)
  colnames(full_cases) <- colnames(cases_matrix)

  leaf_idx <- match(leaves, all_nodes)
  full_cases[leaf_idx, ] <- cases_matrix

  depths <- .compute_depths(tree)
  processing_order <- order(depths, decreasing = TRUE)

  for (idx in processing_order) {
    node <- all_nodes[idx]
    children_idx <- which(tree$parent_id == node & !is.na(tree$parent_id))
    if (length(children_idx) > 0) {
      child_rows <- match(tree$node_id[children_idx], all_nodes)
      if (length(child_rows) == 1) {
        full_cases[idx, ] <- full_cases[child_rows, ]
      } else {
        full_cases[idx, ] <- colSums(full_cases[child_rows, , drop = FALSE])
      }
    }
  }
  full_cases
}


#' @keywords internal
#' Compute the depth of each node in the tree (root = 0).
.compute_depths <- function(tree) {
  depths <- rep(NA_integer_, nrow(tree))
  roots <- which(is.na(tree$parent_id))
  depths[roots] <- 0L

  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (i in seq_len(nrow(tree))) {
      if (is.na(depths[i]) && !is.na(tree$parent_id[i])) {
        parent_idx <- match(tree$parent_id[i], tree$node_id)
        if (!is.na(parent_idx) && !is.na(depths[parent_idx])) {
          depths[i] <- depths[parent_idx] + 1L
          changed <- TRUE
        }
      }
    }
  }
  depths
}
