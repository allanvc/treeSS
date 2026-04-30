## =============================================================================
## vector_inputs.R
##
## Helpers introduced in treeSS 0.1.5 for the vector-based API.
##
## The user supplies parallel vectors (cases, population, region_id, x, y,
## node_id) and a tree (data.frame or two vectors). These helpers:
##   - normalize the tree input to a 2-column data.frame
##   - validate that the parallel vectors are coherent
##   - build the regions table and the leaf x region cases matrix that
##     the C++ core consumes
## =============================================================================


#' @keywords internal
#' Coerce a tree input (either a 2-column data.frame, or two parallel
#' vectors) into a canonical data.frame with columns node_id, parent_id.
.normalize_tree <- function(tree = NULL,
                             tree_node_id = NULL, tree_parent_id = NULL) {
  has_df  <- !is.null(tree)
  has_vec <- !is.null(tree_node_id) && !is.null(tree_parent_id)

  if (has_df && has_vec) {
    stop("Provide either 'tree' (data.frame) or both 'tree_node_id' and ",
         "'tree_parent_id' (vectors), not both.", call. = FALSE)
  }
  if (!has_df && !has_vec) {
    stop("Tree input is missing. Provide 'tree' as a data.frame with ",
         "columns 'node_id' and 'parent_id', or pass 'tree_node_id' and ",
         "'tree_parent_id' as parallel vectors.", call. = FALSE)
  }

  if (has_vec) {
    if (length(tree_node_id) != length(tree_parent_id)) {
      stop("'tree_node_id' and 'tree_parent_id' must have the same length.",
           call. = FALSE)
    }
    tree <- data.frame(node_id   = tree_node_id,
                       parent_id = tree_parent_id,
                       stringsAsFactors = FALSE)
  } else {
    if (!is.data.frame(tree)) {
      stop("'tree' must be a data.frame with columns 'node_id' and ",
           "'parent_id'.", call. = FALSE)
    }
    if (!all(c("node_id", "parent_id") %in% names(tree))) {
      stop("'tree' must contain columns 'node_id' and 'parent_id'.",
           call. = FALSE)
    }
  }
  tree
}


#' @keywords internal
#' Validate that the parallel vectors (cases, population, region_id, x, y,
#' node_id) are coherent: same length, no NAs in keys, etc.
.validate_vectors <- function(cases, population, region_id, x, y,
                              node_id = NULL) {
  n <- length(cases)
  if (length(population) != n) {
    stop("'population' must have the same length as 'cases'.", call. = FALSE)
  }
  if (length(region_id) != n) {
    stop("'region_id' must have the same length as 'cases'.", call. = FALSE)
  }
  if (length(x) != n) {
    stop("'x' must have the same length as 'cases'.", call. = FALSE)
  }
  if (length(y) != n) {
    stop("'y' must have the same length as 'cases'.", call. = FALSE)
  }
  if (!is.null(node_id) && length(node_id) != n) {
    stop("'node_id' must have the same length as 'cases'.", call. = FALSE)
  }
  if (any(is.na(cases))) stop("'cases' contains NA values.", call. = FALSE)
  if (any(cases < 0))    stop("'cases' must be non-negative.", call. = FALSE)
  if (any(is.na(region_id))) stop("'region_id' contains NA values.",
                                  call. = FALSE)
  if (!is.null(node_id) && any(is.na(node_id))) {
    stop("'node_id' contains NA values.", call. = FALSE)
  }
  invisible(TRUE)
}


#' @keywords internal
#' From parallel vectors, build:
#'   - regions: a data.frame with one row per unique region_id
#'              (with population, x, y; region_id is preserved)
#'   - cases_matrix: leaf x region matrix (rows in tree leaf order,
#'                   columns in regions$region_id order)
#'
#' Population/x/y are taken as the FIRST occurrence per region. If they
#' vary across rows of the same region, a warning is issued.
.build_inputs <- function(cases, population, region_id, x, y,
                           node_id, tree) {

  .validate_vectors(cases, population, region_id, x, y, node_id)
  .validate_tree(tree)

  # Order by region_id (first occurrence) to make regions table deterministic
  reg_ids <- unique(region_id)

  # Take population/x/y from first occurrence; warn if inconsistent
  first_idx <- match(reg_ids, region_id)
  reg_pop <- population[first_idx]
  reg_x   <- x[first_idx]
  reg_y   <- y[first_idx]

  # Sanity check: warn if population/x/y vary within region
  for (r in reg_ids) {
    rows <- which(region_id == r)
    if (length(rows) > 1) {
      if (length(unique(population[rows])) > 1L) {
        warning("'population' varies within region_id == ", r,
                "; using the first value.", call. = FALSE)
      }
      if (length(unique(x[rows])) > 1L || length(unique(y[rows])) > 1L) {
        warning("'x' or 'y' varies within region_id == ", r,
                "; using the first values.", call. = FALSE)
      }
    }
  }

  regions <- data.frame(
    region_id  = reg_ids,
    population = reg_pop,
    x          = reg_x,
    y          = reg_y,
    stringsAsFactors = FALSE
  )

  # Build cases matrix
  leaves <- .get_leaves(tree)
  n_leaves  <- length(leaves)
  n_regions <- nrow(regions)

  tree_leaves <- leaves
  if (is.character(tree$node_id) || is.character(node_id)) {
    node_id     <- as.character(node_id)
    tree_leaves <- as.character(tree_leaves)
  }
  if (is.character(regions$region_id) || is.character(region_id)) {
    region_id  <- as.character(region_id)
    regions$region_id <- as.character(regions$region_id)
  }

  unknown_nodes <- setdiff(unique(node_id), tree_leaves)
  if (length(unknown_nodes) > 0) {
    stop("'node_id' contains values that are not leaves of the tree: ",
         paste(utils::head(unknown_nodes, 5), collapse = ", "),
         if (length(unknown_nodes) > 5) ", ..." else "",
         ".", call. = FALSE)
  }

  row_idx <- match(node_id, tree_leaves)
  col_idx <- match(region_id, regions$region_id)

  cases_mat <- matrix(0, nrow = n_leaves, ncol = n_regions)
  rownames(cases_mat) <- as.character(tree_leaves)
  colnames(cases_mat) <- as.character(regions$region_id)

  # Drop rows where cases == 0 to match earlier behaviour
  keep <- cases > 0
  if (any(keep)) {
    key <- paste(row_idx[keep], col_idx[keep], sep = "_")
    if (anyDuplicated(key)) {
      agg <- stats::aggregate(cases[keep],
                              by = list(r = row_idx[keep],
                                        c = col_idx[keep]),
                              FUN = sum)
      cases_mat[cbind(agg$r, agg$c)] <- agg$x
    } else {
      cases_mat[cbind(row_idx[keep], col_idx[keep])] <- cases[keep]
    }
  }

  list(regions = regions, cases_matrix = cases_mat)
}


#' @keywords internal
#' Build a regions table from circular-scan vector inputs (no node_id).
.build_regions_circular <- function(cases, population, region_id, x, y) {
  .validate_vectors(cases, population, region_id, x, y, node_id = NULL)
  n <- length(cases)
  reg_ids <- unique(region_id)
  if (length(reg_ids) != n) {
    stop("For circular_scan, 'cases' must already be aggregated to one ",
         "value per region. Got ", n, " observations but ",
         length(reg_ids), " unique region_id values.\n",
         "If your data is in long format, aggregate first, e.g.:\n",
         "  agg <- aggregate(cases ~ region_id + x + y + population, ",
         "data = ..., FUN = sum)\n",
         "and then pass agg$cases, agg$population, agg$region_id, agg$x, ",
         "agg$y to circular_scan().", call. = FALSE)
  }
  data.frame(
    region_id  = region_id,
    population = population,
    x          = x,
    y          = y,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
#' Inverse of .build_inputs: take a leaf x region cases matrix + regions
#' table + tree, return parallel vectors suitable for treespatial_scan.
#' Used by iterative_scan to feed the next iteration after zero-ing a
#' cluster.
.matrix_to_vectors <- function(cases_matrix, regions, tree) {
  leaves <- .get_leaves(tree)
  n_leaves <- length(leaves)
  n_regions <- nrow(regions)

  pos <- which(cases_matrix > 0, arr.ind = TRUE)
  if (length(pos) == 0) {
    # No cases anywhere - return one dummy zero row per region/leaf so
    # downstream validation doesn't break
    out <- list(
      cases      = rep(0, n_regions),
      population = regions$population,
      region_id  = regions$region_id,
      x          = regions$x,
      y          = regions$y,
      node_id    = rep(leaves[1], n_regions)
    )
    return(out)
  }

  rid_idx <- pos[, "col"]
  nid_idx <- pos[, "row"]

  list(
    cases      = cases_matrix[pos],
    population = regions$population[rid_idx],
    region_id  = regions$region_id[rid_idx],
    x          = regions$x[rid_idx],
    y          = regions$y[rid_idx],
    node_id    = leaves[nid_idx]
  )
}
