#' Aggregate Case Counts from Leaves to All Nodes in a Hierarchical Tree
#'
#' Given case counts at the leaf level (as parallel vectors) and a tree,
#' aggregates case counts upward from child nodes to parent nodes,
#' producing a complete case matrix indexed by all nodes (rows) and
#' regions (columns).
#'
#' @param cases Numeric vector of length \eqn{n}: case counts.
#' @param region_id Vector of region identifiers, length \eqn{n}.
#' @param node_id Vector of leaf identifiers, length \eqn{n}.
#' @param tree A \code{data.frame} with columns \code{node_id} and
#'   \code{parent_id}. As an alternative, pass \code{tree_node_id} and
#'   \code{tree_parent_id}.
#' @param tree_node_id,tree_parent_id Optional. Parallel vectors as an
#'   alternative to \code{tree}.
#'
#' @return A matrix of dimensions \eqn{m \times k} (nodes x regions), with
#'   rows ordered by \code{tree$node_id} and columns by region.
#'
#' @details
#' This function is exposed for inspection and pedagogical use; the scan
#' functions call it internally on the matrix they build from your input
#' vectors.
#'
#' @references
#' Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
#' (2025). A tree-spatial scan statistic.
#' \emph{Environmental and Ecological Statistics}, 32, 953-978.
#' \doi{10.1007/s10651-025-00670-w}
#'
#' @export
#' @examples
#' tree <- data.frame(
#'   node_id   = c(1, 2, 3, 4, 5, 6, 7),
#'   parent_id = c(NA, 1, 1, 2, 2, 3, 3)
#' )
#' # Leaves are 4, 5, 6, 7
#' aggregate_tree(
#'   cases     = c(10, 5, 3, 8,  2, 7, 4, 1,  6, 3, 9, 2),
#'   region_id = rep(1:3, each = 4),
#'   node_id   = rep(c(4, 5, 6, 7), times = 3),
#'   tree      = tree
#' )
aggregate_tree <- function(cases, region_id, node_id,
                           tree = NULL,
                           tree_node_id = NULL, tree_parent_id = NULL) {

  tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)
  .validate_tree(tree)

  if (length(cases) != length(region_id) ||
      length(cases) != length(node_id)) {
    stop("'cases', 'region_id', 'node_id' must have the same length.",
         call. = FALSE)
  }
  if (any(is.na(cases))) stop("'cases' contains NA values.", call. = FALSE)
  if (any(cases < 0))    stop("'cases' must be non-negative.", call. = FALSE)

  leaves <- .get_leaves(tree)
  reg_ids <- unique(region_id)
  n_leaves <- length(leaves)
  n_regions <- length(reg_ids)

  tree_leaves <- leaves
  if (is.character(tree$node_id) || is.character(node_id)) {
    node_id     <- as.character(node_id)
    tree_leaves <- as.character(tree_leaves)
  }
  reg_ids_chr_check <- is.character(reg_ids) || is.character(region_id)
  if (reg_ids_chr_check) {
    region_id <- as.character(region_id)
    reg_ids   <- as.character(reg_ids)
  }

  unknown_nodes <- setdiff(unique(node_id), tree_leaves)
  if (length(unknown_nodes) > 0) {
    stop("'node_id' contains values that are not leaves of the tree: ",
         paste(utils::head(unknown_nodes, 5), collapse = ", "),
         if (length(unknown_nodes) > 5) ", ..." else "",
         ".", call. = FALSE)
  }

  row_idx <- match(node_id, tree_leaves)
  col_idx <- match(region_id, reg_ids)

  cases_mat <- matrix(0, nrow = n_leaves, ncol = n_regions)
  rownames(cases_mat) <- as.character(tree_leaves)
  colnames(cases_mat) <- as.character(reg_ids)

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

  .aggregate_leaves_to_all(cases_mat, tree)
}
