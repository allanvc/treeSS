#' Generate Example Data for Tree-Spatial Scan
#'
#' Creates a synthetic dataset for demonstrating and testing the tree-spatial
#' scan statistic. Returns parallel vectors (\code{cases}, \code{population},
#' \code{region_id}, \code{x}, \code{y}, \code{node_id}) and a \code{tree},
#' matching the input format expected by \code{\link{treespatial_scan}}.
#'
#' @param n_regions Integer. Default 50.
#' @param pop_per_region Numeric. Default 1000.
#' @param cluster_regions Integer vector. Default \code{1:7}.
#' @param cluster_leaves Integer vector. Default \code{c(3, 4)}.
#' @param rr Numeric. Relative risk. Default 2.0.
#' @param Cg Integer. Cases per branch. Default 200.
#' @param seed Integer. Random seed.
#'
#' @return A list with vector components ready to feed into
#'   \code{\link{treespatial_scan}}: \code{cases}, \code{population},
#'   \code{region_id}, \code{x}, \code{y}, \code{node_id}, plus the
#'   \code{tree} (data.frame) and a \code{true_cluster} list describing the
#'   injected cluster.
#'
#' @export
#' @examples
#' ex <- generate_example_data(seed = 42)
#' result <- treespatial_scan(
#'   cases       = ex$cases,
#'   population  = ex$population,
#'   region_id   = ex$region_id,
#'   x           = ex$x,
#'   y           = ex$y,
#'   node_id     = ex$node_id,
#'   tree        = ex$tree,
#'   nsim        = 99
#' )
#' print(result)
generate_example_data <- function(n_regions = 50L, pop_per_region = 1000,
                                  cluster_regions = 1:7,
                                  cluster_leaves = c(3, 4),
                                  rr = 2.0, Cg = 200L, seed = 123L) {
  set.seed(seed)

  sq <- ceiling(sqrt(n_regions))
  grid_x <- rep(seq_len(sq), each = sq)[seq_len(n_regions)]
  grid_y <- rep(seq_len(sq), times = sq)[seq_len(n_regions)]
  grid_x <- grid_x + 0.5 * (grid_y %% 2)

  region_ids <- seq_len(n_regions)
  pop_vec    <- rep(pop_per_region, n_regions)

  tree <- data.frame(
    node_id   = 1:15,
    parent_id = c(NA, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
  )
  n_leaves <- 8
  leaf_node_ids <- 8:15

  N <- n_regions * pop_per_region
  prob_base <- pop_vec / N

  cases_mat <- matrix(0L, nrow = n_leaves, ncol = n_regions)

  for (leaf in seq_len(n_leaves)) {
    if (leaf %in% cluster_leaves) {
      prob_cluster <- prob_base
      prob_cluster[cluster_regions] <- prob_cluster[cluster_regions] * rr
      prob_cluster <- prob_cluster / sum(prob_cluster)
      cases_mat[leaf, ] <- as.integer(stats::rmultinom(1, Cg, prob_cluster))
    } else {
      cases_mat[leaf, ] <- as.integer(stats::rmultinom(1, Cg, prob_base))
    }
  }

  # Convert to long parallel vectors
  pos <- which(cases_mat > 0, arr.ind = TRUE)
  cases_vec   <- cases_mat[pos]
  region_vec  <- region_ids[pos[, "col"]]
  node_vec    <- leaf_node_ids[pos[, "row"]]
  pop_long    <- pop_vec[pos[, "col"]]
  x_long      <- as.numeric(grid_x)[pos[, "col"]]
  y_long      <- as.numeric(grid_y)[pos[, "col"]]

  true_leaf_nodes <- leaf_node_ids[cluster_leaves]
  true_cluster <- list(
    region_ids    = cluster_regions,
    leaf_indices  = cluster_leaves,
    leaf_node_ids = true_leaf_nodes,
    rr            = rr
  )

  list(
    cases        = as.integer(cases_vec),
    population   = pop_long,
    region_id    = as.integer(region_vec),
    x            = x_long,
    y            = y_long,
    node_id      = as.integer(node_vec),
    tree         = tree,
    true_cluster = true_cluster
  )
}
