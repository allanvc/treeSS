#' Tree-Spatial Scan Statistic
#'
#' Performs the tree-spatial scan statistic, combining Kulldorff's circular
#' scan (spatial clusters) with the tree-based scan (hierarchical data
#' mining). Searches all combinations of spatial zones and tree branches
#' to identify pairs \eqn{(z, g)} with significantly more cases than
#' expected.
#'
#' Inputs are passed as parallel vectors of equal length (one entry per
#' (region, tree-leaf) observation). The user is responsible for choosing
#' which column to use as \code{population} (e.g., total population, live
#' births, person-years), making the choice of denominator explicit.
#'
#' \strong{Secondary clusters.} The returned object contains the most likely
#' cluster as well as the full set of evaluated (zone, branch) pairs in
#' \code{secondary_clusters}. To obtain the distinct secondary clusters as
#' described in Section 5.1.1 of Cancado et al. (2025) (filtering out pairs
#' that overlap in regions or branches with already-retained clusters), use
#' \code{\link{filter_clusters}} or \code{\link{get_cluster_regions}} with
#' \code{n_clusters > 1}.
#'
#' @param cases Numeric vector. Number of cases observed for each
#'   (region, leaf) pair. Length \eqn{n}.
#' @param population Numeric vector. Population (or denominator) of the
#'   region for each row. The same value should be repeated across all
#'   rows of a given region; if it varies, the first occurrence per
#'   region is used and a warning is issued.
#' @param region_id Vector of region identifiers. Length \eqn{n}.
#' @param x,y Numeric vectors of region centroid coordinates. Like
#'   \code{population}, these should be constant within region.
#' @param node_id Vector of tree leaf identifiers. Length \eqn{n}. Each
#'   value must match a leaf of the tree.
#' @param tree A \code{data.frame} with columns \code{node_id} and
#'   \code{parent_id}. The root node(s) must have \code{parent_id = NA}.
#'   As an alternative, pass \code{tree_node_id} and \code{tree_parent_id}
#'   as parallel vectors instead of this argument.
#' @param tree_node_id,tree_parent_id Optional. Parallel vectors describing
#'   the tree edges, used as an alternative to \code{tree}. If both
#'   \code{tree} and these vectors are supplied, an error is raised.
#' @param max_pop_pct Numeric. Maximum proportion of total population
#'   allowed inside a zone. Default \code{0.5}.
#' @param nsim Integer. Number of Monte Carlo simulations. Default \code{999}.
#' @param alpha Numeric. Significance level. Default \code{0.05}.
#' @param model Character. \code{"poisson"} (default) or \code{"binomial"}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#' @param n_cores Integer. OpenMP threads for the Monte Carlo loop.
#'   Default \code{1L} (serial).
#'
#' @return An object of class \code{"treespatial_scan"}.
#'
#' @references
#' Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
#' (2025). A tree-spatial scan statistic. \emph{Environmental and Ecological
#' Statistics}, 32, 953--978. \doi{10.1007/s10651-025-00670-w}
#'
#' @seealso \code{\link{circular_scan}}, \code{\link{tree_scan}},
#'   \code{\link{aggregate_tree}}, \code{\link{filter_clusters}},
#'   \code{\link{get_cluster_regions}}, \code{\link{iterative_scan}}
#'
#' @export
#' @examples
#' set.seed(123)
#' n_regions <- 10
#' tree <- data.frame(
#'   node_id   = c(1, 2, 3, 4, 5, 6, 7),
#'   parent_id = c(NA, 1, 1, 2, 2, 3, 3)
#' )
#' # Build vectors: one row per (region, leaf) combination
#' grid <- expand.grid(region_id = 1:n_regions, node_id = c(4, 5, 6, 7))
#' xs   <- runif(n_regions, 0, 10)[grid$region_id]
#' ys   <- runif(n_regions, 0, 10)[grid$region_id]
#' cs   <- rpois(nrow(grid), lambda = 5)
#' cs[grid$node_id == 4 & grid$region_id %in% 1:3] <- rpois(3, 30)
#'
#' result <- treespatial_scan(
#'   cases       = cs,
#'   population  = rep(1000, nrow(grid)),
#'   region_id   = grid$region_id,
#'   x           = xs,
#'   y           = ys,
#'   node_id     = grid$node_id,
#'   tree        = tree,
#'   nsim        = 99
#' )
#' print(result)
treespatial_scan <- function(cases, population, region_id, x, y, node_id,
                             tree = NULL,
                             tree_node_id = NULL, tree_parent_id = NULL,
                             max_pop_pct = 0.5,
                             nsim = 999L, alpha = 0.05,
                             model = c("poisson", "binomial"),
                             seed = NULL, n_cores = 1L) {

  model <- match.arg(model)
  model_int <- if (model == "binomial") 1L else 0L

  tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)

  # Build internal regions table + cases matrix from parallel vectors
  prepared <- .build_inputs(cases, population, region_id, x, y, node_id, tree)
  regions      <- prepared$regions
  cases_matrix <- prepared$cases_matrix

  .validate_regions(regions)

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(regions)
  N <- sum(regions$population)
  max_pop <- max_pop_pct * N
  leaves <- .get_leaves(tree)
  all_nodes <- tree$node_id
  n_nodes <- length(all_nodes)

  # Aggregate cases from leaves to all nodes (uses the cases_matrix built
  # above, not the user's cases vector)
  full_cases <- .aggregate_leaves_to_all(cases_matrix, tree)

  # Total cases per branch (across all regions)
  Cg <- rowSums(full_cases)
  names(Cg) <- as.character(all_nodes)

  # Build candidate spatial zones
  zones <- build_zones(regions, max_pop = max_pop)
  n_zones <- length(zones)

  if (n_zones == 0) {
    stop("No valid candidate zones found.", call. = FALSE)
  }

  # --- C++ backend: observed scan + Monte Carlo ---
  prob <- regions$population / N

  # Prepare CSR structures
  zone_csr <- .zones_to_csr(zones)
  tree_csr <- .tree_to_csr_children(tree)
  leaf_idx_0 <- as.integer(match(leaves, all_nodes) - 1L)  # 0-based
  depths <- .compute_depths(tree)
  proc_order_0 <- as.integer(order(depths, decreasing = TRUE) - 1L)  # 0-based

  cpp_result <- mc_treespatial_cpp(
    as.matrix(full_cases),
    as.numeric(Cg),
    as.numeric(N),
    zone_csr$zone_region_idx,
    zone_csr$zone_ptr,
    zone_csr$zone_pop,
    leaf_idx_0,
    tree_csr$children_idx,
    tree_csr$children_ptr,
    proc_order_0,
    as.numeric(prob),
    as.integer(nsim),
    model_int,
    as.integer(max(1L, n_cores))
  )

  obs_max_llr <- cpp_result$obs_max_llr
  best_g <- cpp_result$obs_best_g  # 1-based
  best_z <- cpp_result$obs_best_z  # 1-based
  sim_llr <- cpp_result$sim_llr

  # P-value
  pvalue <- (sum(sim_llr >= obs_max_llr) + 1) / (nsim + 1)

  # --- Secondary clusters (R-side, single pass, using matrix ops) ---
  zone_mat <- matrix(0L, nrow = n_zones, ncol = n)
  zone_pop_vec <- numeric(n_zones)
  zone_centers <- integer(n_zones)
  zone_sizes <- integer(n_zones)
  zone_region_list <- vector("list", n_zones)

  for (zi in seq_len(n_zones)) {
    idx <- zones[[zi]]$region_idx
    zone_mat[zi, idx] <- 1L
    zone_pop_vec[zi] <- zones[[zi]]$population
    zone_centers[zi] <- zones[[zi]]$center
    zone_sizes[zi] <- length(idx)
    zone_region_list[[zi]] <- idx
  }

  all_pairs <- .build_pairs_df_r(full_cases, Cg, N, zone_mat, zone_pop_vec,
                                  zone_centers, zone_sizes, all_nodes,
                                  model = model)

  # --- Build result ---
  branches <- .get_branches(tree)
  best_leaves <- branches[[as.character(all_nodes[best_g])]]$leaves
  best_zone_idx <- zone_region_list[[best_z]]
  cases_in_zone <- sum(full_cases[best_g, best_zone_idx])
  expected_in_zone <- Cg[best_g] * zone_pop_vec[best_z] / N

  result <- list(
    most_likely_cluster = list(
      node_id = all_nodes[best_g],
      leaf_ids = best_leaves,
      zone_idx = best_zone_idx,
      region_ids = regions$region_id[best_zone_idx],
      cases = cases_in_zone,
      expected = expected_in_zone,
      population = zone_pop_vec[best_z],
      rr = (cases_in_zone / zone_pop_vec[best_z]) / (Cg[best_g] / N),
      llr = obs_max_llr
    ),
    secondary_clusters = all_pairs,
    pvalue = pvalue,
    alpha = alpha,
    nsim = nsim,
    total_cases_by_branch = Cg,
    total_population = N,
    regions = regions,
    tree = tree,
    full_cases = full_cases,
    simulated_llr = sim_llr,
    model = model
  )

  class(result) <- "treespatial_scan"
  result
}

#' @keywords internal
.build_pairs_df_r <- function(full_cases, Cg, N, zone_mat, zone_pop,
                               zone_centers, zone_sizes, all_nodes,
                               model = "poisson") {
  # Compute all (node, zone) cases via single matrix multiply
  cz_mat <- full_cases %*% t(zone_mat)
  n_nodes <- nrow(cz_mat)
  n_zones <- ncol(cz_mat)

  nz_mat <- matrix(rep(zone_pop, each = n_nodes), nrow = n_nodes)
  nz_bar_mat <- N - nz_mat
  cz_bar <- Cg - cz_mat

  # Always compute "expected under H0" - same in Poisson and Binomial as a
  # descriptive statistic (Cg * zone_pop / N). It enters the LLR formula
  # only in the Poisson case, but we report it on the data.frame either way.
  expected <- outer(Cg, zone_pop) / N

  if (model == "binomial") {
    # Binomial LLR
    rate_in  <- cz_mat / pmax(nz_mat, 1e-300)
    rate_out <- cz_bar / pmax(nz_bar_mat, 1e-300)
    mask <- rate_in > rate_out & cz_mat > 0

    # L(Ha): inside
    t1 <- ifelse(mask & cz_mat > 0,
                 cz_mat * log(pmax(cz_mat / nz_mat, 1e-300)), 0)
    t1b <- ifelse(mask & (nz_mat - cz_mat) > 0,
                  (nz_mat - cz_mat) * log(pmax(1 - cz_mat / nz_mat, 1e-300)), 0)
    # L(Ha): outside
    t2 <- ifelse(mask & cz_bar > 0,
                 cz_bar * log(pmax(cz_bar / nz_bar_mat, 1e-300)), 0)
    t2b <- ifelse(mask & (nz_bar_mat - cz_bar) > 0,
                  (nz_bar_mat - cz_bar) * log(pmax(1 - cz_bar / nz_bar_mat, 1e-300)), 0)
    # L(H0)
    p0 <- Cg / N
    t0a <- ifelse(Cg > 0, Cg * log(pmax(p0, 1e-300)), 0)
    t0b <- ifelse((N - Cg) > 0, (N - Cg) * log(pmax(1 - p0, 1e-300)), 0)
    t0 <- t0a + t0b

    llr_mat <- ifelse(mask, t1 + t1b + t2 + t2b - t0, 0)
    llr_mat[llr_mat < 0] <- 0
  } else {
    # Poisson LLR
    expected_bar <- outer(Cg, zone_pop, function(c, p) c - c * p / N)

    mask <- cz_mat > expected & cz_mat > 0 & expected > 0
    t1 <- ifelse(mask & cz_mat > 0,
                 cz_mat * log(pmax(cz_mat, 1e-300) / pmax(expected, 1e-300)), 0)
    t2 <- ifelse(mask & cz_bar > 0,
                 cz_bar * log(pmax(cz_bar, 1e-300) / pmax(expected_bar, 1e-300)), 0)
    llr_mat <- ifelse(mask, t1 + t2, 0)
  }
  llr_mat[Cg == 0, ] <- 0

  # Extract positive pairs
  pos <- which(llr_mat > 0)
  if (length(pos) == 0) {
    return(data.frame(node_id = character(0), center = integer(0),
                      n_regions = integer(0), cases = numeric(0),
                      expected = numeric(0), population = numeric(0),
                      llr = numeric(0), stringsAsFactors = FALSE))
  }

  g_idx <- ((pos - 1) %% n_nodes) + 1
  z_idx <- ((pos - 1) %/% n_nodes) + 1

  df <- data.frame(
    node_id = all_nodes[g_idx],
    center = zone_centers[z_idx],
    n_regions = zone_sizes[z_idx],
    cases = cz_mat[pos],
    expected = expected[pos],
    population = zone_pop[z_idx],
    llr = llr_mat[pos],
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$llr), ]
  # Keep all evaluated pairs (not capped) so that filter_clusters can
  # recover distinct secondary clusters per Cancado et al. (2025) Sec 5.1.1.
  df
}
