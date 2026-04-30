#' Tree-Based Scan Statistic
#'
#' Performs the tree-based scan statistic for detecting clusters in
#' hierarchical data. Uses a Poisson or binomial model with Monte Carlo
#' simulation (implemented in C++ via Rcpp) for significance testing.
#'
#' @param tree A \code{data.frame} with columns \code{node_id} and
#'   \code{parent_id}. Root node(s) must have \code{parent_id = NA}.
#'   Alternatively, pass the tree as parallel vectors via
#'   \code{tree_node_id} and \code{tree_parent_id}.
#' @param tree_node_id,tree_parent_id Optional parallel vectors describing
#'   the tree as an alternative to the \code{tree} data.frame. Both must
#'   have the same length, and the root node(s) must have
#'   \code{tree_parent_id = NA}. Ignored when \code{tree} is supplied.
#' @param cases A numeric vector of case counts at the leaf level.
#' @param population A numeric vector of population at the leaf level, or a
#'   single value. For the binomial model, \code{population} is the number
#'   of trials (cases + controls) per leaf and is required. For the
#'   Poisson model, defaults to \code{1} per leaf if \code{NULL}.
#' @param nsim Integer. Number of Monte Carlo simulations. Default is
#'   \code{999}.
#' @param alpha Numeric. Significance level. Default is \code{0.05}.
#' @param model Character. Likelihood model: either \code{"poisson"}
#'   (default) or \code{"binomial"}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#' @param n_cores Integer. Number of OpenMP threads for the Monte Carlo
#'   loop. Default is \code{1L} (serial). Set higher to parallelize.
#'
#' @return An object of class \code{"tree_scan"} (see package help for
#'   details).
#'
#' @references
#' Kulldorff, M., Fang, Z., & Walsh, S. J. (2003). A tree-based scan
#' statistic for database disease surveillance. \emph{Biometrics}, 59(2),
#' 323–331.
#'
#' @seealso \code{\link{circular_scan}}, \code{\link{treespatial_scan}},
#'   \code{\link{aggregate_tree}}
#'
#' @export
#' @examples
#' tree <- data.frame(
#'   node_id   = c(1, 2, 3, 4, 5, 6, 7, 8),
#'   parent_id = c(NA, 1, 1, 2, 2, 3, 3, 3)
#' )
#' cases <- c(50, 5, 3, 2, 4)
#' pop   <- c(100, 100, 100, 100, 100)
#'
#' result <- tree_scan(tree, cases, population = pop, nsim = 99)
#' print(result)
tree_scan <- function(tree = NULL, cases, population = NULL, nsim = 999L,
                      alpha = 0.05,
                      model = c("poisson", "binomial"),
                      seed = NULL, n_cores = 1L,
                      tree_node_id = NULL, tree_parent_id = NULL) {

  model <- match.arg(model)
  model_int <- if (model == "binomial") 1L else 0L

  tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)
  .validate_tree(tree)

  leaves <- .get_leaves(tree)
  n_leaves <- length(leaves)

  if (length(cases) != n_leaves) {
    stop("Length of 'cases' must equal the number of leaf nodes (",
         n_leaves, ").", call. = FALSE)
  }

  if (is.null(population)) {
    if (model == "binomial") {
      stop("The binomial model requires 'population' (number of trials ",
           "per leaf, i.e., cases + controls). Please supply a numeric ",
           "vector.", call. = FALSE)
    }
    population <- rep(1, n_leaves)
  } else if (length(population) == 1) {
    population <- rep(population, n_leaves)
  }
  if (length(population) != n_leaves) {
    stop("Length of 'population' must equal the number of leaf nodes.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  all_nodes <- tree$node_id
  n_nodes <- length(all_nodes)

  # Aggregate cases and population upward
  leaf_idx <- match(leaves, all_nodes)
  node_cases <- rep(0, n_nodes)
  node_pop <- rep(0, n_nodes)
  node_cases[leaf_idx] <- cases
  node_pop[leaf_idx] <- population

  depths <- .compute_depths(tree)
  proc_order <- order(depths, decreasing = TRUE)

  for (idx in proc_order) {
    node <- all_nodes[idx]
    children_idx_v <- which(tree$parent_id == node & !is.na(tree$parent_id))
    if (length(children_idx_v) > 0) {
      child_match <- match(tree$node_id[children_idx_v], all_nodes)
      node_cases[idx] <- sum(node_cases[child_match])
      node_pop[idx] <- sum(node_pop[child_match])
    }
  }

  C <- sum(cases)
  N <- sum(population)

  # --- C++ backend for MC simulation ---
  tree_csr <- .tree_to_csr_children(tree)
  leaf_idx_0 <- as.integer(leaf_idx - 1L)
  proc_order_0 <- as.integer(proc_order - 1L)

  cpp_result <- mc_treescan_cpp(
    as.numeric(node_cases),
    as.numeric(node_pop),
    as.numeric(C),
    as.numeric(N),
    leaf_idx_0,
    as.numeric(population),
    tree_csr$children_idx,
    tree_csr$children_ptr,
    proc_order_0,
    as.integer(nsim),
    model_int,
    as.integer(max(1L, n_cores))
  )

  obs_llr_vec <- cpp_result$obs_llr
  sim_llr <- cpp_result$sim_llr
  best_idx <- cpp_result$obs_best  # 1-based

  # P-values for all cuts
  pvalues <- numeric(n_nodes)
  for (i in seq_len(n_nodes)) {
    pvalues[i] <- (sum(sim_llr >= obs_llr_vec[i]) + 1) / (nsim + 1)
  }

  # Build cuts table
  expected_vec <- C * node_pop / N

  all_cuts <- data.frame(
    node_id = all_nodes,
    cases = node_cases,
    expected = expected_vec,
    population = node_pop,
    llr = obs_llr_vec,
    pvalue = pvalues,
    stringsAsFactors = FALSE
  )
  all_cuts <- all_cuts[order(-all_cuts$llr), ]
  sig_cuts <- all_cuts[all_cuts$pvalue < alpha, ]

  branches <- .get_branches(tree)
  best_node <- all_nodes[best_idx]
  best_leaves <- branches[[as.character(best_node)]]$leaves

  result <- list(
    most_likely_cluster = list(
      node_id = best_node,
      leaf_ids = best_leaves,
      cases = node_cases[best_idx],
      expected = expected_vec[best_idx],
      population = node_pop[best_idx],
      llr = obs_llr_vec[best_idx]
    ),
    all_cuts = all_cuts,
    significant_cuts = sig_cuts,
    pvalue = pvalues[best_idx],
    alpha = alpha,
    nsim = nsim,
    total_cases = C,
    total_population = N,
    tree = tree,
    simulated_llr = sim_llr
  )

  class(result) <- "tree_scan"
  result
}
