#' Tree-Based Scan Statistic
#'
#' Performs the tree-based scan statistic for detecting clusters in
#' hierarchical data. Uses a Poisson or binomial model with Monte Carlo
#' simulation (implemented in C++ via Rcpp) for significance testing.
#'
#' Inputs are supplied through a single \code{data} data.frame with one row
#' per leaf (rows are matched to the tree by \code{node_id}, so they need
#' not be pre-ordered); counts are summed within leaf.
#'
#' @param data A \code{data.frame} with one row per leaf node. The
#'   \code{cases}, \code{node_id} and (optional) \code{population} arguments
#'   name columns of this data.frame.
#' @param cases Column of \code{data} giving case counts at the leaf level.
#'   Given as an unquoted column name (a string or expression also works).
#' @param node_id Column of \code{data} giving the leaf identifier for each
#'   row. Each value must match a leaf of the tree.
#' @param tree A \code{data.frame} with columns \code{node_id} and
#'   \code{parent_id}. Root node(s) must have \code{parent_id = NA}.
#'   Alternatively, pass the tree as parallel vectors via
#'   \code{tree_node_id} and \code{tree_parent_id}.
#' @param population Optional column of \code{data} giving the population at
#'   the leaf level. For the binomial model, \code{population} is the number
#'   of trials (cases + controls) per leaf and is required. For the
#'   Poisson model, defaults to \code{1} per leaf when omitted.
#' @param tree_node_id,tree_parent_id Optional parallel vectors describing
#'   the tree as an alternative to the \code{tree} data.frame. Both must
#'   have the same length, and the root node(s) must have
#'   \code{tree_parent_id = NA}. Ignored when \code{tree} is supplied.
#' @param nsim Integer. Number of Monte Carlo simulations. Default is
#'   \code{999}.
#' @param alpha Numeric. Significance level. Default is \code{0.05}.
#' @param model Character. Likelihood model: either \code{"poisson"}
#'   (default) or \code{"binomial"}.
#' @param seed Integer or \code{NULL}. Random seed for the Monte
#'   Carlo loop. When non-\code{NULL}, the user's pre-existing RNG
#'   state is saved on entry and restored on exit, so the seed
#'   argument affects only this call and does not leak into
#'   subsequent draws in the user's session.
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
#' # One row per leaf (leaves are 4, 5, 6, 7, 8)
#' leaf_data <- data.frame(
#'   node_id = c(4, 5, 6, 7, 8),
#'   cases   = c(50, 5, 3, 2, 4),
#'   pop     = c(100, 100, 100, 100, 100)
#' )
#'
#' result <- tree_scan(leaf_data, cases = cases, node_id = node_id,
#'                     tree = tree, population = pop, nsim = 99)
#' print(result)
tree_scan <- function(data, cases, node_id, tree = NULL, population = NULL,
                      tree_node_id = NULL, tree_parent_id = NULL,
                      nsim = 999L, alpha = 0.05,
                      model = c("poisson", "binomial"),
                      seed = NULL, n_cores = 1L) {

  ## --- resolve data/column inputs (substitutive interface, treeSS >= 0.2.0) ---
  if (missing(data)) {
    stop("`data` is required: a data.frame with one row per leaf.",
         call. = FALSE)
  }
  data <- as.data.frame(data)
  .env <- parent.frame()
  .cases_in <- .resolve_col(substitute(cases),   data, .env, "cases")
  .node_in  <- .resolve_col(substitute(node_id), data, .env, "node_id")
  .q_pop    <- substitute(population)
  .pop_in   <- if (is.null(.q_pop)) NULL else
    .resolve_col(.q_pop, data, .env, "population")

  model <- match.arg(model)
  model_int <- if (model == "binomial") 1L else 0L

  tree <- .normalize_tree(tree, tree_node_id, tree_parent_id)
  .validate_tree(tree)

  leaves <- .get_leaves(tree)
  n_leaves <- length(leaves)

  ## Match the leaf rows of `data` onto the canonical leaf order of the
  ## tree. Cases are summed within leaf; the denominator is taken from the
  ## first occurrence per leaf.
  node_chr   <- as.character(.node_in)
  leaves_chr <- as.character(leaves)
  unknown <- setdiff(unique(node_chr), leaves_chr)
  if (length(unknown) > 0) {
    stop("'node_id' contains values that are not leaves of the tree: ",
         paste(utils::head(unknown, 5), collapse = ", "),
         if (length(unknown) > 5) ", ..." else "", ".", call. = FALSE)
  }
  if (any(is.na(.cases_in))) stop("'cases' contains NA values.", call. = FALSE)
  if (any(.cases_in < 0))    stop("'cases' must be non-negative.", call. = FALSE)

  fac   <- factor(node_chr, levels = leaves_chr)
  cases <- as.numeric(tapply(as.numeric(.cases_in), fac, sum))
  cases[is.na(cases)] <- 0

  if (is.null(.pop_in)) {
    if (model == "binomial") {
      stop("The binomial model requires a `population` column (number of ",
           "trials per leaf, i.e., cases + controls).", call. = FALSE)
    }
    population <- rep(1, n_leaves)
  } else {
    population <- as.numeric(tapply(as.numeric(.pop_in), fac,
                                    function(z) z[1L]))
    population[is.na(population)] <- 0
  }

  if (!is.null(seed)) {
    .snap__ <- .seed_save_and_set(seed)
    on.exit(.seed_restore(.snap__), add = TRUE)
  }

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
    population_supplied = !is.null(.pop_in),
    tree = tree,
    simulated_llr = sim_llr
  )

  class(result) <- c("tree_scan", "treess")
  result
}
