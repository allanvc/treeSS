#' Print Method for circular_scan Objects
#'
#' @param x An object of class \code{"circular_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.circular_scan <- function(x, ...) {
  cat("Circular Spatial Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total cases:", x$total_cases, "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of regions:", nrow(x$regions), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  cat("  Regions:", paste(mlc$region_ids, collapse = ", "), "\n")
  cat("  Number of regions:", length(mlc$region_ids), "\n")
  cat("  Cases:", mlc$cases, "\n")
  cat("  Expected:", round(mlc$expected, 2), "\n")
  cat("  Population:", mlc$population, "\n")
  cat("  Relative risk:", round(mlc$rr, 4), "\n")
  cat("  Log-LR:", round(mlc$llr, 4), "\n")
  cat("  P-value:", format.pval(x$pvalue, digits = 4), "\n")

  if (!is.null(x$secondary_clusters) && nrow(x$secondary_clusters) > 0) {
    cat("  Secondary candidates:", nrow(x$secondary_clusters),
        "(use filter_clusters() to extract distinct clusters)\n")
  }

  invisible(x)
}

#' Summary Method for circular_scan Objects
#'
#' @param object An object of class \code{"circular_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.circular_scan <- function(object, ...) {
  print.circular_scan(object, ...)
  cat("\nSimulated LLR distribution:\n")
  print(summary(object$simulated_llr))
  invisible(object)
}

#' Print Method for tree_scan Objects
#'
#' @param x An object of class \code{"tree_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.tree_scan <- function(x, ...) {
  cat("Tree-Based Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total cases:", x$total_cases, "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of nodes:", nrow(x$tree), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  cat("  Node ID:", mlc$node_id, "\n")
  cat("  Leaf IDs:", paste(mlc$leaf_ids, collapse = ", "), "\n")
  cat("  Cases:", mlc$cases, "\n")
  cat("  Expected:", round(mlc$expected, 2), "\n")
  cat("  Log-LR:", round(mlc$llr, 4), "\n")
  cat("  P-value:", format.pval(x$pvalue, digits = 4), "\n")

  nsig <- nrow(x$significant_cuts)
  cat("\nSignificant cuts (alpha =", x$alpha, "):", nsig, "\n")

  if (nsig > 0) {
    top <- head(x$significant_cuts, 10)
    cat("\nTop significant cuts:\n")
    print(top, row.names = FALSE)
  }

  invisible(x)
}

#' Summary Method for tree_scan Objects
#'
#' @param object An object of class \code{"tree_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.tree_scan <- function(object, ...) {
  print.tree_scan(object, ...)
  cat("\nAll cuts:\n")
  print(object$all_cuts, row.names = FALSE)
  invisible(object)
}

#' Print Method for treespatial_scan Objects
#'
#' @param x An object of class \code{"treespatial_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.treespatial_scan <- function(x, ...) {
  cat("Tree-Spatial Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of regions:", nrow(x$regions), "\n")
  cat("Number of tree nodes:", nrow(x$tree), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  cat("  Tree node:", mlc$node_id, "\n")
  cat("  Leaf IDs:", paste(mlc$leaf_ids, collapse = ", "), "\n")
  cat("  Regions:", paste(mlc$region_ids, collapse = ", "), "\n")
  cat("  Number of regions:", length(mlc$region_ids), "\n")
  cat("  Cases in (zone, branch):", mlc$cases, "\n")
  cat("  Expected:", round(mlc$expected, 2), "\n")
  cat("  Population:", mlc$population, "\n")
  cat("  Relative risk:", round(mlc$rr, 4), "\n")
  cat("  Log-LR:", round(mlc$llr, 4), "\n")
  cat("  P-value:", format.pval(x$pvalue, digits = 4), "\n")

  invisible(x)
}

#' Summary Method for treespatial_scan Objects
#'
#' @param object An object of class \code{"treespatial_scan"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.treespatial_scan <- function(object, ...) {
  print.treespatial_scan(object, ...)

  cat("\nTotal cases by branch (top 10):\n")
  cg <- sort(object$total_cases_by_branch, decreasing = TRUE)
  print(head(cg, 10))

  cat("\nSimulated LLR distribution:\n")
  print(summary(object$simulated_llr))

  invisible(object)
}
