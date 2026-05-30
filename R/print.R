## .cat_wrapped(): print "  Label: v1, v2, v3, ..." with automatic
## line-wrapping at 65 characters (so the output fits in the R Journal
## two-column PDF format) and optional truncation when there are too
## many values to print in full.
##
## When `length(values) > max_show`, only the first `max_show` are
## shown and a tail of "... and N more" is appended. Pass `max_show
## = -1L` (or any value larger than `length(values)`) to print the
## full list. This mirrors the way `tibble` prints data frames:
## sensible default, with an explicit opt-in for the full output.
##
## Typical use:
##   .cat_wrapped("  Leaf IDs:", mlc$leaf_ids, max_show = 10L)
.cat_wrapped <- function(label, values, sep = ", ", max_show = 10L) {
  n <- length(values)
  if (n == 0L) {
    cat(label, "\n", sep = "")
    return(invisible())
  }
  truncated <- max_show > 0L && n > max_show
  if (truncated) {
    shown <- values[seq_len(max_show)]
    text  <- paste0(paste(shown, collapse = sep),
                    sep, "... and ", n - max_show, " more")
  } else {
    text  <- paste(values, collapse = sep)
  }
  width <- min(getOption("width", 80L), 65L)
  # Indent continuation lines with the same number of spaces as `label`
  # plus one for the space that follows it.
  indent <- paste(rep(" ", nchar(label) + 1L), collapse = "")
  lines <- strwrap(text, width = width,
                   initial = paste0(label, " "),
                   prefix  = indent)
  cat(lines, sep = "\n")
  cat("\n")
  invisible()
}

#' Print Method for circular_scan Objects
#'
#' @param x An object of class \code{"circular_scan"}.
#' @param max_show Integer. Maximum number of region IDs to display
#'   in full before truncating with "... and N more". The default of
#'   \code{10L} keeps console output compact for large clusters; set
#'   \code{max_show = -1L} (or any value larger than the cluster
#'   size) to print every region. Mirrors the convention used by
#'   \pkg{tibble} for printing wide tables.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x} (the input object of class
#'   \code{"circular_scan"}), unchanged. Called for the side effect
#'   of printing a human-readable summary of the scan result to the
#'   console: the total cases and population, the number of regions
#'   scanned, the number of Monte Carlo simulations, and the
#'   contents of the most likely cluster (region IDs, cases,
#'   expected count, population, relative risk, log-likelihood
#'   ratio, and p-value).
#' @export
print.circular_scan <- function(x, max_show = 10L, ...) {
  cat("Circular Spatial Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total cases:", x$total_cases, "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of regions:", nrow(x$regions), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  .cat_wrapped("  Regions:", mlc$region_ids, max_show = max_show)
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
#' Prints the same fields as \code{print.circular_scan()}, followed
#' by a numeric summary of the simulated log-likelihood-ratio
#' distribution under the null. Useful for checking that the Monte
#' Carlo replicates produced a reasonable spread relative to the
#' observed test statistic.
#'
#' @param object An object of class \code{"circular_scan"}.
#' @param ... Further arguments passed to \code{print.circular_scan()}.
#'   In particular, \code{max_show} controls how many region IDs of
#'   the most-likely cluster are printed in full before truncating
#'   with \code{"... and N more"}; see
#'   \code{\link{print.circular_scan}} for details.
#'
#' @return Invisibly returns \code{object} (the input object of
#'   class \code{"circular_scan"}), unchanged. Called for the side
#'   effect of printing the same fields as
#'   \code{print.circular_scan()}, followed by a numeric summary
#'   (min, quartiles, median, mean, max) of the simulated
#'   log-likelihood-ratio distribution under the null.
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
#' @param max_show Integer. Maximum number of leaf IDs to display
#'   in full before truncating with "... and N more". The default of
#'   \code{10L} keeps console output compact when a node spans
#'   hundreds or thousands of leaves (such as the root of a deep
#'   tree); set \code{max_show = -1L} to print every leaf.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x} (the input object of class
#'   \code{"tree_scan"}), unchanged. Called for the side effect of
#'   printing a human-readable summary of the scan result to the
#'   console: the total cases and population, the number of tree
#'   nodes, the number of Monte Carlo simulations, the contents of
#'   the most likely cluster (node ID, leaf IDs, cases, expected
#'   count, log-likelihood ratio, and p-value), and the top
#'   significant cuts at the chosen significance threshold.
#' @export
print.tree_scan <- function(x, max_show = 10L, ...) {
  cat("Tree-Based Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total cases:", x$total_cases, "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of nodes:", nrow(x$tree), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  cat("  Node ID:", mlc$node_id, "\n")
  .cat_wrapped("  Leaf IDs:", mlc$leaf_ids, max_show = max_show)
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
#' Prints the same fields as \code{print.tree_scan()}, followed by
#' the full table of all candidate cuts (not just significant ones)
#' ordered by decreasing log-likelihood ratio.
#'
#' @param object An object of class \code{"tree_scan"}.
#' @param ... Further arguments passed to \code{print.tree_scan()}.
#'   In particular, \code{max_show} controls how many leaf IDs of
#'   the most-likely cluster are printed in full before truncating
#'   with \code{"... and N more"}; see \code{\link{print.tree_scan}}
#'   for details.
#'
#' @return Invisibly returns \code{object} (the input object of
#'   class \code{"tree_scan"}), unchanged. Called for the side
#'   effect of printing the same fields as
#'   \code{print.tree_scan()}, followed by the full table of all
#'   candidate cuts (not just significant ones) ordered by
#'   decreasing log-likelihood ratio.
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
#' @param max_show Integer. Maximum number of leaf IDs and region
#'   IDs to display in full before truncating with "... and N more".
#'   The default of \code{10L} keeps console output compact when the
#'   most likely cluster spans many leaves (for example, when the
#'   root node maximises the likelihood ratio for an aggregated
#'   denominator) or many regions; set \code{max_show = -1L} to
#'   print every value. Mirrors the convention used by \pkg{tibble}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x} (the input object of class
#'   \code{"treespatial_scan"}), unchanged. Called for the side
#'   effect of printing a human-readable summary of the scan result
#'   to the console: the total population, the number of regions
#'   and tree nodes, the number of Monte Carlo simulations, and the
#'   contents of the most likely cluster (tree node, leaf IDs,
#'   region IDs, cases, expected count, population, relative risk,
#'   log-likelihood ratio, and p-value).
#' @export
print.treespatial_scan <- function(x, max_show = 10L, ...) {
  cat("Tree-Spatial Scan Statistic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Total population:", x$total_population, "\n")
  cat("Number of regions:", nrow(x$regions), "\n")
  cat("Number of tree nodes:", nrow(x$tree), "\n")
  cat("Monte Carlo simulations:", x$nsim, "\n\n")

  mlc <- x$most_likely_cluster
  cat("Most likely cluster:\n")
  cat("  Tree node:", mlc$node_id, "\n")
  .cat_wrapped("  Leaf IDs:", mlc$leaf_ids, max_show = max_show)
  .cat_wrapped("  Regions:", mlc$region_ids, max_show = max_show)
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
#' Prints the same fields as \code{print.treespatial_scan()}, followed
#' by the top-10 branches by total case count and a numeric summary
#' of the simulated log-likelihood-ratio distribution under the null.
#' Useful for diagnosing how heavily one branch of the tree dominates
#' the dataset and whether the observed maximum LR is a clear outlier
#' against the simulated null.
#'
#' @param object An object of class \code{"treespatial_scan"}.
#' @param ... Further arguments passed to
#'   \code{print.treespatial_scan()}. In particular, \code{max_show}
#'   controls how many leaf IDs and region IDs of the most-likely
#'   cluster are printed in full before truncating with \code{"...
#'   and N more"}; see \code{\link{print.treespatial_scan}} for
#'   details.
#'
#' @return Invisibly returns \code{object} (the input object of
#'   class \code{"treespatial_scan"}), unchanged. Called for the
#'   side effect of printing the same fields as
#'   \code{print.treespatial_scan()}, followed by the top-10
#'   branches by total case count and a numeric summary (min,
#'   quartiles, median, mean, max) of the simulated
#'   log-likelihood-ratio distribution under the null.
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

#' Print Method for sequential_scan Objects
#'
#' @param x An object of class \code{"sequential_scan"}.
#' @param max_show Integer. Maximum number of region IDs (or leaf IDs,
#'   for tree-only scans) to display in full per iteration before
#'   truncating with \code{"... and N more"}. Default \code{10L}; set
#'   \code{max_show = -1L} to print every value.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x} (the input object of class
#'   \code{"sequential_scan"}), unchanged. Called for the side effect
#'   of printing a human-readable summary of the sequential scan
#'   result to the console: the scan type, the buffer size, the
#'   number of iterations performed, the significance threshold, the
#'   number of significant clusters, and a per-iteration table of
#'   detected clusters.
#' @export
print.sequential_scan <- function(x, max_show = 10L, ...) {
  cat("Sequential Scan (Zhang, Assuncao & Kulldorff, 2010)\n")
  cat(paste(rep("-", 55), collapse = ""), "\n")
  cat("Scan type           :", x$scan_type, "\n")
  cat("Buffer size         :", x$buffer_size,
      "neighbouring region(s) per iteration\n")
  cat("Iterations performed:", x$n_iter, "\n")
  cat("Alpha               :", x$alpha, "\n")
  cat("MC replications/iter:", x$nsim, "\n\n")

  if (nrow(x$clusters) > 0L) {
    n_sig <- sum(x$clusters$significant, na.rm = TRUE)
    cat("Significant clusters: ", n_sig, " of ",
        nrow(x$clusters), "\n\n", sep = "")
    display <- x$clusters
    # Hide list-columns from the tabular print
    if (!isTRUE(max_show < 0L)) {
      display$region_ids <- NULL
      display$leaf_ids   <- NULL
    }
    print(display, row.names = FALSE)
  } else {
    cat("No clusters detected.\n")
  }
  invisible(x)
}


#' Summary Method for sequential_scan Objects
#'
#' @param object An object of class \code{"sequential_scan"}.
#' @param ... Passed to \code{print.sequential_scan()}.
#'
#' @return Invisibly returns \code{object}. Called for the side effect
#'   of printing the same fields as \code{print.sequential_scan()}
#'   followed by, for each iteration, the cluster composition
#'   (region IDs or leaf IDs).
#' @export
summary.sequential_scan <- function(object, ...) {
  print.sequential_scan(object, ...)
  if (nrow(object$clusters) > 0L) {
    cat("\nCluster composition per iteration:\n")
    cl <- object$clusters
    for (k in seq_len(nrow(cl))) {
      cat(sprintf("\n  Iteration %d:\n", cl$iteration[k]))
      if ("region_ids" %in% names(cl)) {
        .cat_wrapped("    Regions :", cl$region_ids[[k]])
      }
      if ("leaf_ids" %in% names(cl)) {
        .cat_wrapped("    Leaves  :", cl$leaf_ids[[k]])
      }
    }
  }
  invisible(object)
}

