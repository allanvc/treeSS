% Generated manually for treeSS 0.1.44; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{print.tree_scan}
\alias{print.tree_scan}
\title{Print Method for tree_scan Objects}
\usage{
\method{print}{tree_scan}(x, max_show = 10L, ...)
}
\arguments{
\item{x}{An object of class \code{"tree_scan"}.}

\item{max_show}{Integer. Maximum number of leaf IDs to display
in full before truncating with "... and N more". The default of
\code{10L} keeps console output compact when a node spans
hundreds or thousands of leaves (such as the root of a deep
tree); set \code{max_show = -1L} to print every leaf.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Invisibly returns \code{x} (the input object of class
\code{"tree_scan"}), unchanged. Called for the side effect of
printing a human-readable summary of the scan result to the
console: the total cases and population, the number of tree
nodes, the number of Monte Carlo simulations, the contents of
the most likely cluster (node ID, leaf IDs, cases, expected
count, log-likelihood ratio, and p-value), and the top
significant cuts at the chosen significance threshold.
}
\description{
Print Method for tree_scan Objects
}
