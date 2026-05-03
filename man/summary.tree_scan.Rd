% Generated manually for treeSS 0.1.41; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{summary.tree_scan}
\alias{summary.tree_scan}
\title{Summary Method for tree_scan Objects}
\usage{
\method{summary}{tree_scan}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"tree_scan"}.}

\item{...}{Further arguments passed to \code{print.tree_scan()}.
In particular, \code{max_show} controls how many leaf IDs of
the most-likely cluster are printed in full before truncating
with \code{"... and N more"}; see \code{\link{print.tree_scan}}
for details.}
}
\description{
Prints the same fields as \code{print.tree_scan()}, followed by
the full table of all candidate cuts (not just significant ones)
ordered by decreasing log-likelihood ratio.
}
