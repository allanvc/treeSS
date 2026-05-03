% Generated manually for treeSS 0.1.43; keep in sync with roxygen comments
% in R/iterative_scan.R. devtools::document() will regenerate this file from
% those comments without manual intervention.
\name{print.iterative_scan}
\alias{print.iterative_scan}
\title{Print Method for iterative_scan Objects}
\usage{
\method{print}{iterative_scan}(x, max_show = 10L, ...)
}
\arguments{
\item{x}{An object of class \code{"iterative_scan"}.}

\item{max_show}{Integer. Accepted for API consistency with
\code{print.treespatial_scan()}, \code{print.tree_scan()}, and
\code{print.circular_scan()}. Currently the only field affected
is the per-iteration \code{leaf_ids} column of \code{x$clusters},
which is omitted from the printed table by default to keep the
summary compact; pass \code{max_show = -1L} to include both
\code{leaf_ids} and \code{region_ids} in the printed table.}

\item{...}{Further arguments passed to or from other methods.}
}
\description{
Print Method for iterative_scan Objects
}
