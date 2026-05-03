% Generated manually for treeSS 0.1.41; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{summary.circular_scan}
\alias{summary.circular_scan}
\title{Summary Method for circular_scan Objects}
\usage{
\method{summary}{circular_scan}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"circular_scan"}.}

\item{...}{Further arguments passed to \code{print.circular_scan()}.
In particular, \code{max_show} controls how many region IDs of
the most-likely cluster are printed in full before truncating
with \code{"... and N more"}; see
\code{\link{print.circular_scan}} for details.}
}
\description{
Prints the same fields as \code{print.circular_scan()}, followed
by a numeric summary of the simulated log-likelihood-ratio
distribution under the null. Useful for checking that the Monte
Carlo replicates produced a reasonable spread relative to the
observed test statistic.
}
