% Generated manually for treeSS 0.1.41; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{summary.treespatial_scan}
\alias{summary.treespatial_scan}
\title{Summary Method for treespatial_scan Objects}
\usage{
\method{summary}{treespatial_scan}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"treespatial_scan"}.}

\item{...}{Further arguments passed to
\code{print.treespatial_scan()}. In particular, \code{max_show}
controls how many leaf IDs and region IDs of the most-likely
cluster are printed in full before truncating with \code{"...
and N more"}; see \code{\link{print.treespatial_scan}} for
details.}
}
\description{
Prints the same fields as \code{print.treespatial_scan()}, followed
by the top-10 branches by total case count and a numeric summary
of the simulated log-likelihood-ratio distribution under the null.
Useful for diagnosing how heavily one branch of the tree dominates
the dataset and whether the observed maximum LR is a clear outlier
against the simulated null.
}
