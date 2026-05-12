% Generated manually for treeSS 0.1.44; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{print.circular_scan}
\alias{print.circular_scan}
\title{Print Method for circular_scan Objects}
\usage{
\method{print}{circular_scan}(x, max_show = 10L, ...)
}
\arguments{
\item{x}{An object of class \code{"circular_scan"}.}

\item{max_show}{Integer. Maximum number of region IDs to display
in full before truncating with "... and N more". The default of
\code{10L} keeps console output compact for large clusters; set
\code{max_show = -1L} (or any value larger than the cluster
size) to print every region. Mirrors the convention used by
\pkg{tibble} for printing wide tables.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Invisibly returns \code{x} (the input object of class
\code{"circular_scan"}), unchanged. Called for the side effect
of printing a human-readable summary of the scan result to the
console: the total cases and population, the number of regions
scanned, the number of Monte Carlo simulations, and the
contents of the most likely cluster (region IDs, cases,
expected count, population, relative risk, log-likelihood
ratio, and p-value).
}
\description{
Print Method for circular_scan Objects
}
