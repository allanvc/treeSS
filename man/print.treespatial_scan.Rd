% Generated manually for treeSS 0.1.41; keep in sync with roxygen comments
% in R/print.R. devtools::document() will regenerate this file from those
% comments without manual intervention.
\name{print.treespatial_scan}
\alias{print.treespatial_scan}
\title{Print Method for treespatial_scan Objects}
\usage{
\method{print}{treespatial_scan}(x, max_show = 10L, ...)
}
\arguments{
\item{x}{An object of class \code{"treespatial_scan"}.}

\item{max_show}{Integer. Maximum number of leaf IDs and region
IDs to display in full before truncating with "... and N more".
The default of \code{10L} keeps console output compact when the
most likely cluster spans many leaves (for example, when the
root node maximises the likelihood ratio for an aggregated
denominator) or many regions; set \code{max_show = -1L} to
print every value. Mirrors the convention used by \pkg{tibble}.}

\item{...}{Further arguments passed to or from other methods.}
}
\description{
Print Method for treespatial_scan Objects
}
