\name{chicago_map}
\alias{chicago_map}
\docType{data}
\title{Chicago Community Area Boundaries}
\description{
Simplified polygon boundaries for the 77 community areas in Chicago.
}
\format{
An \code{sf} data.frame with 77 rows and 3 columns:
\describe{
  \item{NAME}{Community area name.}
  \item{AREA_NUM}{Community area number (integer, 1--77).}
  \item{geometry}{MULTIPOLYGON geometry in WGS84 (EPSG:4326).}
}
}
\source{
City of Chicago Data Portal, Community Area Boundaries.
}
\seealso{\code{\link{chicago_crimes}}}
\examples{
data(chicago_map)
cat("Community areas:", nrow(chicago_map), "\n")
}
\keyword{datasets}
