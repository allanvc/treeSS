\name{london_boroughs_map}
\alias{london_boroughs_map}
\docType{data}
\title{London Borough Boundaries}
\description{
Simplified polygon boundaries for the 33 London boroughs. Use with
\code{\link{london_collisions}} for map-based visualization of cluster results.
}
\format{
An \code{sf} data.frame with 33 rows and 3 columns:
\describe{
  \item{NAME}{Borough name, e.g. \code{"Westminster"}.}
  \item{GSS_CODE}{ONS geography code, e.g. \code{"E09000033"}.}
  \item{geometry}{MULTIPOLYGON geometry in WGS84 (EPSG:4326).}
}
}
\details{
Geometries are simplified (\code{st_simplify(dTolerance = 0.001)}) to
reduce file size. For full-resolution boundaries, download from the
source URL.

To merge with cluster results from \code{\link{get_cluster_regions}},
use: \code{merge(london_boroughs_map, cr, by.x = "NAME", by.y = "borough")}.
}
\source{
London Datastore, Statistical GIS Boundary Files for London
(\url{https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london}).
Contains National Statistics data (c) Crown copyright and database right 2014.
}
\seealso{\code{\link{london_collisions}}}
\examples{
data(london_boroughs_map)
cat("Boroughs:", nrow(london_boroughs_map), "\n")
}
\keyword{datasets}
