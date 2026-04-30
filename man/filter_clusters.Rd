\name{filter_clusters}
\alias{filter_clusters}
\title{Filter Overlapping Secondary Clusters}
\usage{
filter_clusters(result, alpha = NULL)
}
\arguments{
\item{result}{An object of class \code{"treespatial_scan"},
  \code{"circular_scan"}, or \code{"tree_scan"}.}
\item{alpha}{Numeric. Significance level. Default is \code{NULL} (uses
  alpha from the original analysis).}
}
\value{
A \code{data.frame} of distinct significant clusters ordered by
  descending LLR.
}
\description{
Removes overlapping secondary clusters from the output of
\code{\link{treespatial_scan}}, \code{\link{circular_scan}}, or
\code{\link{tree_scan}}.
}
\details{
\strong{treespatial_scan}: two clusters overlap if they have BOTH tree
overlap (ancestor/descendant) AND spatial overlap (same center or
identical regions). Follows Cancado et al. (2025).

\strong{circular_scan}: Kulldorff (1997) criterion --- center of candidate
not in any retained zone, and no retained center in the candidate zone.

\strong{tree_scan}: a secondary cluster is distinct if its node is NOT
an ancestor or descendant of any previously retained node. Follows
Kulldorff et al. (2003).
}
\references{
Kulldorff, M. (1997). A spatial scan statistic.
\emph{Communications in Statistics}, 26(6), 1481-1496.

Kulldorff, M., Fang, Z., & Walsh, S. J. (2003). A tree-based scan
statistic for database disease surveillance. \emph{Biometrics}, 59(2),
323-331.

Cancado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
(2025). A tree-spatial scan statistic. \emph{Environmental and Ecological
Statistics}, 32, 953--978. \doi{10.1007/s10651-025-00670-w}
}
\seealso{
\code{\link{treespatial_scan}}, \code{\link{circular_scan}},
\code{\link{tree_scan}}, \code{\link{iterative_scan}}
}
\examples{
data(london_collisions); data(london_tree)
result <- treespatial_scan(
  cases       = london_collisions$cases,
  population  = london_collisions$population,
  region_id   = london_collisions$region_id,
  x           = london_collisions$x,
  y           = london_collisions$y,
  node_id     = london_collisions$node_id,
  tree        = london_tree,
  nsim        = 99, seed = 42
)
filter_clusters(result)
}
