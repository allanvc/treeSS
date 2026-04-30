#' Build Candidate Spatial Zones Using Circular Windows
#'
#' Constructs the set of candidate spatial zones for the circular scan
#' statistic. Each zone is defined by a circular window centered on a region's
#' centroid, containing all regions whose centroids fall within the circle.
#'
#' @param regions A \code{data.frame} with columns \code{region_id},
#'   \code{population}, \code{x}, and \code{y}. Coordinates \code{x} and
#'   \code{y} represent centroids of each region.
#' @param max_pop A numeric value specifying the maximum population allowed
#'   inside a zone. If \code{NULL} (default), it is set to 50\% of the total
#'   population.
#'
#' @return A \code{list} of zones. Each zone is a list with elements:
#'   \describe{
#'     \item{center}{Integer index of the center region.}
#'     \item{region_idx}{Integer vector of region indices inside the zone.}
#'     \item{population}{Total population inside the zone.}
#'   }
#'
#' @details
#' For each region \eqn{j}, regions are sorted by distance from \eqn{j}. Regions
#' are added incrementally to the zone as long as the total population does
#' not exceed \code{max_pop}. This generates a nested set of circular windows
#' centered on each region.
#'
#' @references
#' Kulldorff, M. (1997). A spatial scan statistic. \emph{Communications in
#' Statistics - Theory and Methods}, 26(6), 1481–1496.
#'
#' @export
#' @examples
#' regions <- data.frame(
#'   region_id = 1:5,
#'   population = c(1000, 2000, 1500, 800, 1200),
#'   x = c(0, 1, 2, 0.5, 1.5),
#'   y = c(0, 0, 0, 1, 1)
#' )
#' zones <- build_zones(regions, max_pop = 3000)
#' length(zones) # number of candidate zones
build_zones <- function(regions, max_pop = NULL) {

  .validate_regions(regions)

  n <- nrow(regions)
  N <- sum(regions$population)

  if (is.null(max_pop)) {
    max_pop <- 0.5 * N
  }

  # Compute distance matrix
  coords <- as.matrix(regions[, c("x", "y")])
  dist_mat <- as.matrix(stats::dist(coords, method = "euclidean"))

  zones <- list()
  zone_id <- 1L

  for (center in seq_len(n)) {
    ordered_idx <- order(dist_mat[center, ])
    zone_regions <- integer(0)
    zone_pop <- 0

    for (i in ordered_idx) {
      new_pop <- zone_pop + regions$population[i]

      if (new_pop <= max_pop) {
        zone_regions <- c(zone_regions, i)
        zone_pop <- new_pop

        zones[[zone_id]] <- list(
          center = center,
          region_idx = zone_regions,
          population = zone_pop
        )
        zone_id <- zone_id + 1L
      }
    }
  }

  zones
}
