#' treeSS: Tree-Spatial Scan Statistic for Cluster Detection
#'
#' Implements the tree-spatial scan statistic for detecting clusters that
#' combine both spatial and hierarchical structures, as proposed by Cançado
#' et al. (2025). The method extends Kulldorff's (1997) circular spatial scan
#' statistic and the tree-based scan statistic (Kulldorff et al. 2003) by
#' searching for anomalies in both geographic regions and branches of
#' hierarchical trees simultaneously.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{treespatial_scan}}}{Performs the tree-spatial scan
#'     statistic, detecting clusters defined by pairs of spatial zones and
#'     tree branches.}
#'   \item{\code{\link{circular_scan}}}{Performs Kulldorff's circular spatial
#'     scan statistic.}
#'   \item{\code{\link{tree_scan}}}{Performs the tree-based scan statistic.}
#' }
#'
#' @section Helper functions:
#' \describe{
#'   \item{\code{\link{build_zones}}}{Constructs candidate spatial zones using
#'     circular windows centered on region centroids.}
#'   \item{\code{\link{aggregate_tree}}}{Aggregates case counts from leaves to
#'     internal nodes of the hierarchical tree.}
#'   \item{\code{\link{filter_clusters}}}{Removes overlapping secondary
#'     clusters from the tree-spatial scan output.}
#' }
#'
#' @references
#' Cançado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
#' (2025). A tree-spatial scan statistic. \emph{Environmental and Ecological
#' Statistics}, 32, 953–978. \doi{10.1007/s10651-025-00670-w}
#'
#' Kulldorff, M. (1997). A spatial scan statistic. \emph{Communications in
#' Statistics - Theory and Methods}, 26(6), 1481–1496.
#'
#' Kulldorff, M., Fang, Z., & Walsh, S. J. (2003). A tree-based scan
#' statistic for database disease surveillance. \emph{Biometrics}, 59(2),
#' 323–331.
#'
#' @docType package
#' @name treeSS-package
#' @useDynLib treeSS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils head
"_PACKAGE"
