## =============================================================================
## Documentation for datasets shipped with treeSS (>= 0.1.4).
##
## Each topical example ships as a pair:
##     <prefix>_<topic>   - combined long-format data.frame (input to scans)
##     <prefix>_tree      - hierarchical tree (node_id, parent_id)
##
## Plus the standalone datasets:
##     fl_deaths             - raw long-format death counts (Florida; user
##                             builds the tree)
##     chicago_map           - sf polygons for Chicago community areas
##     london_boroughs_map   - sf polygons for London boroughs
## =============================================================================


# ==== Rio de Janeiro =========================================================

#' Rio de Janeiro Infant Mortality (2016)
#'
#' Combined long-format dataset of infant mortality in the 92 municipalities
#' of Rio de Janeiro state, Brazil, 2016. Each row represents a (municipality,
#' ICD-10 leaf) combination with at least one death; missing combinations
#' are implicitly zero. Use together with \code{\link{rj_tree}}.
#'
#' @format A data.frame in long format with columns:
#' \describe{
#'   \item{region_id}{Integer identifier of the municipality.}
#'   \item{ibge_code}{6-digit IBGE municipality code.}
#'   \item{name}{Municipality name.}
#'   \item{population}{IBGE municipal population estimate for 2016.}
#'   \item{live_births}{Number of live births in 2016 (SINASC). Used as the
#'     population denominator in Section 5.2 of Cancado et al. (2025).}
#'   \item{x}{Longitude of the municipality centroid.}
#'   \item{y}{Latitude of the municipality centroid.}
#'   \item{node_id}{Character ICD-10 leaf code (4 characters), matching a
#'     leaf of \code{\link{rj_tree}}.}
#'   \item{cases}{Integer count of infant deaths.}
#' }
#'
#' @details
#' To reproduce Section 5.2 of Cancado et al. (2025), use \code{live_births}
#' as the population denominator:
#' \preformatted{
#'   data <- rj_mortality
#'   data$population <- data$live_births
#'   treespatial_scan(data, rj_tree, ...)
#' }
#'
#' Municipality polygons can be obtained via \code{geobr::read_municipality()}.
#'
#' @source
#' Population: IBGE municipal estimates. Live births: DATASUS/SINASC via
#' TabNet. Deaths: DATASUS/SIM microdata via OpenDATASUS
#' (\url{https://opendatasus.saude.gov.br}). Centroids: \code{geobr}.
#'
#' @seealso \code{\link{rj_tree}}
#'
#' @references
#' Cancado, A.L.F., Oliveira, G.S., Quadros, A.V.C., Duczmal, L. (2025).
#' A tree-spatial scan statistic. \emph{Environmental and Ecological
#' Statistics}, 32, 953--978. \doi{10.1007/s10651-025-00670-w}
#'
#' @examples
#' data(rj_mortality)
#' head(rj_mortality)
#' cat("Total deaths:", sum(rj_mortality$cases), "\n")
"rj_mortality"


#' Rio de Janeiro Infant Mortality - ICD-10 Tree
#'
#' ICD-10 hierarchy used by the \code{\link{rj_mortality}} dataset.
#' Three levels: Chapter -> 3-character category -> 4-character subcategory.
#'
#' @format A data.frame with 622 rows and 2 columns:
#' \describe{
#'   \item{node_id}{Character identifier of each ICD-10 node.}
#'   \item{parent_id}{Character identifier of the parent node. \code{NA}
#'     for the root.}
#' }
#'
#' @details
#' The tree has 410 leaves (4-character ICD-10 codes) and 622 total nodes.
#'
#' @seealso \code{\link{rj_mortality}}
#'
#' @examples
#' data(rj_tree)
#' head(rj_tree)
"rj_tree"


# ==== London =================================================================

#' London Road Collisions (2022)
#'
#' Combined long-format dataset of road traffic collisions in the 33 London
#' boroughs, 2022. Each row represents a (borough, collision-category)
#' combination with at least one collision. Use together with
#' \code{\link{london_tree}}.
#'
#' @format A data.frame in long format with columns:
#' \describe{
#'   \item{region_id}{Integer identifier of the borough.}
#'   \item{borough}{Borough name.}
#'   \item{population}{Total collisions in the borough (used as population
#'     denominator).}
#'   \item{x}{Longitude of the borough centroid.}
#'   \item{y}{Latitude of the borough centroid.}
#'   \item{node_id}{Character leaf identifier from \code{\link{london_tree}}.}
#'   \item{cases}{Integer count of collisions.}
#' }
#'
#' @seealso \code{\link{london_tree}}, \code{\link{london_boroughs_map}}
#'
#' @source
#' UK Department for Transport, STATS19 data via the \code{stats19} R package.
#'
#' @examples
#' data(london_collisions)
#' head(london_collisions)
#' cat("Total collisions:", sum(london_collisions$cases), "\n")
"london_collisions"


#' London Road Collisions - Hierarchy
#'
#' Hierarchical classification of road collisions by light conditions,
#' road type, and junction detail.
#'
#' @format A data.frame with 81 rows and 2 columns:
#' \describe{
#'   \item{node_id}{Character identifier of each node.}
#'   \item{parent_id}{Character identifier of the parent node. \code{NA}
#'     for the root.}
#' }
#'
#' @details
#' 3-level tree: Root -> Light conditions (Daylight/Darkness) -> Road type
#' (5 types) -> Junction detail (7 types). 68 leaf categories.
#'
#' @seealso \code{\link{london_collisions}}
#'
#' @examples
#' data(london_tree)
#' head(london_tree)
"london_tree"


#' London Borough Boundaries
#'
#' Simplified polygon boundaries for the 33 London boroughs.
#'
#' @format An \code{sf} data.frame with 33 rows and 3 columns:
#' \describe{
#'   \item{NAME}{Borough name.}
#'   \item{GSS_CODE}{ONS geography code.}
#'   \item{geometry}{MULTIPOLYGON geometry in WGS84 (EPSG:4326).}
#' }
#'
#' @source
#' London Datastore, Statistical GIS Boundary Files for London.
#'
#' @seealso \code{\link{london_collisions}}
#'
#' @examples
#' data(london_boroughs_map)
#' nrow(london_boroughs_map)
"london_boroughs_map"


# ==== Chicago ================================================================

#' Chicago Crime (2023)
#'
#' Combined long-format dataset of crime incidents in the 77 community areas
#' of Chicago, 2023. Use together with \code{\link{chicago_tree}}.
#'
#' @format A data.frame in long format with columns:
#' \describe{
#'   \item{region_id}{Integer identifier of the community area.}
#'   \item{area_number}{Official community area number (1--77).}
#'   \item{name}{Community area name.}
#'   \item{population}{Total incidents in the area. This is a compositional
#'     denominator (sums to the total incident count of the city) and is
#'     useful for analyses that ask "which crime types over-occur in which
#'     areas relative to the citywide profile". Not a population at risk.}
#'   \item{pop_residential}{Residential population of the community area.
#'     Use this as the denominator for incidence-rate analyses (incidents
#'     per resident). Source: U.S. Census Bureau, ACS 2020 5-year estimates,
#'     aggregated to community areas by CMAP Community Data Snapshots.
#'     The 77 values sum to approximately 2.71M, matching the published
#'     Chicago population.}
#'   \item{x}{Longitude of the centroid.}
#'   \item{y}{Latitude of the centroid.}
#'   \item{node_id}{Character leaf identifier from \code{\link{chicago_tree}}.}
#'   \item{cases}{Integer count of crime incidents.}
#' }
#'
#' @seealso \code{\link{chicago_tree}}, \code{\link{chicago_map}}
#'
#' @source
#' Crime incidents: City of Chicago Data Portal, Crimes -- 2001 to Present.
#'
#' Residential population: U.S. Census Bureau ACS 2020 5-year estimates,
#' aggregated to community areas by the Chicago Metropolitan Agency for
#' Planning (CMAP), Community Data Snapshots (2023 release).
#' \url{https://cmap.illinois.gov/data/community-data-snapshots/}
#'
#' @examples
#' data(chicago_crimes)
#' head(chicago_crimes)
#' cat("Total incidents:", sum(chicago_crimes$cases), "\n")
#' cat("Total residents:",
#'     sum(unique(chicago_crimes[, c("area_number", "pop_residential")])$pop_residential),
#'     "\n")
"chicago_crimes"


#' Chicago Crime - Hierarchy
#'
#' Hierarchical classification of crime incidents by type, description and
#' location group.
#'
#' @format A data.frame with 2,841 rows and 2 columns:
#' \describe{
#'   \item{node_id}{Character identifier of each node.}
#'   \item{parent_id}{Character identifier of the parent node. \code{NA}
#'     for the root.}
#' }
#'
#' @details
#' 4-level tree: Root -> Crime Type -> Description -> Location Group.
#' 2,486 leaf nodes. Leaf node IDs follow the pattern
#' \code{"<TYPE> | <DESCRIPTION> | <LOCATION>"}.
#'
#' @seealso \code{\link{chicago_crimes}}
#'
#' @examples
#' data(chicago_tree)
#' head(chicago_tree)
"chicago_tree"


#' Chicago Community Area Boundaries
#'
#' Simplified polygon boundaries for the 77 Chicago community areas.
#'
#' @format An \code{sf} data.frame with 77 rows and 3 columns:
#' \describe{
#'   \item{NAME}{Community area name.}
#'   \item{AREA_NUM}{Community area number (integer, 1--77).}
#'   \item{geometry}{MULTIPOLYGON geometry in WGS84 (EPSG:4326).}
#' }
#'
#' @source
#' City of Chicago Data Portal, Community Area Boundaries.
#'
#' @seealso \code{\link{chicago_crimes}}
#'
#' @examples
#' data(chicago_map)
#' nrow(chicago_map)
"chicago_map"


# ==== Florida (raw) ==========================================================

#' Florida General Mortality Data (2016)
#'
#' Death counts by county and ICD-10 cause of death for the state of Florida,
#' USA, 2016. This is \strong{raw data}: the user builds the ICD-10 tree
#' and renames columns to match treeSS conventions in the analysis script.
#'
#' @format A data.frame with 3,066 rows and 6 columns:
#' \describe{
#'   \item{county_fips}{5-digit FIPS code (character).}
#'   \item{county_name}{County name.}
#'   \item{icd10_code}{ICD-10 cause of death code.}
#'   \item{icd10_desc}{Description of the ICD-10 code.}
#'   \item{deaths}{Number of deaths (integer).}
#'   \item{population}{County population estimate (integer).}
#' }
#'
#' @details
#' To use with \code{\link{treespatial_scan}}, rename and merge with
#' centroids:
#' \preformatted{
#'   # Centroids from tigris
#'   data <- transform(fl_deaths,
#'                     region_id = county_fips,
#'                     node_id   = icd10_code,
#'                     cases     = deaths)
#'   # Add x, y from a centroids table
#'   data <- merge(data, centroids, by = "county_fips")
#'   treespatial_scan(data, fl_tree, ...)
#' }
#'
#' Covers 65 of Florida's 67 counties (2 small counties excluded by CDC
#' suppression rules), 253 ICD-10 codes, and 157,000 total deaths.
#'
#' @source
#' CDC WONDER Compressed Mortality File 1999--2016
#' (\url{https://wonder.cdc.gov/cmf-icd10.html}).
#'
#' @examples
#' data(fl_deaths)
#' head(fl_deaths)
"fl_deaths"
