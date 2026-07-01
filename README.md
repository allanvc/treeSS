# treeSS

**Tree-Spatial Scan Statistic for Cluster Detection**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/treeSS)](https://CRAN.R-project.org/package=treeSS)
[![R-CMD-check](https://github.com/allanvc/treeSS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/allanvc/treeSS/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/allanvc/treeSS/graph/badge.svg)](https://app.codecov.io/gh/allanvc/treeSS)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Implements the tree-spatial scan statistic (Cançado et al., 2025), which
detects clusters that are anomalous in both geographic space and a
hierarchical tree simultaneously. The method searches over circular spatial
zones and branches of a classification tree to find regions where observed
cases significantly exceed expectations under a Poisson or binomial model,
selectable via the `model` argument.

## Installation

```r
# From CRAN (after acceptance)
install.packages("treeSS")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("allanvc/treeSS")
```

## Quick start

```r
library(treeSS)

# Example: London road collisions
data(london_collisions)
data(london_tree)

# The scan functions take a data.frame as the first argument and refer to
# its columns by name. This keeps the choice of denominator,
# coordinates, etc. transparent.
result <- treespatial_scan(
  london_collisions,
  cases       = cases,
  population  = population,
  region_id   = region_id,
  x           = x,
  y           = y,
  node_id     = node_id,
  tree        = london_tree,
  nsim        = 999, seed = 42,
  n_cores     = 4L                 # parallelize the MC over 4 threads
)
print(result)

# Extract cluster membership for visualization
cr <- get_cluster_regions(result, n_clusters = 3, overlap = FALSE)
```

## Included datasets

| Dataset | Country | Domain | Regions | Tree |
|:--------|:--------|:-------|:--------|:-----|
| `rj_mortality` + `rj_tree` | Brazil | Infant mortality | 92 municipalities | ICD-10 (622 nodes) |
| `fl_deaths` | USA | General mortality | 65 counties | raw (built by user) |
| `london_collisions` + `london_tree` | UK | Road collisions | 33 boroughs | Light x Road x Junction (81 nodes) |
| `chicago_crimes` + `chicago_tree` | USA | Crime | 77 community areas | Type x Description x Location (2841 nodes) |
| `rj_map`, `london_boroughs_map`, `chicago_map` | Brazil / UK / USA | Polygon boundaries | 92 / 33 / 77 | -- |

## Key functions

- `treespatial_scan()` — tree-spatial scan (main function)
- `circular_scan()` — Kulldorff's spatial scan
- `tree_scan()` — tree-based scan
- `filter_clusters()` — non-overlapping secondary clusters (Cançado et al. 2025)
- `sequential_scan()` — sequential adjustment for secondary clusters (Zhang, Assunção & Kulldorff 2010)
- `get_cluster_regions()` — cluster membership for any visualization package

## Visualization

The package is visualization-agnostic. `get_cluster_regions()` returns a
data.frame that can be merged with any spatial object for plotting with
ggplot2, leaflet, tmap, or any other mapping package. The bundled `sf`
boundary datasets (`rj_map`, `london_boroughs_map`, `chicago_map`) let the
examples map clusters without any external boundary download. See
`vignette("introduction")` for worked examples with ggplot2 + `rj_map`
(Brazil), leaflet + tigris (USA), and leaflet + `london_boroughs_map`
(London).

## References

Cançado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
(2025). A tree-spatial scan statistic. *Environmental and Ecological
Statistics*, 32, 953–978. doi:10.1007/s10651-025-00670-w
