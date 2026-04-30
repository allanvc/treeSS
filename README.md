# treeSS

**Tree-Spatial Scan Statistic for Cluster Detection**

Implements the tree-spatial scan statistic (Cançado et al., 2025), which
detects clusters that are anomalous in both geographic space and a
hierarchical tree simultaneously. The method searches over circular spatial
zones and branches of a classification tree to find regions where observed
cases significantly exceed expectations under a Poisson model.

## Installation

```r
# From source
install.packages("treeSS_0.1.6.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(treeSS)

# Example: London road collisions
data(london_collisions)
data(london_tree)

# The function takes parallel vectors - extract them from your data.frame
# and pass each one explicitly. This makes the choice of denominator,
# coordinates, etc. transparent.
result <- treespatial_scan(
  cases       = london_collisions$cases,
  population  = london_collisions$population,
  region_id   = london_collisions$region_id,
  x           = london_collisions$x,
  y           = london_collisions$y,
  node_id     = london_collisions$node_id,
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
| `london_boroughs_map`, `chicago_map` | | Polygon boundaries | | -- |

## Key functions

- `treespatial_scan()` — tree-spatial scan (main function)
- `circular_scan()` — Kulldorff's spatial scan
- `tree_scan()` — tree-based scan
- `filter_clusters()` — non-overlapping secondary clusters
- `get_cluster_regions()` — cluster membership for any visualization package

## Visualization

The package is visualization-agnostic. `get_cluster_regions()` returns a
data.frame that can be merged with any spatial object for plotting with
ggplot2, leaflet, tmap, or any other mapping package. See `vignette("introduction")`
for worked examples with ggplot2 + geobr (Brazil), leaflet + tigris (USA),
and leaflet + sf (London).

## References

Cançado, A. L. F., Oliveira, G. S., Quadros, A. V. C., & Duczmal, L.
(2025). A tree-spatial scan statistic. *Environmental and Ecological
Statistics*, 32, 953–978. doi:10.1007/s10651-025-00670-w
