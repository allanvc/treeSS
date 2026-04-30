# treeSS 0.1.18

* Initial release.
* Implements the tree-spatial scan statistic (Cancado et al., 2025).
* Provides `treespatial_scan()` for combined spatial and hierarchical
  cluster detection.
* Provides `circular_scan()` for Kulldorff's circular spatial scan
  statistic.
* Provides `tree_scan()` for the tree-based scan statistic.
* Helper functions: `build_zones()`, `aggregate_tree()`,
  `filter_clusters()`.
* S3 `print()` and `summary()` methods for all scan result classes.
* Monte Carlo simulation for p-value computation using the Poisson model.
