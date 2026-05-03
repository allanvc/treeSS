# treeSS 0.1.43

## CRAN-readiness pass

Testing a a clean `R CMD check --as-cran`.

## DESCRIPTION

- added ORCID from Andre


## Documentation (R/data.R, man/{chicago,london,rj}_tree.Rd)

* Added `\source{}` blocks to all three tree datasets, pointing at
  the corresponding leaf-level dataset and at the `data-raw/`
  build script in the GitHub repo.

## Documentation (R/get_cluster_regions.R, man/get_cluster_regions.Rd)

* Added a worked example to `get_cluster_regions()`. Added `@examples` block.

## Documentation (R/filter_clusters.R)

* Added an `@examples` block to the roxygen comments.

## Bug fixes (R/circular_scan.R, R/tree_scan.R, R/treespatial_scan.R, R/generate_example_data.R)

* The four functions that accept a `seed = ...` argument no longer
  silently overwrite the user's session-level RNG state. Previously,
  calling `treespatial_scan(..., seed = 42)` after a `set.seed(2026)`
  in the user's session would leave the RNG in a state determined by
  the internal Monte Carlo loop, so any subsequent `runif()`,
  `sample()`, etc. was no longer reproducible from the user's
  `set.seed(2026)`. Now the user's pre-existing RNG state is saved on
  entry and restored on exit (whether the function returns normally
  or via an error), so the `seed` argument affects only the result
  of the call. Implementation is in two new internal helpers
  `.seed_save_and_set()` and `.seed_restore()` in `R/utils.R`.

## Print methods (R/iterative_scan.R, man/print.iterative_scan.Rd)

* `print.iterative_scan()` now accepts `max_show` for API
  consistency with the other three print methods. The default
  behavior is unchanged (the table is printed without the
  `region_ids` and `leaf_ids` columns to keep it compact); pass
  `max_show = -1L` to include both columns.

## CRAN submission infrastructure

* Added `cran-comments.md` file.

## README.md

* Updated the install snippet to a CRAN install + a development install
  via `remotes::install_github("allanvc/treeSS")`.

# treeSS 0.1.42

## Documentation (R/print.R, man/*.Rd)

* The `summary()` methods for `circular_scan`, `tree_scan`, and
  `treespatial_scan` now have proper roxygen descriptions and
  explicitly document that the `max_show` argument added in 0.1.39
  is forwarded to the corresponding `print()` method via
  \code{...}. Each summary doc points to the matching print doc
  for the full details.


# treeSS 0.1.39

## Print methods (R/print.R) - truncation by default

The print methods now truncate long `Leaf IDs` and `Regions` lists
by default, in the style of `tibble`. The motivation is the Chicago
example: the most likely cluster turns out to be
the **root** of the FBI crime taxonomy (1900+ leaves), which under
the previous policy printed every single leaf, producing more than
10 pages of console output in the rendered PDF.

* New argument `max_show` on `print.treespatial_scan()`,
  `print.tree_scan()` and `print.circular_scan()`. Default is
  `10L`. When a vector field exceeds this length, only the first
  `max_show` values are shown and a tail of `... and N more` is
  appended. Pass `max_show = -1L` (or any value at least as large
  as the field) to recover the previous full-output behavior.

* The internal `.cat_wrapped()` helper gained the same `max_show`
  argument (default `10L`) and propagates it through the print
  methods.

* No changes to the underlying scan results: only the console / PDF
  rendering of the result objects is affected. The full leaf and
  region IDs are always available on `result$most_likely_cluster$
  leaf_ids` and `result$most_likely_cluster$region_ids` for
  programmatic use.

The choice of default mirrors `tibble`'s behavior: enough to give
the reader a sense of the cluster contents, but not so much that
a single `print()` call dominates the document.


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
