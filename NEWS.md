# treeSS 0.2.6-1

## Documentation

* `README.md` now also carries the CRAN downloads badges (monthly and
  grand-total downloads from the RStudio CRAN mirror, via cranlogs).

# treeSS 0.2.6

* Updated package authorship in `DESCRIPTION`: the package code authors are
  now Allan Quadros and Andre L. F. Cançado. The citation for the underlying
  methodology paper (Cançado, Oliveira, Quadros & Duczmal, 2025) is unchanged.

# treeSS 0.2.5

## Performance

* Zone construction is dramatically faster for larger study areas.
  `build_zones()` previously grew each zone's membership vector and the
  zones list incrementally inside a double loop, causing roughly cubic
  cost in the number of regions; it now precomputes each center's
  inclusion sequence in a single pass and materializes all zones as
  prefix subsets into a preallocated list. On a 300-region study area the
  zone-construction step drops from about 19 s to about 0.14 s
  (~130x), and a full `circular_scan()` at `nsim = 999` from about 20 s
  to about 1.3 s (~15x), measured single-threaded.

* `.zones_to_csr()` no longer concatenates the zone membership index
  incrementally (which was quadratic in the total number of
  zone-region memberships); it now computes pointers with `cumsum()`
  and flattens memberships with a single `unlist()`.

* `.tree_to_csr_children()` replaces the per-node `which()` scan
  (quadratic in the number of tree nodes) with a vectorized
  `match()`/`order()`/`tabulate()` grouping, which is effectively linear.
  This reduces preprocessing time for large classification trees such as
  the 2841-node Chicago crime taxonomy.

* None of these changes affect results: the zone family (including the
  population-cap skip behaviour), the CSR structures passed to the C++
  backend, and therefore all clusters, log-likelihood ratios, and Monte
  Carlo p-values are bit-for-bit identical to 0.2.4 for the same inputs
  and seeds. Verified by snapshot comparison on the Rio de Janeiro,
  Chicago, and London analyses and by `identical()` checks of the old and
  new constructors on synthetic inputs, including cap-skip edge cases.

# treeSS 0.2.4

## Documentation

* `README.md` now carries the standard status badges (CRAN version,
  R-CMD-check, Codecov test coverage, lifecycle, and license).

* The README "Included datasets" table lists the `rj_map` dataset (added
  in 0.2.2) alongside `london_boroughs_map` and `chicago_map` in the
  polygon-boundaries row, and the "Visualization" section now points to the
  bundled `sf` boundary datasets (`rj_map`, `london_boroughs_map`,
  `chicago_map`) instead of the removed `geobr` download for the Brazil
  example.

# treeSS 0.2.3

## Examples and vignettes

* `example_brazil_rj.R` and the `introduction` vignette now map clusters using
  the bundled `rj_map` dataset (added in 0.2.2) instead of downloading
  boundaries with `geobr`, so the Rio de Janeiro material runs without
  `geobr`/`arrow`, consistent with the Chicago and London examples. All
  example scripts use the data-frame/column scan interface.

# treeSS 0.2.2

## New data

* Added the `rj_map` dataset: `sf` polygon boundaries for the 92
  municipalities of Rio de Janeiro state (IBGE Malhas Municipais), with an
  `ibge_code` join key to `rj_mortality`. The Rio de Janeiro example now maps
  clusters with `data(rj_map)`, mirroring `chicago_map` and
  `london_boroughs_map`, so it no longer requires an external polygon
  download (previously `geobr::read_municipality()`).

# treeSS 0.2.1

## Programmatic column selection

* The column arguments of the scan functions now also accept a column name
  held in a variable, e.g. `pop <- "live_births"; treespatial_scan(data,
  population = pop, ...)`. Previously only a bare name (`population`) or a
  literal string (`"live_births"`) worked; a variable holding the name was
  mis-resolved as a length-1 value, which failed with a confusing
  "must have the same length as 'cases'" error. This makes programmatic use
  (looping over datasets/denominators) work as expected. Bare names, literal
  strings, and expressions over columns continue to work unchanged.

# treeSS 0.2.0

## Breaking change: data/column interface

The scan functions no longer take parallel vectors. Each now takes a
**`data` data.frame as its first argument** and refers to its columns by
(unquoted) name. This makes calls shorter, pipe-friendly (`data |>
treespatial_scan(...)`), and removes the repeated `df$column` boilerplate.

* `treespatial_scan()`, `circular_scan()`, `sequential_scan()` and
  `aggregate_tree()` gained a leading `data` argument; their `cases`,
  `population`, `region_id`, `x`, `y` and `node_id` arguments now name
  columns of `data` rather than being vectors. Column arguments accept an
  unquoted name (`cases`), a string (`"cases"`), or an expression on
  columns (`raw_count * weight`).

* `tree_scan()` is now keyed by `node_id`: it takes
  `tree_scan(data, cases, node_id, tree, population = NULL)` with one row
  per leaf. Rows are matched to the tree by `node_id` (so they no longer
  need to be pre-ordered to the tree's leaf order), and counts are summed
  within leaf.

* `sequential_scan()` tree-only mode now also requires a `node_id` column
  (consistent with `tree_scan()`); previously it relied on the row order.

* **Migration.** Wrap your vectors in a `data.frame` and pass column names:

  ```r
  # before (<= 0.1.x)
  treespatial_scan(cases = d$cases, population = d$population,
                   region_id = d$region_id, x = d$x, y = d$y,
                   node_id = d$node_id, tree = tree)

  # now (>= 0.2.0)
  treespatial_scan(d, cases, population, region_id, x, y, node_id,
                   tree = tree)
  ```

  The `tree` argument is unchanged (still a separate `node_id`/`parent_id`
  data.frame, or the `tree_node_id`/`tree_parent_id` vectors). The returned
  objects, their classes, and all `print`/`summary`/`filter_clusters()`/
  `get_cluster_regions()` behaviour are unchanged.

## Internal

* New internal helper `.resolve_col()` performs the (base-R, dependency-free)
  non-standard evaluation that maps column arguments to vectors. The C++
  Monte Carlo core and the statistical results are unchanged.

* Bundled example scripts (`inst/examples/`) and vignettes were updated to
  the new interface.

# treeSS 0.1.52

## Documentation

* All help pages are now generated from roxygen2 comments. The twelve
  `.Rd` files that were previously maintained by hand (the three map/raw
  datasets `chicago_map`, `london_boroughs_map`, `fl_deaths`,
  `filter_clusters()`, and the eight `print`/`summary` methods) have been
  consolidated into the roxygen blocks in `R/data.R`, `R/filter_clusters.R`,
  and `R/print.R`, so `devtools::document()` no longer skips them. Content
  from the hand-written pages (the `st_simplify` note and merge tip for
  `london_boroughs_map`, the Crown-copyright attribution, and the fuller
  `fl_deaths` examples) was preserved. No user-visible change to the
  rendered documentation.

# treeSS 0.1.51

## Unified Monte Carlo implementation across `n_cores`

The three Monte Carlo routines (`mc_treespatial_cpp`, `mc_spatial_cpp`,
`mc_treescan_cpp`) now use a **single native C++ implementation for every
`n_cores >= 1`**. Previously `n_cores = 1` took a separate code path that
used R's `rmultinom()` over `NumericMatrix` objects, while `n_cores > 1`
used a native `std::mt19937` sampler over flat arrays.

* **Results are now invariant to `n_cores`.** Each simulation draws from a
  deterministic per-simulation seed (from R's RNG when `seed` is set), so
  the simulated null distribution and the resulting `p`-value are identical
  for any thread count given a fixed `seed`. `n_cores` changes only
  wall-clock time.
* **Serial and parallel timings are now directly comparable**, because the
  only difference between them is the number of threads, not the algorithm
  or the RNG.
* **Behaviour change:** simulated `p`-values at `n_cores = 1` are no longer
  bit-identical to the pre-0.1.50 serial path (which used R's `rmultinom`).
  Observed statistics, most-likely clusters, and secondary-cluster
  extraction are unaffected. Fix your `seed` to reproduce results.
* Removed the now-unused serial-only helpers `aggregate_up()` and
  `max_llr_all_pairs()`.

## Examples

* `example_chicago.R` now uses the compositional `population` denominator
  (total incidents per area) rather than `pop_residential`. With the
  residential denominator the most likely cluster is a broad-spectrum
  spatial hotspot reported at the tree root; the compositional denominator
  asks which (crime category, area) combinations are over-represented and
  returns a specific branch, which is the tree-spatial use the method is
  designed for.

# treeSS 0.1.50 (2nd CRAN patch)

Small adjustments to the DESCRIPTION file

# treeSS 0.1.49

Small adjustments to the vignettes

# treeSS 0.1.48

## Vignettes restructured

The package now ships two vignettes:

* `vignette("introduction", package = "treeSS")` — Rio de Janeiro
  end-to-end, reproducing Section 5.2 of Cançado et al. (2025). This
  was the previous `introduction` vignette, trimmed to RJ only.

* `vignette("florida", package = "treeSS")` *(new)* — a pedagogical
  walk-through of building the tree-spatial scan inputs from raw data
  using the bundled `fl_deaths` dataset: building the ICD-10 tree
  from the codes that actually appear in the data, downloading
  county polygons + centroids from `tigris`, and assembling the
  parallel-vector input contract that `treespatial_scan()` expects.

The Chicago and London datasets, previously discussed inline in the
`introduction` vignette, are now reserved for the companion software
paper.


# treeSS 0.1.47

## Bug fix: spurious empty facet in sequential-scan map examples

The four bundled plotting examples for `sequential_scan()`
(`example_brazil_rj.R`, `example_chicago.R`, `example_florida.R`)
previously did a left join from the full map polygon set onto the
cluster table. When the shapefile contained polygons not present in
the analysis dataset (3 RJ municipalities missing from the
DATASUS/IBGE 89-municipality subset, for instance), those polygons
emerged with `panel = NA`, which `facet_wrap` rendered as an extra
empty panel labelled "NA".

The examples now cross-join the polygon set with the panel labels
first and then left-join the cluster information by `(id, panel)`,
so every map polygon is drawn in every iteration panel — those that
fall outside the analysis dataset get the `na.value` colour (a
light grey), exactly as intended. No extra "NA" panel is produced.

The `london` example uses `leaflet` rather than `facet_wrap` and was
not affected.


# treeSS 0.1.46

## Removed: `multicluster_scan()`

`multicluster_scan()` (added in 0.1.45 as an adaptation of Li, Wang,
Yang, Li and Lai 2011 to the tree-spatial setting) has been removed.
The function is gone, along with its C++ backend
(`mc_multicluster_treespatial_cpp`, `mc_multicluster_spatial_cpp`),
the `get_cluster_regions.multicluster_scan` S3 method, the
corresponding `print` / `summary` methods, all examples, and the
vignette subsection.

Rationale:

* On real datasets with a concentrated signal (e.g. infant mortality
  in Rio de Janeiro: 622 tree nodes, 5358 zones), the top-K
  candidate pool was dominated by overlapping variants of a single
  geographic neighbourhood, so the fast top-K disjoint-pair search
  could not find a valid pair. The full-pool rescue path was too
  slow to be practical (timing out on `nsim = 999` with 4 cores).

* The factorisation of the joint LLR used by Li et al. (2011) is
  exact under the Poisson model for circular scans; its extension to
  the tree-spatial setting was not formally established.

* `filter_clusters()` (Cançado et al. 2025) and `sequential_scan()`
  (Zhang, Assunção and Kulldorff 2010) together already cover the
  practical secondary-cluster use cases with published, well-studied
  statistical properties.

Users who want joint-cluster detection in the circular case can use
the original implementation from Li et al. (2011) outside this
package.

## Secondary-cluster methods after 0.1.46

The package now offers two clearly-bounded approaches:

* `filter_clusters()` — paper-faithful non-overlap criterion of
  Cançado et al. (2025), Sec. 5.1.1, applied to the single-pass
  candidate pool.

* `sequential_scan()` — sequential adjustment of Zhang, Assunção and
  Kulldorff (2010): detect MLC, remove its regions (with optional
  buffer of nearest neighbours), re-run the scan on the reduced data
  with a fresh Monte Carlo simulation; iterate until the current MLC
  is no longer significant. Each iteration's p-value is correct under
  the conditional argument in the paper, so no multiple-testing
  correction is required.


# treeSS 0.1.45

## Secondary clusters: methods overhaul

Replaced the ad-hoc Holm-Bonferroni `iterative_scan()` with two methods
drawn directly from the published literature on multi-cluster spatial
scan statistics, adapted to the tree-spatial setting. The package now
offers three approaches to secondary-cluster detection, with the
choice driven by which type of shadowing the user wants to remove:

* `filter_clusters()` (unchanged) -- the original non-overlap
  criterion of Cancado et al. (2025) Sec. 5.1.1, applied to the
  single-pass candidate pool.

* `sequential_scan()` (new) -- the sequential adjustment of Zhang,
  Assuncao and Kulldorff (2010), adapted to tree-spatial / circular /
  tree-only inputs. Detects the MLC, removes its regions (and an
  optional `buffer_size` of nearest neighbours) from the dataset, and
  re-runs the scan on the reduced data with a fresh Monte Carlo
  simulation. Iterates until the MLC of the current reduced data is
  no longer significant or `max_iter` is reached. Each iteration's
  p-value is correct under the conditional argument of Section 3 of
  the paper -- no post-hoc multiple-testing correction is applied or
  required.

* `multicluster_scan()` (new) -- the two-cluster joint statistic of
  Li, Wang, Yang, Li and Lai (2011), adapted to tree-spatial and
  circular scans. Builds the alternative as a joint presence of two
  region-disjoint clusters; the joint LLR factorises into the sum of
  the two single-cluster LLRs under Poisson, so the observed maximum
  is found by sweeping the candidate pool. The Monte Carlo for the
  joint statistic runs in C++ (new exports `mc_multicluster_treespatial_cpp`
  and `mc_multicluster_spatial_cpp`) with the same OpenMP backend as
  the other scans, so performance is on par with `treespatial_scan()`.
  The decision rule of Table 2 of the paper is applied: 0, 1, or 2
  significant clusters are reported based on the joint p-value and a
  re-evaluation of the weaker cluster on the reduced dataset.

## Removed

* `iterative_scan()` and its print/summary/`get_cluster_regions`
  methods have been removed. The Holm-Bonferroni "scan + zero cases +
  re-scan" procedure is not part of the published methods we wanted to
  offer; the sequential and multi-cluster scans above cover the
  intended use cases and are grounded in the literature.

* Internal helper `.matrix_to_vectors()` (previously used only by
  `iterative_scan`) has been removed.

## New S3 methods

* `print.sequential_scan()`, `summary.sequential_scan()`
* `print.multicluster_scan()`, `summary.multicluster_scan()`
* `get_cluster_regions.sequential_scan()`,
  `get_cluster_regions.multicluster_scan()`

## Documentation

* `filter_clusters()`, `treespatial_scan()`, and `circular_scan()`
  cross-reference the new methods in `@seealso`.
* The introduction vignette now closes with a section showing all
  three secondary-cluster approaches side by side.
* The four worked examples under `inst/examples/` (Brazil/RJ,
  Chicago, Florida, London) use `sequential_scan()` in place of the
  removed `iterative_scan()` block.

## Tests

* New `tests/testthat/test-sequential-scan.R` covering structure,
  the `max_iter` stopping rule, the buffer mechanism, behaviour
  under H0, and printing.
* New `tests/testthat/test-multicluster-scan.R` covering structure,
  the stronger-versus-weaker ordering, region disjointness of the
  returned pair, the significance decision rule, and printing.
* `tests/testthat/test-get-cluster-regions.R` and
  `tests/testthat/test-binomial.R` updated to drop their references
  to `iterative_scan()`.


# treeSS 0.1.44

## CRAN reviewer feedback

Address the four items requested in the first-round CRAN review.

## DESCRIPTION

* Single-quote software/API names per the CRAN cookbook: `OpenMP` is
  now written as `'OpenMP'` in the package description.
  Reference:
  <https://contributor.r-project.org/cran-cookbook/description_issues.html#formatting-software-names>

* Add DOI links to the two references that were previously cited
  without a link, using the CRAN-mandated
  `authors (year) <doi:...>` form (no space after `doi:`, no space
  inside the angle brackets):
    - Kulldorff (1997) <doi:10.1080/03610929708831995>
    - Kulldorff et al. (2003) <doi:10.1111/1541-0420.00039>
  Reference:
  <https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

## Documentation (R/print.R, R/iterative_scan.R, man/*.Rd)

* Added missing `\value` tags (and the corresponding `@return`
  roxygen blocks) to the seven `print()`/`summary()` method Rd
  files flagged by CRAN. Each documents that the method invisibly
  returns its input object unchanged and is called for its
  printing side effect, with a description of the fields written
  to the console (and, for `summary()` methods, the additional
  fields beyond those of the matching `print()` method):
    - `print.circular_scan.Rd`
    - `print.iterative_scan.Rd`
    - `print.tree_scan.Rd`
    - `print.treespatial_scan.Rd`
    - `summary.circular_scan.Rd`
    - `summary.tree_scan.Rd`
    - `summary.treespatial_scan.Rd`
  Reference:
  <https://contributor.r-project.org/cran-cookbook/docs_issues.html#missing-value-tags-in-.rd-files>

## Bug fixes (R/generate_example_data.R, man/generate_example_data.Rd)

* `generate_example_data()` no longer sets a hardcoded seed within
  the function: the default of the `seed` argument is now `NULL`
  (previously `123L`). When the user does not pass a seed, the
  function draws from the user's session-level RNG state without
  modifying it; when the user passes an explicit integer, the
  existing save-and-restore logic (introduced in 0.1.43) still
  applies. The `\usage{}` block and the `\item{seed}{...}`
  description of the corresponding Rd file have been updated to
  match. The roxygen example
  (`ex <- generate_example_data(seed = 42)`) is unchanged: it
  passes an explicit seed and so remains reproducible.
  Reference:
  <https://contributor.r-project.org/cran-cookbook/code_issues.html#setting-a-specific-seed>


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
* Monte Carlo simulation for p-value computation, with the null resampler matched to the chosen model (multinomial conditional
resampling for Poisson; binomial conditional resampling for binomial).
