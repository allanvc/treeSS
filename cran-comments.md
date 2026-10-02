## Test environments
* local Ubuntu 22.04.5 LTS, R 4.4.1 (R CMD check --as-cran on the
  tarball built from a clean export, 2026-10-02)
* win-builder, R-devel: to be run before submission
* R-hub v2: linux, windows, macos (R-devel): to be run before submission

## R CMD check results

0 errors | 0 warnings | 1 note

* The local check shows only the usual local-machine NOTE
  ("checking for future file timestamps ... unable to verify current
  time"), an artifact of an environment with no network access to a
  time server.
* 125 unit tests pass (testthat), 2 skipped (optional dependencies).

## Comments

* This release adds `plot()` methods for all scan results
  (`plot.treess()` and `plot.tree_scan()`), a common parent class
  `"treess"` with accessor generics (`most_likely_cluster()`,
  `secondary_clusters()`, `pvalue()`), and print improvements. The
  scan algorithms and their results are unchanged from 0.2.6.
* `ggplot2` is added to Suggests; it is only used by `plot.treess()`
  when an `sf` polygon layer is supplied, and all examples and tests
  that need it are conditional on its availability.
* This version is the one described in a manuscript submitted to
  The R Journal, whose code depends on the new `plot()` methods.
