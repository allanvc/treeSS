## Test environments
* local Ubuntu 24.04, R 4.4.2
* win-builder (release and devel)
* macOS builder (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission.

## Notes on the auto-check

The CRAN incoming feasibility check flags the following words in the
DESCRIPTION as possibly misspelled. All are intentional:

* `Cancado`, `Kulldorff`, `Kulldorff's` -- surnames of the authors of
  the cited methodology papers (Cancado et al. 2025; Kulldorff 1997;
  Kulldorff et al. 2003).
* `OpenMP` -- the standard name of the parallelization API used
  optionally by the package.
* `al`, `et` -- components of the standard bibliographic abbreviation
  "et al.".
