## Test environments
* local Ubuntu 22.04.5 LTS, R 4.4.1
* win-builder, R-devel (2026-07-10 r90234 ucrt)
* R-hub v2: linux (R-devel), windows (R-devel), macos (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

* win-builder (R-devel): Status OK, including CRAN incoming feasibility.
* R-hub v2 (linux, windows, macos, all R-devel): all passing.
* The only NOTE observed was on the local machine
  ("checking for future file timestamps ... unable to verify current
  time"), an artifact of the local environment with no network access
  to a time server; it does not appear on win-builder or R-hub.

## Comments

* This version contains performance improvements only (zone
  construction and CSR preprocessing). Results are bit-for-bit
  identical to 0.2.4 for the same inputs and seeds.
