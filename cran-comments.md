## Test environments
* local Ubuntu 22.04.5 LTS, R 4.4.1
* win-builder, R-devel (2026-07-26 r90304 ucrt)
* R-hub v2: linux (R-devel), windows (R-devel), macos (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* win-builder (R-devel): 1 NOTE on CRAN incoming feasibility, flagging
  the URL <https://opendatasus.saude.gov.br> in man/rj_mortality.Rd as
  possibly invalid ("server closed abruptly"). This is the canonical
  URL of the Brazilian Ministry of Health open-data portal, which is
  the source of the rj_mortality dataset; the portal was temporarily
  unavailable at check time (it currently returns server errors for
  all pages). The same URL passed the win-builder incoming checks for
  version 0.2.5 two weeks ago.
* R-hub v2 (linux, windows, macos, all R-devel): all passing.
* The local check showed only the usual local-machine NOTE
  ("checking for future file timestamps ... unable to verify current
  time"), an artifact of an environment with no network access to a
  time server.

## Comments

* This is a metadata-only update shortly after 0.2.5 was accepted:
  it corrects the package authorship in DESCRIPTION (the package code
  authors are Allan Quadros and Andre L. F. Cançado). There are no
  code changes; results are identical to 0.2.5. The citation for the
  underlying methodology paper is unchanged.
* Apologies for the quick resubmission; we wanted the authorship
  metadata on CRAN to be correct as early as possible.
