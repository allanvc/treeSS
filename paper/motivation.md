# Motivation letter

To: The Editors, The R Journal

Dear Editors,

We are pleased to submit the manuscript "*treeSS: Tree-Spatial Scan
Statistics for Cluster Detection in R*" for consideration by The R Journal.

The paper describes the **treeSS** package, the first R package on CRAN
implementing the tree-spatial scan statistic of Cançado, Oliveira, Quadros
and Duczmal (2025, *Environmental and Ecological Statistics*, 32, 953–978,
doi:10.1007/s10651-025-00670-w). The method generalises Kulldorff's
circular spatial scan by searching jointly over circular zones and
branches of a hierarchical classification tree, returning the (zone,
branch) pair maximising the likelihood ratio. Significance is assessed by
Monte Carlo simulation under either a Poisson or a binomial model.

We believe the paper is well suited to the R Journal for three reasons.

First, it fills a concrete gap in the R ecosystem. Existing packages
cover spatial-only scans (`SpatialEpi`), space-time scans
(`scanstatistics`), or wrap external binaries (`rsatscan` for SaTScan).
None of them implements a scan that searches jointly over space and a
classification tree. The reference standalone tools, SaTScan and
TreeScan, do not implement this combined scan either.

Second, the package is supported by three real datasets shipped with it
-- infant mortality in Rio de Janeiro state (Brazil), reported crimes in
Chicago, and road traffic collisions in London -- which make it
straightforward for readers to run the full workflow end-to-end. The Rio
de Janeiro example reproduces Table 7 of Cançado et al. (2025).

Third, the package follows R Journal best practices: snake_case naming,
S3 print/summary methods on every result class, unit-tested with
testthat, vignette-driven, GitHub-hosted, with a pkgdown site and a
Zenodo-archived release. The core search loop is implemented in C++ via
Rcpp and parallelised with OpenMP, with the serial path being
bit-identical to single-threaded execution for full reproducibility.

We confirm that the manuscript has not been published or submitted
elsewhere, that all code and data are open and freely available, and
that the included Rmarkdown source builds in well under ten minutes on a
modern multicore laptop.

This software paper is complementary to the methodological paper of
Cançado et al. (2025): the present submission focuses on package design,
API decisions, performance, and worked applications, deferring
methodological derivations to the original publication.

We would be happy to address any questions or comments from the
editorial board and reviewers.

Yours sincerely,

André L. F. Cançado (corresponding author for methodology)
Allan V. C. Quadros (corresponding author for software / package
maintainer)
Geiziane S. Oliveira
Luiz H. Duczmal
