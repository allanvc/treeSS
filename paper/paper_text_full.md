# treeSS: Tree-Spatial Scan Statistics for Cluster Detection in R

*This is a human-readable rendering of the full paper text, with YAML
header, code chunks and comments stripped. Citations of the form
`@key` and `[@key]` are placeholders that the `rjtools` toolchain
expands to "(Author, Year)" in the compiled PDF/HTML.*

---

## Authors (in order)

1. **Allan V. C. Quadros** — Department of Management, University of
   North Florida (corresponding author / maintainer)
2. **André L. F. Cançado** — University of Brasília
3. **Geiziane S. Oliveira** — Federal University of Minas Gerais
4. **Luiz H. Duczmal** — Federal University of Minas Gerais

## Abstract

We introduce **treeSS**, an R package implementing the tree-spatial scan
statistic of Cançado et al. (2025). The method
generalizes Kulldorff's circular spatial scan statistic by searching jointly
over circular spatial zones and branches of a hierarchical classification
tree, returning the (zone, branch) pair that maximizes the likelihood ratio
under either a Poisson or a binomial model. Statistical significance is
assessed via Monte Carlo simulation. The package targets epidemiological and
surveillance applications where the events to be clustered carry both a
geographic location and a categorical hierarchy, such as deaths classified
by ICD-10, crimes by FBI category, or road collisions by severity. No
previous CRAN package implements this combined scan; existing tools cover
either spatial-only scans (such as **SpatialEpi**), space-time scans (such
as **scanstatistics**), or wrappers around external binaries (such as
**rsatscan** for SaTScan). The core search is implemented in C++ via
**Rcpp** and optionally parallelized with OpenMP. We illustrate the package
on three real datasets shipped with it — infant mortality in Rio de
Janeiro state (Brazil), reported crimes in Chicago community areas, and
road traffic collisions across the 33 boroughs of London — and show that
the package reproduces the results reported in the original methodology
paper.

---

## 1. Introduction

Detecting unusual clusters of events in space is a long-standing problem in
public-health surveillance, criminology, ecology, and demography. The
canonical tool is Kulldorff's circular spatial scan statistic
(Kulldorff & Nagarwalla, 1995; Kulldorff, 1997), implemented in the
standalone SaTScan software (Kulldorff, 2024) and widely used to identify
localized excess risk in disease incidence, mortality, and crime rates.
Within R, the scan statistic is exposed through several packages of varying
scope: **SpatialEpi** (Kim & Wakefield, 2024) provides Kulldorff's circular
scan, **scanstatistics** (Allévius, 2018) focuses on space-time anomaly
detection, and **rsatscan** (Kleinman, 2024) offers a thin wrapper around
the SaTScan binary, requiring an external installation.

A complementary line of work treats the events themselves as carrying
*hierarchical* structure rather than spatial structure. The tree-based scan
statistic of Kulldorff, Fang, and Walsh (2003) was developed for
adverse-event surveillance, where events (such as drug-related diagnoses)
are classified into a tree (such as ICD chapters, subchapters, and leaf
codes), and the goal is to identify a *branch* of the tree showing excess
risk relative to its expectation. The reference implementation, TreeScan
(Kulldorff, 2024), is again a standalone binary; no CRAN package implements
this method directly.

In many applied settings the two structures coexist. Infant deaths are
indexed by both *municipality* (a spatial location) and *ICD-10 cause* (a
position in a disease tree). Reported crimes are indexed by both
*neighborhood* and *offense category* (a hierarchy of theft, assault,
property crime, and so on). Road traffic collisions are indexed by both
*borough* and *severity ladder*. A purely spatial scan would aggregate over
all causes and miss a cause-specific cluster confined to a region; a purely
tree-based scan would aggregate over the whole study area and miss the
geographic concentration of a particular cause. The cluster of interest, in
each case, is a (region set, tree branch) pair.

Cançado et al. (2025) introduced a *tree-spatial scan statistic* that
searches for exactly such pairs. For each circular spatial zone *z* in the
sense of Kulldorff (1997) and each branch *g* of a user-supplied tree, the
statistic computes the likelihood ratio comparing a model where events in
(*z*, *g*) follow an elevated Poisson (or binomial) rate against the null of
a uniform rate, and reports the pair maximizing this ratio. Significance is
assessed by Monte Carlo simulation under the null. The method was validated
on infant mortality data from the state of Rio de Janeiro, where it
identified a previously undocumented cluster of deaths from an ICD-10 leaf
code ("J18.9 — Other bacterial pneumonias") concentrated in the north of
the state.

This paper introduces **treeSS**, an R package implementing this method. To
our knowledge it is the first CRAN package providing a tree-spatial scan;
**scanstatistics** does not address the tree dimension, **rsatscan**
delegates to a binary that does not implement this combined scan, and the
standalone TreeScan software does not perform the spatial search. The
package additionally provides standalone implementations of the circular
spatial scan (Kulldorff, 1997) and the tree-based scan (Kulldorff et al.,
2003) under a unified interface, so that the three scan types — spatial-
only, tree-only, and tree-spatial — can be compared on the same data with
the same input contract. The likelihood-ratio search and the Monte Carlo
loop are implemented in C++ via **Rcpp** (Eddelbuettel & François, 2011)
and optionally parallelized with OpenMP (Dagum & Menon, 1998), yielding a
4–8× speedup on multicore machines. Three real datasets are shipped with
the package — infant mortality in Rio de Janeiro (DATASUS; IBGE), reported
crime in Chicago (City of Chicago), and road traffic collisions in London
(UK Department for Transport) — to allow users to run the full workflow
end-to-end without preparing data themselves.

The remainder of this paper is structured as follows. Section 2 summarizes
the methodology to the level needed for the rest of the paper, deferring
derivations and proofs to Cançado et al. (2025). Section 3 describes the
package design, with particular attention to the choice of input contract,
the class hierarchy of scan results, and the helper functions for cluster
extraction and visualization. Section 4 walks through a complete worked
example reproducing the results of Cançado et al. (2025) on the Rio de
Janeiro infant mortality data. Section 5 gives two further applications, on
the Chicago crime and London collisions datasets, to illustrate the breadth
of the method. Section 6 discusses implementation details, performance, and
limitations. Section 7 compares **treeSS** with related packages and
standalone tools. We conclude in Section 8.

---

## 2. Methodology overview

This section sketches just enough of the tree-spatial scan statistic to make
the rest of the paper self-contained. For full derivations, the choice of
significance procedure, and simulation studies, we refer the reader to
Cançado et al. (2025).

**Setup.** Let *i* = 1, ..., *m* index the regions of a study area, with
populations *n_i* and centroids (*x_i*, *y_i*). Let *T* be a rooted tree
whose leaves *ℓ* enumerate the categories that observed events fall into,
and whose internal nodes represent aggregations of those categories. We
denote by *L*(*g*) ⊆ {*ℓ*} the set of leaves in the subtree rooted at node
*g*. The data are integer counts *c_{i,ℓ}* of events occurring in region
*i* and classified into leaf *ℓ*. Writing *C* = Σ *c_{i,ℓ}* for the total
event count and *N* = Σ *n_i* for the total population, the count rolled up
to a tree branch *g* is *C_g* = Σ_*i* Σ_{*ℓ*∈*L*(*g*)} *c_{i,ℓ}*, and the
count restricted to both a zone *z* ⊆ {1,...,*m*} and a branch *g* is
*c_{z,g}* with population *n_z* = Σ_{*i*∈*z*} *n_i*.

**Spatial zones.** Following Kulldorff (1997), the family of candidate
zones 𝒵 consists of circular zones: for each region *i* (the center) and
each *k* = 1, 2, ..., the zone *z_{i,k}* is the union of the *k* regions
nearest to *i* in Euclidean distance between centroids, capped at a maximum
population fraction (default 50% of *N*). This produces of order O(*m*²)
candidate zones with limited combinatorial blow-up.

**Likelihood ratio.** Under a Poisson model, the log likelihood ratio for a
single (zone, branch) pair (*z*, *g*) is

    Λ(z, g) = c_{z,g} log(c_{z,g} / E[c_{z,g}])
            + (C_g - c_{z,g}) log((C_g - c_{z,g}) / (C_g - E[c_{z,g}]))

if *c_{z,g}* > E[*c_{z,g}*], and 0 otherwise, where E[*c_{z,g}*] = *C_g* ·
*n_z* / *N* is the expected count in (*z*, *g*) under uniform risk. Under a
binomial model — in which *n_i* counts trials rather than population at
risk — the corresponding expression replaces Poisson rates by binomial
proportions; we omit the formula for brevity (see Cançado et al. 2025, eq.
4). The tree-spatial scan statistic is

    Λ* = max over (z ∈ 𝒵, g ∈ T) of Λ(z, g),

and the *most likely cluster* is the maximizing pair.

**Inference.** No tractable null distribution is available for Λ*.
Following Kulldorff (1997), significance is assessed by Monte Carlo
simulation: under H₀, the observed leaf-by-region count matrix is
re-sampled from a multinomial distribution conditional on the row totals
*C_ℓ* (the leaf marginals) and on the population shares *n_i* / *N*, the
statistic is recomputed on each replicate, and the *p*-value of Λ* is the
rank of the observed value among the simulated maxima.

**Secondary clusters.** The most likely cluster is rarely the only signal
of interest. Cançado et al. (2025), Section 5.1.1, retains additional
(zone, branch) pairs in decreasing order of Λ, declaring a candidate
*distinct* from those already retained when its branch is neither an
ancestor nor a descendant of any retained branch, *or* its zone has no
region in common with any retained zone. Section 3 below describes how this
is exposed in **treeSS** through the `filter_clusters()` and
`get_cluster_regions()` helpers. We additionally provide a sequential
*conditional* procedure inspired by Zhang, Assunção and Kulldorff (2010) in
`iterative_scan()`, which is not part of the original paper and is offered
as a complementary tool for users who prefer Holm–Bonferroni-corrected
sequential testing.

The three scans — circular spatial, tree, tree-spatial — share this
template. The spatial-only case takes a degenerate single-leaf tree and
reduces to Kulldorff (1997). The tree-only case takes a degenerate
single-region geography and reduces to Kulldorff et al. (2003). The package
exposes all three via a unified interface, described next.

---

## 3. Package design

The interface of **treeSS** was shaped by three observations from the
applied workflows we have in mind. First, scan-statistic inputs are
heterogeneous: per-region counts and populations, per-region centroids,
per-event leaf labels, and a tree structure that lives outside the
case-level table. Second, the choice of denominator (resident population
vs. live births vs. number of trials) is a substantive scientific decision
in surveillance work, and the package should make that decision explicit
rather than infer it from column names. Third, users will mix the three
scan types — spatial-only, tree-only, and tree-spatial — on the same data,
so the three functions should share an input contract and a result schema.
The subsections below describe how these observations translate into
concrete API choices.

### Input contract: parallel vectors

A scan-statistic call has many distinct inputs: per-region case counts and
populations, per-region centroids, the leaf classification of each event,
the tree structure itself, the model, the number of simulations, and the
significance threshold. Three input designs were considered:

1. **A single `data.frame`** with reserved column names. This is the style
   used by SaTScan input files. It is concise but ambiguous when the user's
   data has columns with different names: the user must rename, or the
   package must guess.
2. **A formula interface**, in the style of `lm()`. This works well when
   inputs are scalars per row, but a tree-spatial scan needs tree
   structure as a side input that does not fit into a formula naturally.
3. **Parallel vectors**, one per logical role.

We chose the third. Every scan function in **treeSS** accepts the
case-level inputs as named, parallel vectors:

    treespatial_scan(
      cases, population, region_id, x, y, node_id,
      tree = ..., ...
    )

Each role is named explicitly. The user picks which column of their data is
the denominator (for instance `live_births` rather than the resident
`population`, when computing infant mortality rates). There is no implicit
column-name convention to memorise, and partial inputs are caught at call
time rather than mid-scan. The cost is verbosity: a typical call has six to
eight named arguments. We consider that an acceptable price for an
explicit, auditable interface, especially in surveillance contexts where
the choice of denominator is a substantive scientific decision.

The tree itself can be supplied either as a two-column `data.frame` with
columns `node_id` and `parent_id` (root having `parent_id = NA`), or as two
parallel vectors `tree_node_id` and `tree_parent_id`. The two forms are
interconvertible, so users can use whichever fits their preprocessing.

### Function hierarchy

The package exposes three scan functions, an iterative wrapper, two
post-processing helpers, and one tree-construction helper:

| Function | Purpose |
|---|---|
| `treespatial_scan()` | Tree-spatial scan; main entry point |
| `circular_scan()` | Kulldorff's circular spatial scan, no tree |
| `tree_scan()` | Tree-based scan, no spatial dimension |
| `iterative_scan()` | Sequential conditional scan with Holm–Bonferroni correction |
| `filter_clusters()` | Distinct secondary clusters per Cançado et al. (2025), Sec. 5.1.1 |
| `get_cluster_regions()` | Tidy region-by-cluster table for plotting |
| `aggregate_tree()` | Roll counts up the tree to internal nodes |
| `build_zones()` | Construct the family of circular zones |

The first three functions return objects of S3 classes `"treespatial_scan"`,
`"circular_scan"`, and `"tree_scan"`, respectively; `iterative_scan()`
returns an object of class `"iterative_scan"`. Every class has `print()`
and `summary()` methods, and `get_cluster_regions()` is an S3 generic with
a default method (single-pass scans) and a dedicated method for
`"iterative_scan"` objects, allowing the same downstream visualization
pipeline regardless of which scan was used.

### Choice of model: Poisson and binomial

All scan functions accept `model = c("poisson", "binomial")`, defaulting to
Poisson. The choice changes the form of the likelihood ratio (Section 2)
but not the search procedure or the Monte Carlo step. The Poisson model is
appropriate when the denominator is a *population at risk* and case rates
are small, which is the typical surveillance setting. The binomial model is
appropriate when the denominator is a *number of trials* (cases-plus-
controls), as in case-control studies of rare adverse events.

### Parallel Monte Carlo

The Monte Carlo step is the dominant cost: each replicate re-runs the full
O(|𝒵|·|*T*|) likelihood evaluation. All scan functions expose `n_cores`
(default `1L`); when `n_cores > 1`, replicates are distributed across
OpenMP threads in C++. Each thread carries its own `std::mt19937` generator
seeded deterministically from the user-supplied `seed`. Two reproducibility
properties hold:

* `n_cores = 1` is bit-identical to the pre-OpenMP serial implementation.
* `n_cores > 1` is reproducible for any fixed pair `(seed, n_cores)`, but
  differs slightly from the serial path because R's RNG cannot be shared
  across threads safely; the simulated null distribution comes from the
  thread-local Mersenne Twister streams.

For publication-critical work where the simulated *p*-values must be
bit-identical between machines, we recommend `n_cores = 1L`. For
exploratory work, `n_cores = 4` typically yields a 4–6× speedup (see
Section 6).

### Result objects

Each scan returns a list with a stable schema. For `treespatial_scan()`:

| Element | Description |
|---|---|
| `most_likely_cluster` | List with `node_id`, `region_ids`, `cases`, `expected`, `population`, `rr`, `llr` |
| `secondary_clusters` | Data.frame of all (zone, branch) pairs evaluated, sorted by LLR |
| `pvalue` | Monte Carlo *p*-value for the most likely cluster |
| `nsim`, `alpha`, `model` | Hyperparameters passed in by the user |
| `simulated_llr` | Vector of length `nsim` with the simulated null distribution |
| `regions`, `tree`, `full_cases` | Inputs preserved for downstream use |

The `secondary_clusters` table makes the post-processing functions
`filter_clusters()` and `get_cluster_regions()` cheap: no additional Monte
Carlo simulation is required to extract distinct clusters of arbitrary
rank. The user pays for the simulation once, then explores the result
interactively.

### Datasets

The package ships three real datasets, each accompanied by its tree and
(where applicable) a polygon layer for visualization:

| Dataset | Country | Domain | Regions | Tree (nodes) |
|---|---|---|---|---|
| `rj_mortality` | Brazil | Health | 89 | ICD-10 (622) |
| `chicago_crimes` | USA | Crime | 77 | Crime taxonomy (2841) |
| `london_collisions` | UK | Road safety | 33 | Conditions (81) |

Each dataset is paired with its tree object (e.g., `rj_tree`); the regions
column reports the number of spatial units.

They were chosen to span three different domains, three different
denominator conventions, and three different tree sizes, so that users
encountering the method can find an example structurally similar to their
own.

A fourth dataset, `fl_deaths`, ships in raw long form — 65 Florida
counties, 253 ICD-10 codes, ~157,000 deaths in 2016 — and is intended as a
pedagogical exercise in which the user constructs the tree and the case
matrix from scratch. It is exercised in
`inst/examples/example_florida.R`.

---

## 4. Worked example: replicating Cançado et al. (2025)

We reproduce Section 5.2 of Cançado et al. (2025), the analysis of infant
deaths in Rio de Janeiro state in 2016. The dataset, `rj_mortality`,
contains counts of infant deaths in 89 of the state's 92 municipalities,
classified by 4-character ICD-10 leaf codes. The tree, `rj_tree`, is the
ICD-10 hierarchy with 622 nodes and 410 leaves. The denominator for an
infant-mortality rate analysis is `live_births`, the number of registered
live births in 2016 (DATASUS/SINASC).

    library(treeSS)
    data(rj_mortality)
    data(rj_tree)

The scan is one function call:

    result_rj <- treespatial_scan(
      cases       = rj_mortality$cases,
      population  = rj_mortality$live_births,
      region_id   = rj_mortality$region_id,
      x           = rj_mortality$x,
      y           = rj_mortality$y,
      node_id     = rj_mortality$node_id,
      tree        = rj_tree,
      max_pop_pct = 0.50,
      nsim        = 999, seed = 2016, n_cores = 4L
    )
    print(result_rj)

The most likely cluster is the ICD-10 leaf code **P209** ("Other bacterial
pneumonias of newborn") concentrated in 18 contiguous municipalities of the
north fluminense region. Twenty-seven deaths are observed against an
expectation of approximately 3.34 under uniform risk — a relative risk of
~8× — with a Monte Carlo *p*-value below 1/(nsim+1). This matches Table 7
of Cançado et al. (2025): the same 18 municipalities (by IBGE code), 27
observed deaths, expected ≈ 3.37, log-likelihood ratio LR = 38.48 (we
obtain ≈ 38.83). The minor LR discrepancy is explained in Cançado et al.
(2025) by a SIM database update post-2016 and by the discontinuation of TCU
population estimates after 2017; the conclusions are identical.

### Distinct secondary clusters

The original paper reports more than one cluster per scan, using the
Section 5.1.1 procedure: pairs are retained in decreasing Λ as long as they
are *distinct* from those already retained. **treeSS** exposes this
directly via `filter_clusters()`:

    fc <- filter_clusters(result_rj)
    head(fc[, c("node_id", "n_regions", "cases", "expected", "llr", "pvalue")], 5)

No additional Monte Carlo simulation is performed; the secondary clusters
are read from `result_rj$secondary_clusters` and filtered for distinctness
in O(*k*²) time, where *k* is the number of clusters requested. The top
filtered candidates correspond to the rows of Table 7 of Cançado et al.
(2025).

### Visualizing the cluster

`get_cluster_regions()` returns a tidy data.frame mapping each region to
each retained cluster, ready to merge with a polygon layer. The package
does not ship `sf` polygons for Rio de Janeiro — they are easy to fetch on
demand from **geobr**:

    library(ggplot2); library(geobr); library(sf)

    rj_map <- read_municipality(code_muni = "RJ", year = 2016)
    rj_map$code6 <- as.integer(substr(rj_map$code_muni, 1, 6))
    ibge_lookup <- unique(rj_mortality[, c("region_id", "ibge_code")])

    # Most likely cluster (Figure 1)
    cr_rj_1 <- merge(
      get_cluster_regions(result_rj, n_clusters = 1L, overlap = FALSE),
      ibge_lookup, by = "region_id"
    )
    rj_p1 <- merge(rj_map, cr_rj_1, by.x = "code6", by.y = "ibge_code",
                   all.x = TRUE)

    ggplot(rj_p1) +
      geom_sf(aes(fill = factor(cluster)), color = "gray40", alpha = 0.80) +
      scale_fill_manual(values = c("1" = "#C44E52"), na.value = "gray95",
                        na.translate = FALSE, name = "Cluster") +
      theme_minimal()

The result is **Figure 1**: the 18 contiguous municipalities of the
north fluminense region, shaded red, against the rest of the state in
gray.

The same `get_cluster_regions()` call with `n_clusters = 2, overlap = TRUE`
and a `facet_wrap(~ panel)` produces a two-panel plot showing the most
likely cluster alongside the top *distinct* secondary cluster of Cançado
et al. (2025), Section 5.1.1, with no additional Monte Carlo simulation —
the secondary cluster is read from `result_rj$secondary_clusters`. This is
**Figure 2**.

A complete annotated worked example is shipped at
`system.file("examples", "example_brazil_rj.R", package = "treeSS")`.

---

## 5. Two further applications

The second and third datasets exercise the same workflow on different
domains, denominator conventions, and tree sizes.

### Reported crimes in Chicago, 2023

`chicago_crimes` records ~250,000 crime incidents across the 77 community
areas of Chicago in 2023, classified by a 2841-node tree of FBI primary
type, secondary description, and location category. The dataset ships with
two denominator columns: `population`, the total *incidents* in each area
(a compositional denominator), and `pop_residential`, the area's
*residential population* from ACS 2020 5-year estimates. Either is a
legitimate denominator and they answer different questions: incidents-per-
resident vs. share-of-citywide-incidents. We use `pop_residential` for an
incidence-rate analysis.

    data(chicago_crimes); data(chicago_tree)

    result_chi <- treespatial_scan(
      cases       = chicago_crimes$cases,
      population  = chicago_crimes$pop_residential,
      region_id   = chicago_crimes$region_id,
      x           = chicago_crimes$x,
      y           = chicago_crimes$y,
      node_id     = chicago_crimes$node_id,
      tree        = chicago_tree,
      max_pop_pct = 0.25,
      nsim        = 999, seed = 2023, n_cores = 4L
    )
    print(result_chi)

The most likely cluster is a contiguous group of community areas on
Chicago's South Side — Englewood, West Englewood, Auburn Gresham, Roseland
and immediate neighbors — driven by an over-occurrence of violent and
property crimes against residents. The **treeSS** result is consistent
with the spatial pattern that Chicago's annual public-safety reports
describe; what is novel is that **treeSS** identifies the *combination*
of crime category and geographic extent jointly, rather than running
separate spatial scans per category. **Figure 3** shows the primary
cluster; **Figure 4** shows the top two distinct clusters retained by
the procedure of Cançado et al. (2025), Section 5.1.1.

### Road traffic collisions in London, 2022

`london_collisions` records all police-reported road collisions in the 33
London boroughs in 2022 (UK STATS19), classified by a tree combining
lighting condition, road type, and junction type. The tree is small (81
nodes), the geography is small (33 boroughs), and the dataset illustrates
the binomial model: each collision is classified as *slight* or *serious-
or-fatal*, and the binomial scan asks where the *proportion* serious-or-
fatal departs from the citywide baseline.

    data(london_collisions); data(london_tree); data(london_boroughs_map)

    result_ldn <- treespatial_scan(
      cases       = london_collisions$cases,
      population  = london_collisions$population,
      region_id   = london_collisions$region_id,
      x           = london_collisions$x,
      y           = london_collisions$y,
      node_id     = london_collisions$node_id,
      tree        = london_tree,
      max_pop_pct = 0.40,
      nsim        = 999, seed = 2022, n_cores = 4L,
      model       = "poisson"
    )
    print(result_ldn)

The cluster identifies a contiguous group of outer boroughs with elevated
collisions on dark, unlit single-carriageway roads at non-signalized
junctions — a road-engineering pattern that pure-spatial scans on raw
collision counts would miss because total collisions are dominated by
central, traffic-heavy boroughs. **Figure 5** shows the primary cluster;
**Figure 6** shows the top two distinct clusters per Cançado et al.
(2025), Section 5.1.1.

For comparison, **Figure 7** shows the result of the *iterative* scan
procedure, in which the data are conditioned on each detected cluster
before a fresh Monte Carlo test for the next, with Holm–Bonferroni
correction across iterations. The two procedures answer different
questions: the Cançado et al. (2025) filter ranks all (zone, branch)
pairs found by a single scan and selects the most distinct, while the
iterative procedure performs a sequence of independent significance
tests. Both are valid; the iterative procedure tends to detect smaller
numbers of more strongly significant clusters.

    iter_ldn <- iterative_scan(
      cases       = london_collisions$cases,
      population  = london_collisions$population,
      region_id   = london_collisions$region_id,
      x           = london_collisions$x,
      y           = london_collisions$y,
      node_id     = london_collisions$node_id,
      tree        = london_tree,
      max_pop_pct = 0.40,
      max_iter    = 5L, alpha = 0.05,
      nsim        = 999, seed = 2022, n_cores = 4L
    )

A complete annotated workflow with the iterative variant is shipped at
`system.file("examples", "example_london.R", package = "treeSS")`.

---

## 6. Implementation and performance

### Architecture

R orchestrates and validates; C++ does the heavy lifting. The R code
(approximately 1,800 lines across 13 files) handles input normalization,
tree validation, zone construction, output formatting, S3 methods, and all
post-processing. Three Rcpp-exported functions in `src/treescan_core.cpp`
(~700 lines of C++) handle the inner loop:

* `mc_treespatial_cpp()` — Monte Carlo for the tree-spatial scan;
* `mc_spatial_cpp()` — Monte Carlo for the circular scan;
* `mc_treescan_cpp()` — Monte Carlo for the tree scan.

Each takes the prepared inputs (a leaf-by-region case matrix, the zone
membership matrix, the tree structure, and per-tree-node case rollups) and
returns the simulated null distribution along with the observed log-
likelihood ratio.

The likelihood-ratio evaluation for all (zone, branch) pairs is vectorized
as a single matrix multiplication: **C**_z = **F** **Z**ᵀ, where **F** is
the |*T*|×*m* matrix of branch case rollups and **Z** is the |𝒵|×*m*
binary zone-membership matrix. This avoids an explicit double loop over
zones and branches in the inner Monte Carlo iteration.

### Complexity and scaling

For preprocessing, building the family of circular zones requires O(*m*
log *m*) work per center (sorting other regions by distance) and produces
O(*m*²) candidate zones, capped by `max_pop_pct`. The tree rollup is
linear in the tree size, O(|*T*|). The Monte Carlo loop runs `nsim`
independent replicates, each of which costs O(|*T*|·|𝒵|) floating-point
operations dominated by the matrix multiply. The package holds two
matrices in memory: the leaf-by-region case matrix (O(*ℓ*·*m*)) and the
zone-membership matrix (O(|𝒵|·*m*)). For the largest dataset shipped
(`chicago_crimes`, *m* = 77, |*T*| = 2841, *ℓ* ≈ 1900), peak memory is
well under 1 GB.

### Wall-clock benchmarks

Indicative wall-clock times for the three shipped datasets at `nsim = 999`
on a typical multicore laptop:

| Dataset | Regions | Tree nodes | Serial (s) | 4 threads (s) | Speedup |
|---|---|---|---|---|---|
| `rj_mortality` | 89 | 622 | ~120 | ~30 | ~4× |
| `chicago_crimes` | 77 | 2841 | ~480 | ~110 | ~4–5× |
| `london_collisions` | 33 | 81 | ~5 | ~2 | ~2–3× |

London is small enough that thread-startup overhead dominates. The Chicago
benchmark is the dominant cost in the test suite of the package and the
only one for which we recommend `nsim = 999` with `n_cores = 4L` rather
than the serial path. The **treeSS** Rmarkdown vignette runs end-to-end in
well under ten minutes on a 2024 laptop using `n_cores = 4`, comfortably
below the R Journal's reproducibility budget.

### Validation

The package ships 48 unit tests covering input normalization, tree
validation, zone construction, Poisson and binomial likelihood ratios on
small toy inputs with known analytical answers, parallel-vs-serial
equivalence at `n_cores = 1`, and full reproducibility under fixed seeds.
Continuous integration on GitHub Actions exercises the package on Linux,
macOS and Windows with R-release and R-devel.

### Limitations

Three limitations are worth flagging.

* **Circular zones inherit Kulldorff's elongation bias.** When the true
  cluster has irregular shape, circular zones will either include
  uninvolved regions or split the cluster across several zones. The
  package does not yet support irregularly shaped zones, although this is
  a planned extension. Users with strong prior belief in elongated
  clusters can pre-filter the data spatially.
* **Degenerate trees.** A tree consisting of a single node coincides with
  no tree at all and reduces the tree-spatial scan to the circular scan.
  We recommend using `circular_scan()` directly in that case.
* **Small `nsim`.** Setting `nsim < 99` produces a *p*-value resolution
  coarser than 0.01 and is discouraged. The default of 999 gives a
  resolution of 0.001 with negligible runtime overhead at typical problem
  sizes.

---

## 7. Comparison with related packages

| Tool | Scan type | Implementation | On CRAN |
|---|---|---|---|
| **treeSS** | tree-spatial, tree, circular | R + C++ (Rcpp, OpenMP) | yes (this paper) |
| **scanstatistics** | space-time anomaly detection | R | yes |
| **rsatscan** | inherits from SaTScan | R wrapper around binary | yes |
| **SpatialEpi** | circular Kulldorff | R | yes |
| **SaTScan** | circular, elliptic, space-time | C++ standalone | no |
| **TreeScan** | tree-only, tree-time | C++ standalone | no |

**SpatialEpi** provides Kulldorff's circular spatial scan in pure R. For
purely spatial problems on small to moderate study areas it is a mature
and convenient choice; **treeSS** does not aim to replace it. However, on
the spatial-only subset of problems, `treeSS::circular_scan()` returns
equivalent point estimates and Monte Carlo *p*-values, with the benefit of
the OpenMP-parallelized C++ inner loop for large `nsim`.

**scanstatistics** addresses *space-time* anomaly detection (Allévius,
2018), a different problem in which the tree dimension is absent. The two
packages are complementary rather than competitive.

**rsatscan** is a thin R wrapper around the SaTScan binary (Kulldorff,
2024); users must install SaTScan separately and the package writes input
files to disk, calls the binary, and parses the output. SaTScan is the
reference implementation of Kulldorff's spatial and space-time scans and
offers shape options (elliptic, irregular) that **treeSS** does not.
However, SaTScan does not implement a tree-spatial scan: there is no way
to pass a tree of categories alongside the spatial inputs.

The standalone TreeScan software (Kulldorff, 2024) implements the tree-
based scan and the tree-time scan but, like SaTScan, does not perform a
joint search over space and tree. There is no R interface to TreeScan.

To our knowledge **treeSS** is therefore the first implementation of any
kind — R or otherwise — of the tree-spatial scan statistic of Cançado et
al. (2025).

---

## 8. Discussion

We have described **treeSS**, an R package that implements the tree-
spatial scan statistic of Cançado et al. (2025) alongside the classical
circular spatial scan and the tree-based scan, under a unified parallel
vector input contract. The package fills a concrete gap in the R
ecosystem: no previous CRAN package implements a scan that searches
jointly over circular spatial zones and branches of a hierarchical
classification tree. The core search is implemented in C++ via Rcpp and
optionally parallelized with OpenMP, and the package ships three real
datasets that allow users to run a complete workflow without preparing
data themselves.

We see three priorities for future work. First, **irregular zones**: the
elongation bias of circular zones is well-documented and a Duczmal-style
graph-cut zone family would extend **treeSS** to that setting with no
change to the methodology of Cançado et al. (2025) — only the zone-
construction step. Second, a **space-time-tree extension**: combining the
tree dimension of Kulldorff et al. (2003) with the space-time dimension of
Allévius (2018) is a natural three-way generalization relevant to outbreak
surveillance. Third, **native `sf` integration**: currently users supply
centroids as (*x*, *y*) columns; accepting an `sf::sf` object directly
would make the package more idiomatic for modern spatial-R workflows.

In its current form **treeSS** is, we believe, ready to support applied
work. The Rio de Janeiro reproduction in Section 4 and the Chicago and
London applications in Section 5 illustrate the breadth of problems for
which a tree-spatial perspective gives non-trivial signal that pure-
spatial or pure-tree analyses miss. We hope the package will lower the
barrier to such analyses and encourage the broader adoption of joint scan
statistics in surveillance practice.

## Acknowledgements

A.L.F.C., G.S.O. and L.H.D. acknowledge the support of CNPq and CAPES.
A.V.C.Q. thanks the maintainers of the open-source R ecosystem on which
this work is built, in particular the authors of **Rcpp**, **SpatialEpi**,
**geobr** and the **rjtools** team. The DATASUS, IBGE, City of Chicago,
and UK Department for Transport open data programs made the applied
examples in this paper possible.
