# R Journal paper for `treeSS`

This directory contains the source material for the R Journal paper
describing the **treeSS** package. It is excluded from the package
build via `.Rbuildignore` and lives only in the GitHub repository.

## Files

| File | Purpose |
|---|---|
| `treeSS.Rmd` | Paper source in `rjtools::rjournal_*_article` format |
| `treeSS.bib` | BibTeX references cited in the paper |
| `motivation.md` | Cover/motivation letter for submission |
| `README.md` | This file |

## Building locally

The paper is built with the `rjtools` package (CRAN), which produces
both PDF and HTML in the format the R Journal expects.

```r
install.packages("rjtools")
library(rjtools)

# If pandoc errors, pin a known-good version first:
# pandoc::pandoc_activate(version = "3.1.6")

rmarkdown::render("treeSS.Rmd", output_format = "all")
```

This produces `treeSS.pdf` and `treeSS.html`.

## Pre-submission checks

```r
rjtools::initial_check_article(".")
```

Individual checks of interest:

- `check_proposed_pkg()` -- treeSS must be on CRAN (see "Submission
  checklist" below).
- `check_packages_available()` -- all `\CRANpkg{}` references must
  resolve to CRAN/Bioc.
- `check_title()` -- title in title case.
- `check_section()` -- section headings in sentence case.
- `check_spelling()`.

## Submission checklist

These are the requirements at upload time, drawn from
<https://journal.r-project.org/submissions.html> and
<https://journal.r-project.org/R_package_guidelines.html>:

- [ ] `treeSS` is on CRAN.
- [ ] Paper is at most 20 pages.
- [ ] Reproducibility: `treeSS.Rmd` knits in under 10 minutes
      (cache enabled; main scan uses `n_cores = 4`).
- [ ] Cover letter (`motivation.md`) describes why the paper is
      suitable for the R Journal.
- [ ] `_Rpackages.txt` lists every package needed to reproduce.
- [ ] All references in `treeSS.bib` are actually cited; no orphans.
- [ ] Repository tagged with the version corresponding to the
      submitted paper.
- [ ] pkgdown site live at <https://allanvc.github.io/treeSS>.
- [ ] Zenodo DOI minted from a tagged GitHub release.

## Submission

Upload zip via <https://forms.gle/Eqkf6cFJM3mjuxZUA>. Zip must contain:

- `treeSS.pdf`
- `treeSS.Rmd`, `treeSS.bib`, any `.tex`/`.sty`/figure files needed
  to rebuild
- `_Rpackages.txt`
- `motivation.md` (cover letter)
- supplementary files if any

Zip size limit: 10 MB.
