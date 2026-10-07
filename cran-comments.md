## R CMD check results

0 errors | 0 warnings | 0 notes

### The note in the pretest

The pretest reported one note, a possibly invalid URL in `README.md`: the
codecov badge's link target, `https://codecov.io/gh/person-c/easybio`, now
answers 301 and points at `https://app.codecov.io/gh/person-c/easybio`.
`README.md` links to the new URL. The badge image is left alone because it
is not affected: `codecov.io` still serves `/graph/badge.svg` directly,
with no redirect to follow.

## Submission summary

This is a minor release, 1.2.3 to 1.3.0. The headline change is the
annotation database: the built-in reference moves from CellMarker 2.0 to
CellMarker 3.0 (418,139 entries, of which the five columns the package
uses are kept).

### Changes users will notice

* `available_tissue_type()` takes the tissue class to look in, and
  `match_ref()` warns when its two tissue filters select no reference entry
  at all. `tissue_class` and `tissue_type` read like a hierarchy but are two
  labels recorded per entry and combined with AND, so a class and a type
  that never occur together used to return no candidate without saying so.
  The single-cell vignette and example script now bound the search for the
  PBMC example (`tissue_class = c("Blood", "Bone marrow")`) and explain why.
* Every exported function and argument now uses `snake_case`. The previous
  camelCase names are kept as deprecated aliases that emit a lifecycle
  warning and are scheduled for removal in 1.4.0 — nothing was dropped, so
  existing code keeps working. `matchCellMarker2()` became `match_ref()`, because it also
  accepts a user-supplied reference and the old name implied otherwise.
  For the same reason the two tables returned by `prepare_tcga()` report
  `expr_count` and `expr_fpkm`: a list field cannot be aliased the way a
  function can, so the old `exprCount`/`exprFpkm` names still resolve but
  warn, and they too go in 1.4.0.
* `match_ref()` now ranks candidate cell types by `uniqueN` (the number of
  distinct matching markers), with `N` as the tie-breaker, where it
  previously ranked by `N` alone: a single heavily reported marker could
  otherwise outweigh a cell type matching on dozens of markers.
* It gains `min_pct`, which drops markers whose detection rate in the
  cluster they were found for is too low, and reports that detection rate
  of each matched marker in the new `pct_with` column.
* `plot_possible_cell()` gains `value = "pct"`, which shows how much of a
  candidate's evidence is actually detected rather than how much matched.
* The startup message is rewritten. The one in 1.2.3 announced "significant
  breaking changes in single-cell annotation workflow" without naming any of
  them, and was printed at every attach, non-interactive ones included. The new
  text names the two changes a user will actually meet — the database moving to
  CellMarker 3.0, which happened underneath the same function names and the
  same arguments, and the snake_case renames — and is printed only in an
  interactive session, so `R CMD check`, CI and `Rscript` stay quiet. It is a
  `packageStartupMessage()`, so `suppressPackageStartupMessages()` silences it.

### Package-level

* The internal database is now stored xz-compressed, 1.24 MB against
  1.62 MB. It is decompressed once at install time, so this costs nothing
  at run time.
* No new hard dependency; `BiocParallel` joins `Suggests` because a
  vignette now passes it to `fgsea()` explicitly (see below).

## Test environments

* local: Windows 11, R 4.5.0 (ucrt)
* GitHub Actions, on every push and pull request: ubuntu-latest with R
  release and R devel, macos-latest, windows-latest
  (`.github/workflows/R-CMD-check.yaml`), and the same workflow runs
  `R CMD check --as-cran`, `lintr::lint_package()` and the coverage report.

## Notes

The `fgsea()` call in `vignettes/example_limma.Rmd` passes
`BPPARAM = BiocParallel::SerialParam()`. Left to itself, fgsea hands the
multilevel step to a BiocParallel worker, which intermittently fails to
find fgsea's own compiled function and aborts the vignette build; the
example is small enough that running it on one thread costs nothing. This
is also why `BiocParallel` is named in `Suggests`.
