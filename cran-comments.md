## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are explained below. Neither is a defect in the package.

### 1. "Found the following (possibly) invalid URLs ... Status: 403"

The CellMarker URL is correct, and the host refuses requests that do not
look like a browser:

```
$ curl -s -o /dev/null -w '%{http_code}' https://bio-bigdata.hrbmu.edu.cn/CellMarker/
403
$ curl -s -o /dev/null -w '%{http_code}' -A 'Mozilla/5.0' https://bio-bigdata.hrbmu.edu.cn/CellMarker/
200
```

The link opens normally in a browser; only the automated checker is turned
away. The same URL is cited from `DESCRIPTION` and `man/easybio-package.Rd`,
which is where the note points.

### 2. "unable to verify current time"

Local to the machine this check ran on, which has no route to a time
service. It does not appear on machines that can reach one.

## Submission summary

This is a minor release, 1.2.2 to 1.3.0. The headline change is the
annotation database: the built-in reference moves from CellMarker 2.0 to
CellMarker 3.0 (418,139 entries, of which the five columns the package
uses are kept).

### Changes users will notice

* Every exported function and argument now uses `snake_case`. The previous
  camelCase names are kept as deprecated aliases that emit a lifecycle
  warning and are scheduled for removal in 1.4.0 — nothing was dropped, so
  existing code keeps working. `matchCellMarker2()` became `match_ref()`, because it also
  accepts a user-supplied reference and the old name implied otherwise.
* `match_ref()` now ranks candidate cell types by `uniqueN` (the number of
  distinct matching markers), with `N` as the tie-breaker, where it
  previously ranked by `N` alone: a single heavily reported marker could
  otherwise outweigh a cell type matching on dozens of markers.
* It gains `min_pct`, which drops markers whose detection rate in the
  cluster they were found for is too low, and reports that detection rate
  of each matched marker in the new `pct_with` column.
* `plot_possible_cell()` gains `value = "pct"`, which shows how much of a
  candidate's evidence is actually detected rather than how much matched.

### Package-level

* The internal database is now stored xz-compressed, 1.24 MB against
  1.62 MB. It is decompressed once at install time, so this costs nothing
  at run time.
* The startup message announcing the 1.3.0 changes was removed. It had no
  version check, so it also fired for users who installed 1.3.0 fresh.
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
