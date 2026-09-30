# To do

1. **Dependency Review**: Considering the removal of the `GEOquery` package, as it may no longer be necessary for preparing GEO series data.
2. **Plot Customization**: All `plot*` functions will be enhanced to allow greater user customization.
3. **New S3 Class**: Development of a new S3 class, similar to `dgeList`, which will offer improved customization and performance.

# Version 1.3.0 Changes

- all exported functions follow the snake_case naming style; camelCase names are kept as deprecated aliases and will be removed in the next version (lintr config now enforces snake_case).
- `list2dt()` and `list2graph()` are the exception and keep their original names: the `2`-for-"to" idiom matches base R (`list2DF()`, `list2env()`) and both names satisfy the snake_case linter, so renaming them would have broken a released name for no benefit. The briefly considered `list_to_dt()`/`list_to_graph()` never reached a release.
- function parameters renamed to snake_case as well (e.g., `avg_log2fc_threshold`, `top_cell_n`, `tissue_class`, `tissue_type`, `min_count`, `min_unique_n`, `ignore_case`, `group_column`, `sample_info`, `feature_info`, `data_text`, `gsea_param`, `ticks_size`, `fgsea_res`).
- `prepare_tcga()` returns `expr_count` and `expr_fpkm` instead of `exprCount` and `exprFpkm`, so every field of its two tables is snake_case. The old names still work and warn; they will stop working in 1.4.0.
- upgrade the built-in annotation database from CellMarker 2.0 to CellMarker 3.0 (418,139 entries; only the columns used by the package are kept, see `data-raw/cellmarker.R`).
- `tuneParameters()` now stores the annotation in the meta.data column "CellMarker3.0" instead of "CellMarker2.0".
- `matchCellMarker2()` renamed to `match_ref()` because it also supports custom reference datasets. The old name is deprecated and will be removed in the next version.
- the built-in annotation dataset is renamed from `cellMarker2` to `cellMarker3` (internal; reflects the CellMarker 3.0 source).
- `match_ref()` ranks candidate cell types by `uniqueN` (breadth of matched markers) with `N` as tie-breaker instead of by `N` alone; a single heavily reported marker could otherwise outweigh a cell type matching on dozens of markers.
- `match_ref()` reports the detection rate (`pct.1`) of each matching marker in the new `pct_with` column, aligned with `ordered_symbol` and `NA` when the input has no `pct.1`. It is reported for auditing and does not affect the ranking; see the `min_pct` argument below if you want detection rate to influence the annotation.
- `match_ref()` gains a `min_pct` argument: markers whose detection rate in the cluster they were found for (`pct.1`) is below it are dropped before matching, so an annotation cannot rest on genes that are barely detected. It defaults to `NULL` (no filtering, i.e. the previous behaviour) and is ignored, with a message, if the input has no `pct.1` column. Note that `n` is applied after this gate, so the markers used are not a subset of the unfiltered ones. `Seurat::FindAllMarkers(min.pct = )` cannot replace it: that argument gates which genes are tested in either population, so a gene detected in a handful of cells can still end up as a positive marker.
- `plot_possible_cell()` gains a `value` argument: `"N"` (the default, unchanged), `"uniqueN"` (the measure `match_ref()` ranks candidates by) or `"pct"`, which draws the share of matching markers that are actually detected at `min_pct` as tiles labelled with `uniqueN` — the evidence check for an annotation. `value = "pct"` errors when the input carries no detection rate. `min_unique_n` now includes candidates matching exactly that many markers; it previously required more than that.
- the deprecated aliases report "deprecated as of 1.3.0" when called and are scheduled for removal in 1.4.0.

# Version 1.2.3 Changes

- dont run `Seurat::DotPlot` in vignette to avoid upstream upgrade's affects

# Version 1.2.2 Changes

- fix vignette name in package startup message

# Version 1.2.1 Changes

- make each function do a simple task instead of integrating all to a function(`check_marker`, `plotSeuratDot`).
- support using custom dataset to annotate and allow user to set the threshold(`matchCellMarker2`).
- use more convenient input in `finsert`.
- guess user's typo input (`get_marker`).
- update docs and hints.
- set the minimal `data.table` version to 1.15.0
- support generating plots for each cell type in `plotSeuratDot`.


# Version 1.1.1 Changes

- support argument `tissueType` and `tissueClass` in single cell related function.
- fix some typo errors in example-single-cell.R
- update some miscellaneous functions.

# Version 1.1.0 Changes

### Enhancements

- **`Artist` Class**: Now includes a default data argument for data exploration. All commands and results are saved, allowing users to revisit previous analyses.
- **Bulk RNA-seq Vignettes**: Updated for better usability.
- **`get_marker()` Updates**: Added extra messages when searching for undefined cells in the database, improving clarity for users.

### Bug Fixes

- **`prepare_geo()` Fix**: Non-character ID columns in GPL data are now converted to character to prevent misinterpretation of numeric IDs. Expression data can now be read directly from supplementary files without needing to download them locally.
- **`plotORA()` Fix**: Fixed an issue where the plot would fail to display a legend when mapping a variable to the `fill` aesthetic.



# Version 1.0.1 Changes

- **NEWS.md File**: Added to document changes and updates.
- **macOS Compatibility**: Fixed an error in `uniprot_id_map()` on the latest macOS during example usage.
- **Vignette Optimization**: Vignettes have been updated and optimized for improved performance and usability.
