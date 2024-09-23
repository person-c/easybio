## To do

1. **Dependency Review**: Considering the removal of the `GEOquery` package, as it may no longer be necessary for preparing GEO series data.
2. **Plot Customization**: All `plot*` functions will be enhanced to allow greater user customization.
3. **New S3 Class**: Development of a new S3 class, similar to `dgeList`, which will offer improved customization and performance.

# New Developed features

- support argument `tissueType` and `tissueClass`.

## Version 1.1.0 Changes

### Enhancements
- **`Artist` Class**: Now includes a default data argument for data exploration. All commands and results are saved, allowing users to revisit previous analyses.
- **Bulk RNA-seq Vignettes**: Updated for better usability.
- **`get_marker()` Updates**: Added extra messages when searching for undefined cells in the database, improving clarity for users.

### Bug Fixes
- **`prepare_geo()` Fix**: Non-character ID columns in GPL data are now converted to character to prevent misinterpretation of numeric IDs. Expression data can now be read directly from supplementary files without needing to download them locally.
- **`plotORA()` Fix**: Fixed an issue where the plot would fail to display a legend when mapping a variable to the `fill` aesthetic.



## Version 1.0.1 Changes
- **NEWS.md File**: Added to document changes and updates.
- **macOS Compatibility**: Fixed an error in `uniprot_id_map()` on the latest macOS during example usage.
- **Vignette Optimization**: Vignettes have been updated and optimized for improved performance and usability.
