## R CMD check results

0 errors | 0 warnings | 0 note


# Version 1.1.0(minor changes)

#### Enhancements

- **`Artist` Class**: Now includes a default data argument for data exploration. All commands and results are saved, allowing users to revisit previous analyses.
- **Bulk RNA-seq Vignettes**: Updated for better usability.
- **`get_marker()` Updates**: Added extra messages when searching for undefined cells in the database, improving clarity for users.
- Add citation.

#### Bug Fixes

- **`prepare_geo()` Fix**: Non-character ID columns in GPL data are now converted to character to prevent misinterpretation of numeric IDs. Expression data can now be read directly from supplementary files without needing to download them locally.
- **`plotORA()` Fix**: Fixed an issue where the plot would fail to display a legend when mapping a variable to the `fill` aesthetic.

## easybio 1.0.1(patch)

* Added a `NEWS.md` file to track changes to the package.
* Avoid error of function `uniprot_id_map()` in latest mac in examples.
* update vignettes.

## 20240903 

- add return tag in Artist.R and import R6 package
- update description and vignettes
- improve some functions

## 20240827

- add link in description
- return value in  Artist.R `plotGSE` `plotORA`, `plotRank`
- remove writing in default patt in and return original user's settings in inst/example-single-cell.R; R/misc.R; R/sc.R
- update doc and vignettes
