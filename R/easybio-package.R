#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @import ggplot2
.datatable.aware <- TRUE # nolint: object_name_linter.

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
#
# All imports are declared here rather than next to the functions that use
# them: roxygen's @import is package-wide wherever it is written, so repeating
# it in each file that called ggplot2() only looked like it belonged to that
# file. Packages used in a handful of places (lifecycle, httr2, xml2) are
# called with :: instead and need no import at all, and the base packages are
# named function by function -- importing stats, utils, graphics and grDevices
# whole pulled in four namespaces for a dozen functions.
##
## R6 is imported even though Artist is its only use: a declared dependency
## that nothing is imported from is a NOTE in R CMD check, and the class is
## built once, so there is nothing to gain from keeping R6:: on the call.
## usethis namespace: start
#' @importFrom checkmate assert_character assert_data_frame assert_names assert_number assert_string assert_subset
#' @importFrom grDevices rainbow
#' @importFrom graphics boxplot lines par plot title
#' @importFrom R6 R6Class
#' @importFrom stats density formula model.matrix na.omit reorder setNames t.test wilcox.test
#' @importFrom utils adist combn head
## usethis namespace: end
NULL

#' Example marker data from Seurat::FindAllMarkers()
#'
#' The data were obtained by the seurat PBMC workflow.
#' exact script for this data is available as system.file("example-single-cell.R", package="easybio")
#' @docType data
#' @name pbmc.markers
NULL

#' Example DEGs data from Limma-Voom workflow for TCGA-CHOL project
#'
#' The data were obtained by the limma-voom workflow
#' @docType data
#' @name CHOL_DEGs
NULL
