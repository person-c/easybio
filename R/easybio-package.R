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

# --- Startup notice ---

# The notice is shown at every attach, not once: a message that appears a
# single time is easy to scroll past, and the one change here that has no other
# way of announcing itself is the database, which changed underneath the same
# function names and the same arguments. The camelCase names do warn at the
# point of use, but only for someone who is still calling them.
#
# The console test is what keeps R CMD check, CI and Rscript quiet -- they
# attach the package as well, and there is nobody there to read a notice.
#
# The text names the release it describes, so 1.4.0 has to rewrite or drop it
# along with the aliases it mentions.
.easybio_notice <- c(
  "* The built-in annotation database is now CellMarker 3.0",
  "* Exported functions and arguments are snake_case; the camelCase names",
  "  still work and warn, and are removed in 1.4.0",
  '* The workflow: vignette("example-sc-seq-workflow", package = "easybio")'
)

.onAttach <- function(libname, pkgname) {
  if (!interactive()) {
    return(invisible())
  }

  header <- paste("easybio", as.character(utils::packageVersion(pkgname)))
  packageStartupMessage(paste(c(header, .easybio_notice), collapse = "\n"))
}
