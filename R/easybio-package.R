#' @keywords internal
"_PACKAGE"

#' @import data.table
.datatable.aware <- TRUE

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom checkmate assert_data_frame assert_string assert_subset
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

.onAttach <- function(libname, pkgname) {
  msg <- c(
    "easybio has been updated with significant breaking changes in single-cell annotation workflow.",
    "To learn the new workflow, please run:",
    '  vignette("example-single-cell-annotation", package = "easybio")'
  )
  packageStartupMessage(paste(msg, collapse = "\n"))
}
