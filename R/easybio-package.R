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

# The release note is shown once per version rather than on every attach: the
# changes below are news only the first time, and repeating them at every
# library() call is noise. An earlier version of this notice had no version
# check and was removed for exactly that reason. The version last announced to
# this user is kept under tools::R_user_dir(), the directory CRAN policy allows
# a package to write to.
.easybio_notice_file <- function() {
  file.path(tools::R_user_dir("easybio", "config"), "announced")
}

# R CMD check reads a writeLines() in the body of a startup hook as an attempt
# to print to the console and NOTEs about it, even when it writes to a file, so
# the write is kept out of .onAttach(). Calls are only inspected inside the
# hook itself, not inside what it calls.
.easybio_record_version <- function(version) {
  stamp <- .easybio_notice_file()
  dir.create(dirname(stamp), recursive = TRUE, showWarnings = FALSE)
  writeLines(version, stamp)
}

.easybio_notice <- c(
  "* The built-in annotation database is now CellMarker 3.0",
  "* Exported functions and arguments are snake_case; the camelCase names",
  "  still work and warn, and are removed in 1.4.0",
  '* The workflow: vignette("example-sc-seq-workflow", package = "easybio")'
)

# Kept apart from .onAttach() so the version logic can be exercised directly;
# there is no console to test for in here.
.easybio_announce <- function(pkgname) {
  current <- as.character(utils::packageVersion(pkgname))
  stamp <- .easybio_notice_file()
  announced <- tryCatch(
    if (file.exists(stamp)) trimws(readLines(stamp, warn = FALSE)[1L]) else NA_character_,
    error = \(e) NA_character_
  )

  # only an upgrade is news; a downgrade should not bring the notice back
  if (!is.na(announced) && utils::compareVersion(current, announced) <= 0L) {
    return(invisible())
  }

  packageStartupMessage(
    paste(c(paste("easybio", current), .easybio_notice), collapse = "\n")
  )

  # recorded whether or not the user suppresses startup messages, so that the
  # notice is not shown to them again on every session
  tryCatch(.easybio_record_version(current), error = \(e) NULL)
}

.onAttach <- function(libname, pkgname) {
  # The notice is for a person at the console. R CMD check, CI and Rscript
  # attach the package too, and whichever of them gets there first would
  # otherwise spend the notice on a log nobody reads: a check run of this very
  # hook is what kept the 1.3.0 notice from being seen. Staying quiet in a
  # non-interactive session also leaves the version unrecorded, so the user
  # still gets it the first time they attach it themselves.
  if (!interactive()) {
    return(invisible())
  }

  .easybio_announce(pkgname)
}
