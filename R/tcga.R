#' TCGA helper function: prepare TCGA data for limma or survival analysis.
#'
#' @param data data from TCGA biolinks.
#'
#' @return list
#' @import data.table
#' @export
dprepare.tcga <- function(data) {
  sample_info <- data@colData
  sample_info$OS <- fcoalesce(sample_info$days_to_death, sample_info$days_to_last_follow_up)

  genes_info <- as.data.frame(data@rowRanges)
  count <- data@assays@data$unstranded
  rownames(count) <- rownames(genes_info)
  colnames(count) <- rownames(sample_info)

  # tumor smaple data
  tumor_sample_index <- sample_info$sample_type %ilike% "Tumor"
  count2 <- count[, tumor_sample_index]
  sample_info2 <- sample_info[, tumor_sample_index]

  list(
    all = list(
      count = count,
      sample_info = sample_info,
      genes_info = genes_info
    ),
    tumor = list(
      count = count2,
      sample_info = sample_info2,
      genes_info = genes_info
    )
  )
}
