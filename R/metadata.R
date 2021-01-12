
validateMetadataTable <- function(tbl) {
  checkColumns(c("Sample.IDs"), colnames(tbl), "metadata table")
}

#' @export
read_metadata <- function(filename) {
  tbl <- read_csv(filename)
  validateMetadataTable(tbl)
  return(tbl)
}

#' @export
reduce_metadata_by_expression <- function(expression) {
  metadata[metadata$Sample.IDs %in% colnames(expression),]
}

#'
#'
#' @export
getSamplesForFactor <- function(metadata, factor = NULL, level = "all"){
  if(is.null(factor) || level == "all") {
    return(unique(metadata$Sample.IDs))
  } else {
    return(unique(metadata$Sample.IDs[metadata[[factor]] == level]))
  }
}
