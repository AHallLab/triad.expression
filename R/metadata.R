
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

getSamplesForVarietyAndFactor <- function(metadata, varietyCol, variety, factorCol = NULL, level = "all") {
  varietyMetadata <- metadata[metadata[[varietyCol]] == variety, ]
  if(is.null(factor) || level == "all") {
    return(unique(varietyMetadata$Sample.IDs))
  } else {
    return(unique(varietyMetadata$Sample.IDs[varietyMetadata[[factorCol]] == level]))
  }
}

getSamplesForFactors <- function(metadata, factorList) {
  conditions <- lapply(names(factorList), function(factor) {
    level <- factorList[[factor]]
    return(metadata[[factor]] == level)
  })
  if(length(conditions) > 1) {
    conditions <- do.call("&", conditions)
  } else {
    conditions <- unlist(conditions)
  }
  return(metadata$Sample.IDs[conditions])
}


