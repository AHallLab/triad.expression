
#'
#'
#' @export
overlapping_genes_and_regions <- function(triadData, regions) {
  factorList <- attributes(triadData)$factorList
  sharedFactors <- factorList[na.omit(match(colnames(regions), factorList))]
  message(paste0("Shared factors: ", sharedFactors))

  triadData <- triadData[!(is.na(triadData$chr) | is.na(triadData$start) | is.na(triadData$end)), ]

  overlaps <- as.data.frame(findOverlaps(makeGRangesFromDataFrame(triadData, keep.extra.columns = T), makeGRangesFromDataFrame(regions, keep.extra.columns = T)))

  t <- triadData[overlaps$queryHits, ]
  r <- regions[overlaps$subjectHits, ]
  colnames(r) <- paste0("region.", colnames(r), sep = "")
  overlapping <- cbind(t, r)

  factorFilter <- purrr::reduce(lapply(sharedFactors, function(factor) {
    overlapping[[factor]] == overlapping[[paste0("region.", factor, sep = "")]]
  }), "&")

  overlapping <- overlapping[factorFilter, ]
  return(overlapping)
}


