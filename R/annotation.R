
#'
#' @export
load_gene_annotations_for_varieties <- function(fileList = varietyInfo$triad_location_paths) {
  outList <- lapply(names(fileList), function(variety) {
    file <- fileList[[variety]]
    tbl <- as_tibble(read.delim(file))
    dtbl <- distinct(tbl)
    if(nrow(tbl) > nrow(dtbl)) {
      warning(paste0("Duplicate entries were detected in ", variety, " these will be removed"))
    }
    dtbl$variety <- variety
    dtbl$subgenome <- substr(dtbl$Chr, 5, 5)
    checkColumns(c("Chr", "Start", "End", "Strand", "Gene", "TriadGrp", "variety"), colnames(dtbl), "triad location")
    return(dtbl)
  })
  names(outList) <- names(fileList)
  return(outList)
}
