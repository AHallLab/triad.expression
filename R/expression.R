
#' Normalize the triad mean expression data produced by triad_expression_mean_by_factors.
#'
#' @export
normalize_triad_expression_means <- function(meanExpression, factorList = NULL) {
  if(is.null(factorList)) {
    if(!is.null(attributes(meanExpression)$factorList)){
      factorList <- attributes(meanExpression)$factorList
    } else {
      # Auto-detect the factor-columns.
      tblNames <- colnames(meanExpression)
      factorList <- tblNames[!(tblNames %in% c("mean", "gene", "samples", "subgenome", "group_id"))]
      warning(paste0("Detected factor columns. Setting factorList to ", paste0(factorList, collapse = ", ")))
    }
  }

  result <- meanExpression %>%
    group_by(across(all_of(c("group_id", factorList)))) %>%
    mutate(triad_sum = sum(mean), normalised_mean = mean / sum(mean)) %>%
    ungroup
  attributes(result)$factorList <- factorList
  return(result)
}
