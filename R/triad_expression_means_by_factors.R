#' Compute mean expression for each triad for each combination of factor levels.
#'
#' @param expression A data frame of expression values, in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param homology A data frame of metadata about each triad. Include which genes make them up on the A, B, and D genomes.
#' @param factorList A list of column names that exist in `metadata` (as strings) to use as a factors.
#'
#' @examples
#' data(expression_data)
#' data(expression_metadata)
#' data(triad_homology)
#'
#' meanExpressionByVarietyAndTissue <- triad_expression_means_by_factors(expression_data, expression_metadata, triad_homology, c("High.level.variety", "High.level.tissue"))
#'
#' @export
triad_expression_means_by_factors <- function(expression, metadata, homology, factorList) {
  expression_means <- gene_expression_means_by_factors(expression, metadata, factorList)
  triads_flat <- reshape_homology_triads(homology)
  result <- expression_means %>% left_join(triads_flat, by = "gene") %>% na.omit
  attributes(result)$factorList <- factorList
  return(result)
}
