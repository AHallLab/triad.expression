#' Compute mean expression for each gene for every combination of factor levels.
#'
#' @param expression A data frame of expression values in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param homology A data frame of metadata about each triad. Include which genes make them up on the A, B, and D genomes.
#' @param factorList A list of column names that exist in `metadata` (as strings) to use as a factors.
#'
#' @export
gene_expression_means_by_factors <- function(expression, metadata, factorList) {
  factorLevels <- lapply(factorList, function(factor) { unique(metadata[[factor]]) })
  names(factorLevels) <- factorList
  factorCombos <- expand.grid(factorLevels)
  # do.call(rbind, apply(factorCombos, 1, function(x){
  #   tbl <- gene_expression_means_for_factor_levels(expression, metadata, as.list(x))
  #   for(f in factorList){
  #     tbl[[f]] <- x[f]
  #   }
  #   return(tbl)
  # }))
  bind_rows(apply(factorCombos, 1, function(x){
    tbl <- gene_expression_means_for_factor_levels(expression, metadata, as.list(x))
    for(f in factorList){
      tbl[[f]] <- x[f]
    }
    return(tbl)
  }))
}
