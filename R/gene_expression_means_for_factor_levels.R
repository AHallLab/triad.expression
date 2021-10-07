#' Compute mean expression for each gene, for a single combination of factor levels.
#'
#' @param expression A data frame of expression values in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param factorLevels A List that represents a dictionary of factor and level pairs.
#'
#'
gene_expression_means_for_factor_levels <- function(expression, metadata, factorLevels) {
  filteredExpression <- filterExpressionForFactors(expression, metadata, factorLevels)
  mean_values <- filteredExpression %>% select(-gene) %>% rowMeans
  tibble(
    mean = mean_values,
    gene = filteredExpression$gene,
    samples = ncol(filteredExpression) - 1
  )
}
