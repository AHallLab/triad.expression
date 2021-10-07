
getSamplesForFactors <- function(metadata, factorList) {
  selections <- purrr::reduce(
    purrr::imap(factorList, ~ metadata[[.y]] == .x), `&`
  )
  return(metadata$Sample.IDs[selections])
}

#' Filter a gene expression data set to specific combination of factor levels, based on the metadata describing the samples.
#'
#' @param expression A data frame of expression values in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param factorLevels A List that represents a dictionary of factor and level pairs.
#'
filterExpressionForFactors <- function(expression, metadata, factorList) {
  samples <- getSamplesForFactors(metadata, factorList)
  samplePresent <- samples %in% colnames(expression)
  samples_in_expression <- samples[samplePresent]
  if(length(samples_in_expression) < length(samples)) {
    samples_not_in_expression <- samples[!samplePresent]
    warning(
      paste0("Samples ",
             paste0(samples_not_in_expression, collapse = ", "),
             " are not present in the expression dataset."
      )
    )
  }
  #return(subset(expression, select = c("gene", samples_in_expression)))
  expression %>% select(gene, all_of(samples_in_expression))
}

#filterExpressionForFactors2 <- function(expression, metadata, factorList) {
#  expression %>%
#    select(gene, any_of(getSamplesForFactors(metadata, factorList)))
#}

