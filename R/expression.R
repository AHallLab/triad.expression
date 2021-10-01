
# filter_expression_for_factor_level <- function(expression, metadata, factor, level) {
#   samples <- getSamplesForFactor(metadata, factor = factor, level = level)
#   samples_in_expression <- samples[samples %in% colnames(expression)]
#   return(expression[, c("gene", samples_in_expression)])
# }
#
# filterExpressionForVarietyAndFactorLevel <- function(expression, metadata, variety, varietyCol, level, factorCol) {
#   samples <- getSamplesForVarietyAndFactor(metadata, varietyCol, variety, factorCol, level)
#   samples_in_expression <- samples[samples %in% colnames(expression)]
#   return(expression[, c("gene", samples_in_expression)])
# }

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

#' Filter a gene expression data set to specific combination of factor levels, based on the metadata describing the samples.
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
  return(subset(expression, select = c("gene", samples_in_expression)))
}

#' Compute mean expression for each gene, for a single combination of factor levels.
#'
#' @param expression A data frame of expression values in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param factorLevels A List that represents a dictionary of factor and level.
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

#' Compute mean expression for each gene for every combination of factor levels.
#'
#' @param expression A data frame of expression values in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param homology A data frame of metadata about each triad. Include which genes make them up on the A, B, and D genomes.
#' @param factorList A string of a metadata column name to use as a factor.
#'
#' @export
gene_expression_means_by_factors <- function(expression, metadata, factorList) {
  factorLevels <- lapply(factorList, function(factor) { unique(metadata[[factor]]) })
  names(factorLevels) <- factorList
  factorCombos <- expand.grid(factorLevels)
  do.call(rbind, apply(factorCombos, 1, function(x){
    tbl <- gene_expression_means_for_factor_levels(expression, metadata, as.list(x))
    for(f in factorList){
      tbl[[f]] <- x[f]
    }
    return(tbl)
  }))
}

#' Compute mean expression for each triad for each combination of factor levels.
#'
#' @param expression A data frame of expression values, in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param homology A data frame of metadata about each triad. Include which genes make them up on the A, B, and D genomes.
#' @param factorList A string of a metadata column name to use as a factor.
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
