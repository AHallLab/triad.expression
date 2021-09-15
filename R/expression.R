
filter_expression_for_factor_level <- function(expression, metadata, factor, level) {
  samples <- getSamplesForFactor(metadata, factor = factor, level = level)
  samples_in_expression <- samples[samples %in% colnames(expression)]
  return(expression[, c("gene", samples_in_expression)])
}

filterExpressionForVarietyAndFactorLevel <- function(expression, metadata, variety, varietyCol, level, factorCol) {
  samples <- getSamplesForVarietyAndFactor(metadata, varietyCol, variety, factorCol, level)
  samples_in_expression <- samples[samples %in% colnames(expression)]
  return(expression[, c("gene", samples_in_expression)])
}

# NEWNEWNEW
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

gene_expression_means <- function(expression) {
  mean_values <- expression %>% select(-.data$gene) %>% rowMeans
  tibble(
    mean = mean_values,
    gene = expression$gene,
    samples = ncol(expression) - 1
  )
}

#NEWNEWNEW
gene_expression_means_for_factor_levels <- function(expression, metadata, factorLevels) {
  filteredExpression <- filterExpressionForFactors(expression, metadata, factorLevels)
  mean_values <- filteredExpression %>% select(-gene) %>% rowMeans
  tibble(
    mean = mean_values,
    gene = filteredExpression$gene,
    samples = ncol(filteredExpression) - 1
  )
}

#NEWNEWNEW
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
#' NEWNEWNEW
#'
#' @param expression A data frame of expression values, in which each column is a sample, and each row is a single gene.
#' @param metadata A data frame containing metadata for each sample.
#' @param homology A data frame of metadata about each triad. Include which genes make them up on the A, B, and D genomes.
#' @param factorList A string of a metadata column name to use as a factor.
#'
#' @export
triad_expression_means_by_factors <- function(expression, metadata, homology, factorList) {
  expression_means <- gene_expression_means_by_factors(expression, metadata, factorList)
  triads_flat <- reshape_homology_triads(homology)
  return(expression_means %>% left_join(triads_flat, by = "gene") %>% na.omit)
}

#' Normalize the triad mean expression data produced by triad_expression_mean_by_factors.
#'
#' @export
normalize_triad_expression_means <- function(meanExpression, factorList = NULL) {
  if(is.null(factorList)){
    # Auto-detect the factor-columns.
    tblNames <- colnames(meanExpression)
    factorList <- tblNames[!(tblNames %in% c("mean", "gene", "samples", "subgenome", "group_id"))]
    warning(paste0("Detected factor columns. Setting factorList to ", paste0(factorList, collapse = ", ")))
  }

  meanExpression %>%
    group_by(across(all_of(c("group_id", factorList)))) %>%
    mutate(triad_sum = sum(mean), normalised_mean = mean / sum(mean)) %>%
    ungroup
}
