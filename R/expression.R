
#' @export
read_expression <- function(file) {
  return(read_csv(file))
  # TODO: Add validation of table.
}

#' @export
filter_expression_for_factor_level <- function(expression, metadata, factor, level) {
  samples <- getSamplesForFactor(metadata, factor = factor, level = level)
  samples_in_expression <- samples[samples %in% colnames(expression)]
  return(expression[, c("gene", samples_in_expression)])
}

#' @export
gene_expression_means <- function(expression) {
  mean_values <- expression %>% select(-gene) %>% rowMeans
  tibble(
    mean = mean_values,
    gene = expression$gene,
    samples = ncol(expression) - 1
  )
}

#' @export
gene_expression_means_by_factor_levels <- function(expression, metadata, factor) {
  factorLevels <- c("all", unique(metadata[[factor]]))

  exprMeans <- do.call(rbind, lapply(factorLevels, function(level){
    expressionForLevel <- filter_expression_for_factor_level(expression, metadata, factor, level)
    exprMean <- gene_expression_means(expressionForLevel)
    exprMean$level <- level
    return(exprMean)
  }))

  exprMeans <- rbind(
    exprMeans,
    exprMeans %>%
      filter(level != "all") %>%
      group_by(gene) %>%
      summarise(
        mean = mean(mean),
        samples = n_distinct(level),
        level = "all_levels")
  )
}

#' @export
triad_expression_mean_by_factor_levels <- function(metadata, homology, expression, factor) {
  expression_means <- gene_expression_means_by_factor_levels(expression, metadata, factor)
  triads_flat <- reshape_homology_triads(homology)
  return(expression_means %>% left_join(triads_flat, by = "gene"))
}

#' @export
normalize_mean_triad_expression <- function(meanExpression) {
  meanExpression %>%
    group_by(group_id, level) %>%
    mutate(triad_sum = sum(mean), normalised_mean = mean / sum(mean))
}
