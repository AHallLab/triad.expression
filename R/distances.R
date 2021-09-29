#' Make a matrix from a table of normalised mean expression values
#'
#' Putting the data into a matrix form makes it easy to do distances
#' from centroids etc.
mk_expression_matrix <- function(normalizedMeanExpression) {
  triadMatrix <- acast(as.data.frame(normalizedMeanExpression),
                       group_id ~ subgenome,
                       value.var = "normalised_mean")
  triadMatrix[stats::complete.cases(triadMatrix), ]
}


#' Compute the distances of each triad from a set of centroids.
#'
#' @export
centroid_distances <- function(normalizedMeanExpression, factorList = NULL) {
  if(is.null(factorList)) {
    if(!is.null(attributes(normalizedMeanExpression)$factorList)){
      factorList <- attributes(normalizedMeanExpression)$factorList
    } else {
      # Auto-detect the factor-columns.
      tblNames <- colnames(normalizedMeanExpression)
      factorList <- tblNames[!(tblNames %in% c("mean", "gene", "samples", "subgenome", "group_id", "triad_sum", "normalised_mean"))]
      warning(paste0("Detected factor columns. Setting factorList to ", paste0(factorList, collapse = ", ")))
    }
  }

  # Central, A, B, or D dominant, or A, B, or D suppressed.
  expectationCentroids <- t(matrix(c(0.33, 0.33, 0.33, 1.0, 0.0, 0.0, 0.0,
                                     1.00, 0.00, 0.00, 0.0, 1.0, 0.0, 0.5,
                                     0.50, 0.50, 0.00, 0.5, 0.5, 0.5, 0.0), nrow = 3))
  colnames(expectationCentroids) <- c("A", "B", "D")
  rownames(expectationCentroids) <- c("Central",
                                      "A.dominant", "B.dominant", "D.dominant",
                                      "A.suppressed", "B.suppressed","D.suppressed")

  result <- normalizedMeanExpression %>%
    group_by(across(all_of(factorList))) %>%
    group_modify(~ {
      triadMatrix <- mk_expression_matrix(.x)
      expectationDistances <- rdist(triadMatrix, expectationCentroids)
      colnames(expectationDistances) <- rownames(expectationCentroids)
      rownames(expectationDistances) <- rownames(triadMatrix)
      mins <- apply(expectationDistances, 1, which.min)
      minDesc <- colnames(expectationDistances)[mins]
      tibble(cbind(as.data.frame(expectationDistances), triadMatrix)) %>%
        mutate(group_id = as.numeric(rownames(triadMatrix)),
               clust = mins,
               clust.description = minDesc) %>%
        inner_join(.x, by = "group_id")
    }) %>% ungroup

  attributes(result)$factorList <- factorList
  return(result)
}

#' @export
join_distances_and_annotation <- function(distances, annotation, by = c("group_id" = "group_id", "High.level.variety" = "variety", "subgenome" = "subgenome")){
  return(distances %>%
           left_join(annotation, by = by) %>%
           rename(cs.gene = .data$gene.x, variety.gene = .data$gene.y))
}

limit_to_regions <- function() {

}



squash_expectation_distances <- function(data) {
  data %>%
    select(group_id, clust, clust.description) %>%
    distinct
}
