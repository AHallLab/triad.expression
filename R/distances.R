#' @importFrom fields rdist
#' @import reshape2

#' @export
mk_expression_matrix <- function(normalizedMeanExpression) {
  triadMatrix <- acast(as.data.frame(normalizedMeanExpression),
                       group_id ~ subgenome,
                       value.var = "normalised_mean")
  triadMatrix[complete.cases(triadMatrix), ]
}

#' @export
distances_from_centroids <- function(normalizedMeanExpression, lvl) {
  normalizedMeanExpression <- normalizedMeanExpression %>%
    filter(level == lvl)

  # Make a matrix for triads from input table
  triadMatrix <- mk_expression_matrix(normalizedMeanExpression)

  # Compute centroids of triads.
  centroid <- t(as.matrix(colMeans(triadMatrix,)))

  # Create a centeroid matrix where each row corresponds to an expected case:
  # central, A, B, or D dominant, or A, B, or D suppressed.
  expectationCentroids <- t(matrix(c(0.33, 0.33, 0.33, 1.0, 0.0, 0.0, 0.0,
                                     1.00, 0.00, 0.00, 0.0, 1.0, 0.0, 0.5,
                                     0.50, 0.50, 0.00, 0.5, 0.5, 0.5, 0.0), nrow = 3))
  colnames(expectationCentroids) <- c("A", "B", "D")
  rownames(expectationCentroids) <- c("Central",
                                      "A.dominant", "B.dominant", "D.dominant",
                                      "A.suppressed", "B.suppressed","D.suppressed")

  # Compute and rank distances of triads from centroid.
  distances <- rdist(triadMatrix, centroid)
  rankedDistances <- rank(distances) / length(distances)

  # Compute distances of triads from expectation centroids.
  expectationDistances <- rdist(triadMatrix, expectationCentroids)
  colnames(expectationDistances) <- rownames(expectationCentroids)
  rownames(expectationDistances) <- rownames(triadMatrix)

  mins <- apply(expectationDistances, 1, which.min)

  list(tibble(as.data.frame(cbind(triadMatrix, expectationDistances))) %>%
         mutate(group_id = rownames(triadMatrix),
                clust = mins,
                clust.description = colnames(expectationDistances)[mins]),

  normalizedMeanExpression %>%
    inner_join(tibble(distance = distances,
                      ranked_distance = rankedDistances,
                      group_id = as.integer(rownames(triadMatrix)))))

}


