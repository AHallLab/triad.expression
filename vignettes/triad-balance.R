## ---- warning=FALSE-----------------------------------------------------------
library(triad.expression)

data(expression_data)
data(expression_metadata)
data(triad_homology)

## -----------------------------------------------------------------------------
meanExpressionByVariety2 <- triad_expression_means_by_factors(
  tibble::as_tibble(expression_data),
  expression_metadata,
  tibble::as_tibble(triad_homology),
  "High.level.variety"
)

## -----------------------------------------------------------------------------
meanExpressionByVarietyAndTissue <- triad_expression_means_by_factors(
  tibble::as_tibble(expression_data),
  expression_metadata,
  triad_homology,
  c("High.level.variety", "High.level.tissue")
)

## -----------------------------------------------------------------------------
normalizedMeans <- normalize_triad_expression_means(
  meanExpressionByVarietyAndTissue
)

## -----------------------------------------------------------------------------
eBalances <- triad_balance(normalizedMeans)

