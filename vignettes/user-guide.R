## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=FALSE-----------------------------------------------------------
library(triad.expression)
library(dplyr) # dplyr is going to be useful too.

## -----------------------------------------------------------------------------
?expression_data

## -----------------------------------------------------------------------------
data(expression_data)
data(expression_metadata)
data(triad_homology)

## -----------------------------------------------------------------------------
meanExpressionByVariety <- triad_expression_means_by_factors(expression_data, expression_metadata, triad_homology, "High.level.variety")

## -----------------------------------------------------------------------------
meanExpressionByVarietyAndTissue <- triad_expression_means_by_factors(expression_data, expression_metadata, triad_homology, c("High.level.variety", "High.level.tissue"))

## -----------------------------------------------------------------------------
normalizedMeans <- normalize_triad_expression_means(meanExpressionByVarietyAndTissue)

## -----------------------------------------------------------------------------
cDistances <- triad_balance(normalizedMeans)

## -----------------------------------------------------------------------------
data(expression_data)
data(expression_metadata)
data(triad_homology)

meanExpressionByVariety <- triad_expression_means_by_factors(expression_data, expression_metadata, triad_homology, "High.level.variety")

normalizedMeans <- normalize_triad_expression_means(meanExpressionByVariety)

cDistances <- triad_balance(normalizedMeans)

## -----------------------------------------------------------------------------
data(gene_locations)
distWithLocation <- join_distances_and_annotation(cDistances, gene_locations)


## ---- fig.width=7-------------------------------------------------------------
distWithLocation %>%
  filter(chr == "chr2A")  %>%
  group_by(group_id) %>%
  filter(n_distinct(clust.description) > 1) %>%
  loom_plot("start", "High.level.variety", "group_id", "clust.description", yLab = "Variety")

## ---- fig.width=7-------------------------------------------------------------
distWithLocation %>%
  filter(chr == "chr2A")  %>%
  group_by(group_id) %>%
  filter(
    "suppressed" %in% substr(unique(clust.description), 3, 100L),
    "dominant" %in% substr(unique(clust.description), 3, 100L)
  ) %>%
  loom_plot("start", "High.level.variety", "group_id", "clust.description", yLab = "Variety")


