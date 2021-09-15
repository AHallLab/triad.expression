# Analysis of the 170521_panTrans_commonRef data.
# Ben J. Ward - 2021.

library(dplyr)
library(triad.expression)

data(triad_homology)
data(expression_data)
data(expression_metadata)
data(gene_locations)

pvdist <- per_variety_expectation_distances(expression_metadata,
                                            expression_data,
                                            triad_homology,
                                            varietyCol = "High.level.variety",
                                            fact = "High.level.tissue",
                                            level = "all")



pvdist <- per_variety_expectation_distances(metadata, expression, homology, varietyCol = "High.level.variety", fact = "High.level.tissue", level = "all")
pvjoin <- pv_join_dists_and_annotation(pvdist, varietyGeneAnnotations)

chr1A.plot <- par_coord_plot(pvjoin, "chr1A", varieties, showSuppression = T)
chr2B.plot <- par_coord_plot(pvjoin, "chr2B", varieties, showSuppression = T)

par_coord_plot(pvjoin, "chr1B", varieties)
par_coord_plot(pvjoin, "chr1D", varieties)

# par_coord_plot(pvjoin, "chr2A", varieties)
# Alien-stuff.
# Overlay the categories.


chromosomes <- apply(expand.grid("chr", 1:7, c("A", "B", "D")), 1, function(x) { paste0(x, collapse="") })

combinations <- expand.grid(names(pvjoin), names(pvjoin), chromosomes) %>% filter(Var1 != Var2)
allPlots <- apply(combinations, 1, function(x) {pv_xy_plot(pvjoin, x[[1]], x[[2]], x[[3]])})





squashedAnalyses <- lapply(varieties, function(v) {
  squashedAnalysis <- triad.expression::squashExpressionAnalysis(perVarietyExpectationDistances, v)
  lapply(chromosomes, function(chr) {
    squashedLocations <- triad.expression::squashGeneLocations(geneLocations, chr, v)
    joinedSquashed <- squashedLocations %>%
      left_join(squashedAnalysis, by = c("TriadGrp" = "group_id"))
  })
})


triads_string <- squashedTableToString(joinedSquashed, "_", "::")


woot <- function(expressionData, annotationData, var, chromosome) {
  squashedClassifications <- expressionData %>%
    filter(variety == var) %>%
    select(group_id, clust, clust.description, subgenome)

  squashedAnnotation <- annotationData %>%
    filter(variety == var, TriadGrp %in% squashedClassifications$group_id) %>%
    select(Chr, subgenome, Start, End, Strand, Gene, TriadGrp)

  joined <- squashedAnnotation %>%
    full_join(squashedClassifications, by = c("TriadGrp" = "group_id", "subgenome" = "subgenome"))
}







count_classification_diffs <- function(data, pair) {
  a <- data %>%
    filter(variety == pair[1]) %>%
    group_by(group_id) %>%
    summarise(classification = unique(clust.description)) %>%
    arrange(group_id)

  b <- data %>%
    filter(variety == pair[2]) %>%
    group_by(group_id) %>%
    summarise(classification = unique(clust.description)) %>%
    arrange(group_id)

  return(sum(a$classification != b$classification) / length(a$classification))
}

pairwise_classification_dist <- function(tbl, triads) {
  # Limit the full dataset to only include listed `triads`.
  filtbl <- tbl %>% filter(group_id %in% triads)
  # Sanity check the triad entries are complete and stuff is ok.
  subgenomeCheck <- filtbl %>%
    group_by(variety, group_id) %>%
    summarise(subgenome_check = length(unique(subgenome)))
  # This should always be true.
  stopifnot(subgenomeCheck$subgenome_check == 3)

  varieties <- unique(filtbl$variety)
  message(varieties)
  pairs <- combn(varieties, 2)
  diffCounts <- apply(pairs, 2, function(pair) { count_classification_diffs(filtbl, pair) })

  return(data.frame(individualA = pairs[1,], individualB = pairs[2,], dist = diffCounts))
}

find_triads_for_window <- function(annotationTbl, assembly, chromosome, winStart, winEnd) {
  genesInWindow <- annotationTbl[(annotationTbl$Start <= winEnd) &
                                 (annotationTbl$End >= winStart) &
                                   annotationTbl$race == assembly & annotationTbl$Chr == chromosome, ]
  return(unique(genesInWindow$TriadGrp))
}

pairwise_classification_dist_for_window <-
  function(classificationTbl,
           annotationTbl,
           assembly,
           chromosome,
           winStart,
           winEnd) {

    tbl <- pairwise_classification_dist(
      classificationTbl,
      find_triads_for_window(annotationTbl, assembly, chromosome, winStart, winEnd))

    tbl$chromosome <- chromosome
    tbl$windowStart <- winStart
    tbl$windowEnd <- winEnd
    return(tbl)
}

pairwise_classification_dist_for_window(combo, geneLocations, "arina", "chr1B", 1, 1000000)

find_triads_for_window(geneLocations, "arina", "chr1B", 1, 1000000)
nim <-
  pairwise_classification_dist(
  combo,
  find_triads_for_window(geneLocations, "arina", "chr1B", 1, 1000000))





