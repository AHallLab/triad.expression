# Read in input data sources
metadata <- read_metadata("metadata.csv")
homology <- read_homology("HCTriads.csv")
expression <- read_csv("expression_data_count.csv")

# Do some pre-filtering, may be custom per project / lab / whatever

genes_to_use <- (function(){
  expressed_genes <- read.csv("./expressed_genes_tpmsOver0.5AtLeast8Samples.csv")
  colnames(expressed_genes) <- c("gene", "tpm", "count")
  genes_to_use <- data.frame(gene = expressed_genes$gene)
  return(genes_to_use)
})()

homology <- homology %>%
  filter(synteny == "segmental homeologs" &
           cardinality_abs == '1:1:1' &
           (A %in% genes_to_use$gene | B %in% genes_to_use$gene | D %in% genes_to_use$gene))

expression <- expression %>%
  filter(gene %in% c(as.character(homology$A),
                     as.character(homology$B),
                     as.character(homology$D)))

# Get mean expression for each triad, for each level of a factor
meanExpression <- triad_expression_mean_by_factor_levels(metadata, homology, expression, "High.level.tissue")

normalizedMeanExpression <- normalize_mean_triad_expression(meanExpression)

controidDistances <- distance_from_centroid(normalizedMeanExpression, "all")
