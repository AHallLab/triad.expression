# triad.expression

## Description

triad.expression is an R package used internally by the Hall group to do certain analyses of expression with so called Triads.

## Installation

triad.expression is not published on CRAN so you will have to install triad.expression yourself from this GitHub repository. But that's ok, it's easy enough.

The easiest way to do this is to use the `install_github` command from the R devtools package like so:

```R
devtools::install_github()
```

If you do not have devtools installed, no problem, that is available through CRAN and so using

```R
install.packages("devtools") 
```

in an R session should install it for you.

Once installed, use the package with `library`.

```R
library(triad.expression)
library(dplyr) # dplyr is going to be useful too.
```







### Locating triad genes in haplotype blocks

The first step is to take the locations of triad genes, and place them into haplotype blocks.

So let's open the triad gene location data, and haplotype block location data...

```R
geneLocations <- load_gene_locations("triad_locations")
haplotypeLocations <- load_haplotype_locations("haplotype_locations.tsv")
```

Ok, once the gene and haplotype block locations are known, we can combine the two to place the triad genes
into haplotype blocks...

```R
placedTriadGenes <- place_genes_into_haplotypes(haplotypeLocations, geneLocations)
```

Now the triad genes placed into blocks, and the output from the `distance_from_centroids` function, can be combined to make a table summarizing how many triads of different types (Central, A.dominant, etc...) are present in each haplotype block.

```R
blockExpressionCounts <- count_triad_expressions_for_all_blocks(placedTriadGenes, centroidMat)
```

### Triad classifications in chromosomes as sequences

If you have annotation information for the triads of interest, as well as triads classified
according to their expression profile, these can be exported as text sequences for each chromosome,
these can then be aligned in, for example MAFFT, and visualised to see changes in order or expression.

Below is an example.

```R
# Load raw data

expression <- read.table("newest_data/expression.csv", header = T)
homology <- read.table("newest_data/homology.csv", header = T)
metadata <- read.table("newest_data/metadata.csv", sep = "\t", header = T)
geneLocations <- load_gene_locations("triad_locations")
geneLocations$subgenome <- substr(geneLocations$Chr, 5, 5)

# Compute expression categories for the triad genes for several wheat varieties.

perVarietyExpectationDistances <- lapply(unique(metadata$High.level.variety), function(v) {
  meanExpression <- triad_expression_mean_by_factor_levels(dplyr::filter(metadata, High.level.variety == v), homology, expression, "High.level.tissue")
  normalizedMeanExpression <- normalize_mean_triad_expression(meanExpression)
  normalizedMeanExpression$variety <- v
  centroidDistances <- expectation_distances(normalizedMeanExpression, "all")
  return(centroidDistances)
})
perVarietyExpectationDistances <- do.call(rbind, perVarietyExpectationDistances)

# Tidying gene locations "race" column to be consistent with perVarietyExpectationDistances's
# variety column.
geneLocations$variety <- geneLocations$race
geneLocations$variety[geneLocations$variety == "arina"] <- "ARI"
geneLocations$variety[geneLocations$variety == "jagger"] <- "JAG"
geneLocations$variety[geneLocations$variety == "julius"] <- "JUL"
geneLocations$variety[geneLocations$variety == "lancer"] <- "LER"
geneLocations$variety[geneLocations$variety == "landmark"] <- "LAN"
geneLocations$variety[geneLocations$variety == "mace"] <- "MAC"
geneLocations$variety[geneLocations$variety == "norin61"] <- "NOR"
geneLocations$variety[geneLocations$variety == "stanley"] <- "STA"

```
