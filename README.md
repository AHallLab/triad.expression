# triad.expression

## Description
triad.expression is an R package used internally by the Hall group to do certain analyses of expression with
so called Triads.

## Installation

triad.expression is not published on CRAN so you will have to install triad.expression yourself from this
GitHub repository.

The easiest way to do this is to use the `install_github` command from the R devtools package like so:

```R
devtools::install_github()
```

If you do not have devtools installed, no problem, that is available through CRAN and so using

```R
install.packages("devtools") 
```

in an R session should install it for you.


## Example use

### Computing expression categories of triads.

First things first, load the package:

```R
library(triad.expression)
```

The first steps require you to load data from a csv file.
3 tables are nessesery, expression data table, a homology table, and a metadata table.

The format and requirements of these tables are given in detail in the package documentation.

You load them using the `read_metadata`, `read_expression`, and `read_homology` functions.

However, this package also has some examples ones included, these can just be loaded with the `data`
command, so let's just do that:

```R
data(expression)
data(homology)
data(metadata)
```

Next you want to compute the mean expression values for each gene, for every level of a factor.

For example, in the example metadata there is a column called "High.level.tissue", with a few levels.
For each level of that column ("roots", "leaves/shoots", etc.) compute the mean expression levels
for all the genes, using the samples for that level.

This is accomplished with the `triad_expression_mean_by_factor_levels` function.
It takes as input, your metadata table, expression table, and homology table.

```R
meanExpression <- triad_expression_mean_by_factor_levels(metadata, homology, expression, "High.level.tissue")
```

The resulting table contains the mean expression of each gene, for every factor level, and for the entire dataset; see the level column has an "all" category, as well as "roots" etc.

Next we have to normalize the data, which is done with `normalize_mean_triad_expression`.

```R
normalizedMeanExpression <- normalize_mean_triad_expression(meanExpression)
```

Once you have normalized values, you can compute the distance of each triad from the data's centroid,
and from the centroids of a set of 7 expected triad expression categories:

- Central
- A dominant
- B dominant
- C dominant
- A suppressed
- B suppressed
- D suppressed

This is done with the `distances_from_centroids` function, which must be given the data, and a factor level.
The example below computes the distances for the triad expression derived from the root samples.
Doing this for multiple factor levels would allow you to see any differences between tissues for a given triad.

```R
centroidDistances <- distances_from_centroids(normalizedMeanExpression, "roots")

centroidMat <- centroidDistances[[1]]
```

Distances can be plotted in a figure of a ternary plot an accompanying boxplots with the `make_expression_figure` function:

```R
make_expression_figure(centroidMat)
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


