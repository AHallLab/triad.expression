---
title: 'triad.expression: Analysis and visualization of Wheat gene triad expression patterns'
tags:
  - R
  - bioinformatics
  - gene expression
  - wheat
  - genomics
authors:
  - name: Ben J. Ward^[corresponding author]
    orcid: 0000-0001-6337-5238
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Earlham Institute
   index: 1
date: 29 September 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The `triad.expression` R package seeks to provide methods for analyzing and visualizing
patterns of gene expression in sets of homologous genes in the Wheat genome called Triads.
A Triad is a set of three homologous (related in position and function) genes where a
gene exists in the A, B, and D sub genomes of Wheat.
Each of the three genes in a triad may be expressed at different levels in different
tissues or in different Wheat varieties.
`triad.expression` allows you to compute mean gene expression levels of triads across
combinations of factors such as experimental condition, or tissue, allowing the
detection of changes in gene expression between levels of those factors.
Using the mean expression values computed for triads `triad.expression` can categorize
triads into one of seven expression categories, visualize expression levels of triads
as ternary plots, and visualize triad expression categories across chromosomes and
between factor levels such as variety.


# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Features of `triad.expression`

All these features are documented in more detail in the package vignette.

## Included data

The package has some data files included that have been used in past projects.
All of them or (more likely) some of them can be swapped out for your own novel data-sets.

### Triad homology
Triads are made up of three homologous genes from the A, B, and D sub-genomes of wheat.  
This package makes use of, and includes, a data set that describes these homologies.
To load it for use, simply use `data(triad_homology)`.

### Gene locations
Genomic locations (chromosome, start position, end position, strand etc.) of the genes in each triad.
This data is included in this package for several varieties of wheat.
This data is loaded for use through the `data(gene_locations)` command.

### Gene expression data
Analyses in the package make use of gene expression data for the triads.
This data includes both the raw expression values obtained from experiment,
as well as a table of metadata about the samples and the experiments they came from.

The metadata is accessed using `data(expression_metadata)`, and expression level data is
accessed using `data(expression_data)`.

## Analyses

## Computing mean gene expression values and classifying triads
`triad.expression` can compute the mean expression values for each triad gene, for
every level of combination of factors, through the use of a function called `triad_expression_means_by_factors`.
This function takes the expression data sets, and the triad homology information, and
a list of metadata column names to use as the factors.
For example, in the `expression_metadata` included with the package, there is a column
called "High.level.variety", with a level for each variety of wheat. This can be used
as a factor.

```R
library(triad.expression)

data(expression_data)
data(expression_metadata)
data(triad_homology)

meanExpressionByVarietyAndTissue <-
    triad_expression_means_by_factors(expression_data,
                                      expression_metadata,
                                      triad_homology,
                                      c("High.level.variety", "High.level.tissue"))
```

Once these values have been computed, you can normalize the data, which is done with the
`normalize_triad_expression_means` function.

```R
normalizedMeans <- 
    normalize_triad_expression_means(meanExpressionByVarietyAndTissue)
```

Once you have normalized mean triad expression data, using the `centroid_distances` function,
you can compute each triads distance from the centroids, and from the centroids, assign to each
triad one of 7 expression categories:

- Central
- A dominant
- B dominant
- C dominant
- A suppressed
- B suppressed
- D suppressed

```R
cDistances <- centroid_distances(normalizedMeans)
```

## Plotting

### Ternary plots

Once you have computed centroid distances `triad.expression` helps you produce ternary plots
through the `make_expression_ternplot` function.
Some filtering of the computed centroid distances data may be necessary, as in the example below.

```R
library(dplyr)
cDistances %>%
  # We only want to plot the data for aerial organ
  # tissue and the ARI variety.
  filter(High.level.tissue == "aerial organs",
         High.level.variety == "ARI") %>%
  make_expression_ternplot
```

The resulting ternary plots look like the example in \autoref{fig:exampletern}:

![An example ternary plot.\label{fig:exampletern}](example_ternaryplot.png)


### Loom plots

Structural rearrangements of the genome that have occurred during Wheat's evolutionary history,
may be linked to changes in a Triad's expression pattern. In addition so could other spatial factors
such as where a triad exists in the genome.

To help visualize such phenomena `triad.expression` provides a function called `loom_plot`.
Depending on how you preprocess the data you can visualize a few different phenomena with it.

The first step is typically to merge your computed triad distances with the gene_locations data provided in
`triad.expression`, or your own:

```R
data(gene_locations)
distWithLocation <- join_distances_and_annotation(cDistances, gene_locations)
```

From here, this data can be pre-processed and passed to the `loom_plot` function.
In this case, the data is limited to only one chromosome, 2A, and to the root tissue.
Then it is filtered to exclude triads that don't ever differ in expression pattern between
the different Wheat varieties.

The `loom_plot` function accepts a variable to use as the x-axis, y-axis, for drawing the connecting lines, and for color.
In the example blow the x-axis is the start position of every triad gene, High.level.variety is the y-axis, group_id (triad)
is used to draw the connecting lines, and clust.description (expression pattern) is used for color.

```R
distWithLocation %>%
  # Limit to Chromosome 2A, and the root tissue.
  filter(chr == "chr2A", High.level.tissue == "root")  %>%
  # Exclude triads that never differ in expression between
  # the different wheat varieties.
  group_by(group_id) %>%
  filter(n_distinct(clust.description) > 1) %>%
  # Plot!
  loom_plot("start",
            "High.level.variety",
            "group_id",
            "clust.description",
            yLab = "Variety")
```

The resulting loom plot looks like the example in \autoref{fig:exampleloom}:

![An example loom plot.\label{fig:exampleloom}](example_loomplot.png)


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge the work of Dr Ricardo H Ramírez González from the John Innes Centre, who's supplementary scripts xxx formed the initial inspiration for this package.

# References
