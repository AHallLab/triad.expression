# triad.expression

## Description

The `triad.expression` R package seeks to provide methods for
analyzing and visualizing information such as gene expression balance
in sets of homologous genes in the Wheat genome called Triads.

A Triad is a set of three homologous (related in position and
function) genes where a gene exists in the A, B, and D sub genomes of
Wheat.

Each of the three genes in a triad may be expressed at different
levels in different tissues or in different Wheat varieties. This is
sometimes referred to as "Triad Balance".

`triad.expression` allows you to compute mean gene expression levels 
of triads across combinations of factors such as experimental
condition, or tissue, allowing the detection of changes in gene 
expression between levels of those factors.
From such data "Triad Balance" can then be computed, and visualized 
using ternary plots or loom plots.

### Features Summary

- Example data sets
- Compute mean gene expression across factors
- Calculate and assign "Triad Balance"
- Visualize Triad data across chromosomes with "Loom Plots"
- Visualize Triad Balance of all triads with ternary plots

## Installation

triad.expression is not currently published on CRAN, so you will have
to install `triad.expression` yourself from this GitHub repository.
But that's ok, it's easy enough.

The easiest way to do this is to use the `install_github` command
from the R devtools package like so:

```R
devtools::install_github()
```

If you do not have devtools installed, no problem, that is available
through CRAN and so using

```R
install.packages("devtools") 
```

in an R session should install it for you.

Once installed, use the package with `library`.

```R
library(triad.expression)
library(dplyr) # dplyr is going to be useful too.
```








