#' @importFrom fields rdist
#' @importFrom reshape2 acast melt
#' @importFrom readr read_csv read_table
#' @import dplyr
#' @import ggplot2
#' @import ggtern
#' @import ggpubr
#' @import GenomicRanges
NULL


#' Wheat expression experiment data for 210 samples
#'
#' A data table containing expression levels of 51108 genes for 210 samples,
#' derived from RNA-seq experiments.
#'
#' @format A data frame with 51108 rows and 210 variables:
#' \describe{
#'   \item{gene}{name of the gene}
#'   \item{Sample_*}{A column of gene expression values for the sample corresponding to the column name}
#'   ...
#' }
#'
"expression"
