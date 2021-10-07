#' @importFrom fields rdist
#' @importFrom reshape2 acast melt
#' @importFrom readr read_csv read_table
#' @importFrom ggpubr ggarrange
#' @import dplyr
#' @import ggplot2
#' @import GenomicRanges
#' @import RColorBrewer
#' @import ggtern
NULL

#' An example gene expression data-set.
#'
#' A dataset containning gene expression results of 102,442 genes
#' from 232 samples.
#'
#' Each row represents a single gene.
#' The column 'gene', has the name of the gene and all other columns
#' are separate samples for which the levels of expression of each gene
#' were measured.
#'
#' Note: The data comes from many varieties of wheat, but in this data-set,
#' the gene names are all Chinese Spring based, for all varieties.
#' So for example, values in the data-set for gene "TraesCS3A02G200800",
#' from a sample taken from the variety ARI, simply correspond to the gene
#' in ARI that is the homologue of the Chinese Spring gene "TraesCS3A02G200800".
#'
#' @docType data
#' @usage data(expression_data)
#' @format A data frame with 102442 rows and 233 variables.
#' @keywords datasets
"expression_data"

#' An example gene annotation data-set
#'
#'
#' Each row corresponds to a single gene's location in a single genome of a variety of wheat.
#'
#' \itemize{
#'   \item chr. chromosome the gene exists in
#'   \item start. the base position in the genome, of the first base in the gene
#'   \item end. the base position in the genome, of the final base in the gene
#'   \item strand. which strand the gene exists on. Either '+' or '-'
#'   \item gene. the variety specific name of the gene
#'   \item group_id. the ID of the triad the gene belongs in
#'   \item variety. the variety of wheat
#'   \item subgenome. the subgenome of wheat the gene exists in ('A', 'B', or 'D')
#' }
#'
#' @docType data
#' @usage data(gene_locations)
#' @format A tibble with 318649 rows and 8 variables
#' @keywords datasets
"gene_locations"

#' An example region dataset
#'
#' An example dataset of some genomic regions identified in this paper [todo].
#'
#' The variables are:
#'
#' \itemize{
#'   \item chr. chromosome the region exists in
#'   \item start. the base position in the genome of the first base of the region
#'   \item end. the base position in the genome of the final base of the region
#'   \item High.level.variety. The wheat variety this region was identified in
#' }
#'
#' @docType data
#' @usage data(wicker_te_regions)
#' @format A data frame with 320 rows and 4 variables
#' @keywords datasets
"wicker_te_regions"




varietyInfo <- new.env()

varietyInfo$tags <- c("ARI", "JAG", "JUL", "LER", "LAN", "MAC", "NOR", "STA")

varietyInfo$triad_location_files <- list(
  ARI = "arina_noChrUn_rhb_triads.tsv",
  JAG = "jagger_noChrUn_rhb_triads.tsv",
  JUL = "julius_noChrUn_rhb_triads.tsv",
  LER = "lancer_noChrUn_rhb_triads.tsv",
  LAN = "landmark_noChrUn_rhb_triads.tsv",
  MAC = "mace_noChrUn_rhb_triads.tsv",
  NOR = "norin61_noChrUn_rhb_triads.tsv",
  STA = "stanley_noChrUn_rhb_triads.tsv"
)

#   list(
#   ARI = system.file("extdata", "triad_locations", "arina_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   JAG = system.file("extdata", "triad_locations", "jagger_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   JUL = system.file("extdata", "triad_locations", "julius_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   LER = system.file("extdata", "triad_locations", "lancer_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   LAN = system.file("extdata", "triad_locations", "landmark_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   MAC = system.file("extdata", "triad_locations", "mace_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   NOR = system.file("extdata", "triad_locations", "norin61_noChrUn_rhb_triads.tsv", package = "triad.expression"),
#   STA = system.file("extdata", "triad_locations", "stanley_noChrUn_rhb_triads.tsv", package = "triad.expression")
# )
