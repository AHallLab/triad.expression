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
