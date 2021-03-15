
#'
#' @export
load_gene_locations <- function(folder) {
  files <- list.files("triad_locations")
  files <- files[grepl(".tsv", files)]
  triadTable <- do.call(rbind, lapply(files, function(x) {
    tbl <- as_tibble(read.delim(paste0(folder, '/', x)))
    dtbl <- distinct(tbl)
    if(nrow(tbl) > nrow(dtbl)) {
      warning(paste0("Duplicate entries were detected in ", x, " these will be removed"))
    }
    dtbl$race <- unlist(strsplit(x, "_"))[1]
    return(dtbl)
  }))
  checkColumns(c("Chr", "Start", "End", "Strand", "Gene", "TriadGrp", "race"), colnames(triadTable), "triad location")
  return(triadTable)
}

#'
#' @export
load_haplotype_locations <- function(filepath) {
  haplotypeBlocks <- read.table(filepath, header = TRUE)
  # Limit haplotype blocks to just the ones that appear in the triad locations data.
  return(haplotypeBlocks)
}

place_genes_in_blocks_for_assembly <- function(blocks, triads, assm) {
  blocks <- blocks[blocks$assembly == assm, ]
  triads <- triads[triads$race == assm, ]
  blockRanges <- makeGRangesFromDataFrame(blocks, keep.extra.columns = TRUE)
  triadRanges <- makeGRangesFromDataFrame(triads,
                                          start.field = "Start",
                                          end.field = "End",
                                          strand.field = "Strand",
                                          keep.extra.columns = TRUE)

  # Identify triads that overlap a haplotype block.
  overlaps <- as.data.frame(findOverlaps(triadRanges, blockRanges, ignore.strand = T))

  t <- triads[overlaps$queryHits, ]
  b <- blocks[overlaps$subjectHits, ]#c("start",
                                      #"end",
                                      #"block_no",
                                      #"chr_length",
                                      #"merged_block")]
  #colnames(b) <- c("block_start", "block_end", "block_no", "chr_length", "merged_block")
  r <- cbind(t, b)
  return(r)
}

#'
#' @export
place_genes_into_haplotypes <- function(haplotypeBlockTable, geneLocationTable) {
  # Limit datasets to the genomes / races present in both gene data, and haplotype data.
  genomes <- intersect(unique(haplotypeBlockTable$assembly), unique(geneLocationTable$race))
  haplotypeBlockTable <- haplotypeBlockTable[haplotypeBlockTable$assembly %in% genomes, ]
  geneLocationTable <- geneLocationTable[geneLocationTable$race %in% genomes, ]

  placedGenes <- do.call(rbind, lapply(genomes, function(genome) {
    place_genes_in_blocks_for_assembly(haplotypeBlockTable, geneLocationTable, genome)
  }))

  return(placedGenes)
}

#'
#' @export
count_triad_expressions_for_block <- function(gene_haplo_tbl, triad_expression_matrix, block) {
  triad_groups <- gene_haplo_tbl[gene_haplo_tbl$block_no == block, "TriadGrp"]
  block_expressions <- triad_expression_matrix[triad_expression_matrix$group_id %in% triad_groups, ]
  block_expression_tbl <- table(block_expressions$clust.description)
  df <- data.frame(block = block,
                   Central = block_expression_tbl["Central"],
                   A.dominant = block_expression_tbl["A.dominant"],
                   B.dominant = block_expression_tbl["D.dominant"],
                   D.dominant = block_expression_tbl["D.dominant"],
                   A.suppressed = block_expression_tbl["A.suppressed"],
                   B.suppressed = block_expression_tbl["B.suppressed"],
                   D.suppressed = block_expression_tbl["D.suppressed"], row.names = NULL)
  df[is.na(df)] <- 0
  return(df)
}

#'
#' @export
count_triad_expressions_for_all_blocks <- function(gene_haplo_tbl, triad_expression_matrix) {
  blocks <- unique(gene_haplo_tbl$block_no)
  return(do.call(rbind, lapply(blocks, function(block) {
    count_triad_expressions_for_block(gene_haplo_tbl, triad_expression_matrix, block)
  })))
}

