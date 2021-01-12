#' @import tibble
#' @import readr

#' @export
validateHomologyTable <- function(tbl) {
  checkColumns(c("group_id", "A", "B", "D"), colnames(tbl), "homology table")
}

#' @export
read_homology <- function(filename) {
  tbl <- read_csv(filename)
  validateHomologyTable(tbl)
  return(tbl)
}

reshape_homology_triads <- function(homology) {
  do.call(rbind, lapply(c("A", "B", "D"), function(subgenome){
    tibble(gene = homology[[subgenome]], subgenome = subgenome, group_id = homology$group_id)
  }))
}

get_means_df <- function(metadata, tpms, triads, type = "High.level.tissue"){

  # BJW, a subfunction that can be run for any `type` and `factor`.
  # When combined with use of lapply, removes a lot of code repetition and bloat where bugs may lurk.
  meanMaker <- function(type, fac) {
    samples <- getSamplesForFactor(metadata, type = type, factor = fac)
    # In the origional function, there's an issue here when trying to subset `tpms` with `samples`,
    # as if `samples` contains any names not existing in the cols of `tpms`, this throws an error.
    # My fix is to ensure to only keep samples that exist as colnames in tpms.
    # i.e. samples[samples %in% colnames(tpms_for_triads)]
    samples_in_tpms <- samples[samples %in% colnames(tpms)]
    if(length(samples_in_tpms) > 1){
      meanvals <- rowMeans(tpms[, samples_in_tpms])
    } else {
      meanvals <- tpms[, samples_in_tpms]
    }

    out <- data.frame(
      value = meanvals,
      factor = fac,
      gene = rownames(tpms),
      samples = length(samples_in_tpms),
      stringsAsFactors = FALSE
    )

    return(out)
  }

  # Use meanMaker on tpms for `type` and `fac` set to "all", followed by using
  # meanMaker on tpms for each level of column `type` in the metadata.
  # Combine all tables produced with do.call and rbind.
  type_factors <- levels(metadata[,type])
  values <- do.call(rbind, c(
    list(meanMaker("all", "all")),
    lapply(type_factors, function(fac) meanMaker(type, fac))
  ))

  tmp <- values %>%
    filter(factor != "all") %>%
    group_by(gene) %>%
    summarise(
      value = mean(value),
      samples = n_distinct(factor),
      factor = "all_means"
    )

  values <- rbind(values, tmp)

  triads_flat <- reshape_triad_groups(triads)

  vals <- values %>% left_join(triads_flat, by = "gene")

  return(vals)
}
