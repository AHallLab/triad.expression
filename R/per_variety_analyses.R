
#'
#' @export
pv_xy_plot <- function(pvAnnotatedDistances, varX, varY, chrom) {
  data <- do.call(rbind, pvAnnotatedDistances[c(varX, varY)]) %>%
    select(.data$group_id, .data$clust, .data$gene, .data$subgenome, .data$variety, .data$Chr, .data$Start, .data$End, .data$Strand, .data$Gene) %>%
    filter(.data$Chr == chrom) %>%
    arrange(.data$group_id)
  plotTbl <- data %>%
    group_by(.data$group_id) %>%
    summarize(diff = .data$clust[1] != .data$clust[2], x = .data$Start[1], y = .data$Start[2]) %>%
    na.omit()
  pvPlot <- ggplot(plotTbl, aes(x = .data$x, y = .data$y, colour = .data$diff, shape = .data$diff, alpha = .data$diff)) +
    geom_point()
  return(
    list(
      plt = pvPlot,
      data = data,
      plot_tbl = plotTbl
    )
  )
}
