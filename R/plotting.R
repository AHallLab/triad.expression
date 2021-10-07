
#'
#'
#' @export
make_expression_boxplot <- function(data) {
  boxplotFont <- ggplot2::element_text(face = "bold", size = 12, colour = "Black")

  boxplotTheme <- ggplot2::theme_light() +
    theme(legend.position = "none",
          axis.text.x = boxplotFont,
          axis.text.y = boxplotFont,
          axis.title.y = boxplotFont,
          plot.title = ggplot2::element_text(face = "bold",
                                    size = 12,
                                    hjust = 0.5,
                                    colour = "black"))

  dat.m <- reshape2::melt(data, id.vars = 'clust.description', measure.vars = c('A', 'B', 'D'))

  dat.m$clust.description <- factor(
    dat.m$clust.description,
    levels = c("A.dominant", "B.dominant", "D.dominant", "Central",
               "A.suppressed", "B.suppressed", "D.suppressed")
  )

  palette <- c("#999999", "#9DD584", "#30BFEB", "#FFC929", "#469C3B", "#3D2374", "#F18931")

  ggarrange(ggplot(dat.m, aes(x = .data$variable, y = .data$value, fill = .data$clust.description)) +
    geom_boxplot() +
    boxplotTheme +
    labs(title = NULL,
         x = NULL,
         y = "Relative contribution (%)") +
    scale_y_continuous(labels = seq(0, 100, by = 25)) +
    scale_fill_manual(values = palette) +
    facet_wrap(~clust.description, nrow = 2))
}

#'
#'
#' @export
make_expression_ternplot <- function(data) {
  palette <- c("#999999", "#9DD584", "#30BFEB", "#FFC929", "#469C3B", "#3D2374", "#F18931")

  data <- data %>%
    select(.data$A, .data$B, .data$D, .data$clust.description) %>%
    distinct()

  ggtern(data, aes(x = .data$A, y = .data$B, z = .data$D)) +
    geom_point(aes(colour = .data$clust.description), alpha = 0.5) +
    theme_bw() +
    theme(legend.text = element_text(size = 12),
          axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
          axis.title.y = element_text(face = "bold", size = 12, colour = "black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.key = element_rect(colour = "white", fill = NA)) +
    scale_colour_manual(values = palette) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    xlab("A") + ylab("B") + zlab("D") +
    theme_arrownormal() +
    theme_ticklength(major = unit(5.0,'mm'), minor = unit(2.5,'mm'))
}


#'
#'
#' @export
make_expression_figure <- function(data) {
  grid.arrange(make_expression_ternplot(data),
               make_expression_boxplot(data), ncol = 2)
}


#'
#' @export
loom_plot <- function(distWithAnnotation, positionVar, yVar, geneVar, colourVar, colourLines = F, xLab = "Position (Mb)", yLab = "Subgenome", colourLab = "Expression Balance", arrangeWrap = T) {
  pFrame <- distWithAnnotation[,c(positionVar, yVar, geneVar, colourVar)]
  colnames(pFrame) <- c("xx", "yy", "gg", "cc")
  pFrame$cc <- as.factor(pFrame$cc)

  if(colourLines) {
    pathGeom <- geom_path(aes(color = .data$cc))
  } else {
    pathGeom <- geom_path()
  }

  p <- ggplot(pFrame, aes(x = .data$xx, y = .data$yy, group = .data$gg)) +
                   pathGeom +
                   geom_point(aes(color = .data$cc)) +
                   xlab(xLab) +
                   ylab(yLab) +
                   scale_color_brewer(name = colourLab, palette = "Set1") +
                   scale_x_continuous(labels = function(x) format(x / 1000000))

  if(arrangeWrap) {
    return(ggarrange(p))
  } else {
    return(p)
  }
}

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
