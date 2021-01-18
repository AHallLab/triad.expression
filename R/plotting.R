
#'
#'
#' @export
make_expression_boxplot <- function(data) {
  boxplotFont <- element_text(face = "bold", size = 12, colour = "Black")

  boxplotTheme <- theme_light() +
    theme(legend.position = "none",
          axis.text.x = boxplotFont,
          axis.text.y = boxplotFont,
          axis.title.y = boxplotFont,
          plot.title = element_text(face = "bold",
                                    size = 12,
                                    hjust = 0.5,
                                    colour = "black"))

  dat.m <- melt(data, id.vars = 'clust.description', measure.vars = c('A', 'B', 'D'))

  dat.m$clust.description <- factor(
    dat.m$clust.description,
    levels = c("A.dominant", "B.dominant", "D.dominant", "Central",
               "A.suppressed", "B.suppressed", "D.suppressed")
  )

  palette <- c("#999999", "#9DD584", "#30BFEB", "#FFC929", "#469C3B", "#3D2374", "#F18931")

  ggarrange(ggplot(dat.m, aes(x = variable, y = value, fill = clust.description)) +
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

  data$clust.description <- factor(
    data$clust.description,
    levels = c("A.dominant", "B.dominant", "D.dominant", "Central",
               "A.suppressed", "B.suppressed", "D.suppressed")
  )

  ggtern(data, aes(x = A, y = B, z = D)) +
    geom_point(aes(colour = clust.description), alpha = 0.5) +
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
