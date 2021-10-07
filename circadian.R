library(dplyr)

data("triad_homology")
rawTable <- read.csv("inst/extdata/circadian/q_values.csv", header = T)
csLocations <- read.table("inst/extdata/cs_Gene_loc.csv", header = T) # Is actually a tsv.
triadsFlat <- triad.expression:::reshape_homology_triads(triad_homology)

joinedTable <- csLocations %>%
  select(gene, chr, start, end, strand, variety) %>%
  inner_join(rawTable, by = c("gene" = "CycID")) %>%
  inner_join(triadsFlat, by = "gene") %>%
  mutate(
    classification = case_when(
      meta2d_BH.Q > 0.05 ~ "arrhythmic",
      meta2d_BH.Q < 0.01 ~ "rhythmic"
    )
  ) %>%
  group_by(group_id) %>%
  mutate(
    n_arrythmic = sum(classification == "arrhythmic", na.rm = T),
    n_rhythmic = sum(classification == "rhythmic", na.rm = T)
  ) %>%
  arrange(group_id)

joinedTable %>%
  filter(n_arrythmic > 0, grepl("1", chr)) %>%
  loom_plot("start", "subgenome", "group_id", "n_arrythmic", T, yLab = "Subgenome")

runTbl <- joinedTable %>%
  filter(n_arrythmic > 0, grepl("1", chr)) %>%
  group_by(subgenome) %>%
  arrange(start) %>%
  group_map(~ {
    o <- compute_runs(.x, "n_arrythmic")
    o$subgenome <- .y$subgenome
    o
  }) %>% bind_rows







runTbl %>%
  filter(runLength > 1) %>%
  ggplot(aes(x = startBp, y = subgenome)) +
  geom_point(aes(color = as.factor(value))) +
  scale_x_continuous(labels = function(x) format(x / 1000000))

loom_plot(runTbl2, "startBp", "subgenome", "group_id", "n_arrythmic", T, yLab = "Subgenome", colourLab = "n_arrythmic")



runTbl2 <- joinedTable %>%
  filter(n_arrythmic > 0, grepl("1", chr)) %>%
  group_by(subgenome) %>%
  arrange(start) %>%
  group_map(~ {
    o <- merge_runs(.x, compute_runs(.x, "n_arrythmic"))
    o$subgenome <- .y$subgenome
    o
  }) %>% bind_rows








runTbl2 %>%
  filter(n_arrythmic.runLength > 2) %>%
  loom_plot("start", "subgenome", "group_id", "n_arrythmic", T, yLab = "Subgenome", colourLab = "n_arrythmic")





ggpubr:::ggarrange(
  joinedTable %>%
    filter(n_arrythmic > 0, grepl("1", chr)) %>%
    ggplot(aes(x = start, y = subgenome, group = group_id, color = as.factor(n_arrythmic))) +
    geom_path() +
    geom_point(aes(color = as.factor(n_arrythmic))) +
    scale_x_continuous(labels = function(x) format(x / 1000000)) +
    xlab("Position (Mb)") +
    ylab("Subgenome")
)

