
# How the raw data files were processed.

# Expression data

prepExpressionData <- function(filepath) {
  data <- read.table(filepath, header = T)
  gene <- rownames(data)
  rownames(data) <- NULL
  return(cbind(gene, data))
}

# expression_data <- prepExpressionData("inst/extdata/170521_pantrans_commonref_expression.csv")
# save(expression_data, file = "data/expression_data.RData")

# Triad Homology


validateHomologyTable <- function(tbl) {
  checkColumns(c("group_id", "A", "B", "D"), colnames(tbl), "homology table")
}


prepHomology <- function(filepath) {
  tbl <- read.table(filepath, header = T)
  validateHomologyTable(tbl)
  return(tbl)
}

#triad_homology <- prepHomology("inst/extdata/170521_pantrans_commonref_homology.csv")
