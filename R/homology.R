#' @export
validateHomologyTable <- function(tbl) {
  checkColumns(c("group_id", "A", "B", "D"), colnames(tbl), "homology table")
}

#' @export
read_homology <- function(filename = triadInfo$homology_path) {
  tbl <- read.table(filename, header = T)
  validateHomologyTable(tbl)
  return(tbl)
}

reshape_homology_triads <- function(homology) {
  do.call(rbind, lapply(c("A", "B", "D"), function(subgenome){
    tibble(gene = homology[[subgenome]], subgenome = subgenome, group_id = homology$group_id)
  }))
}
