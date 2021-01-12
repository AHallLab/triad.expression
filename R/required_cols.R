

checkColumns <- function(required, present, tbl_type) {
  missing <- required[!(required %in% present)]
  if(length(missing)){
    stop(
      paste0("a ",
             tbl_type,
             " requires columns: ",
             paste0(required, collapse = ", "),
             "\n",
             "missing column(s): ",
             paste0(missingColumns, collapse = ", ")
      )
    )
  }
}
