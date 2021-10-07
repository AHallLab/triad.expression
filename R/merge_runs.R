merge_runs <- function(data, runData) {
  runVariable <- attributes(runData)$runVariable
  newRunDataNames <- paste(
    runVariable,
    c("runID", "runStartBp", "runEndBp", "runValue", "runStart", "runEnd", "runLength"),
    sep="."
  )

  if(!(runVariable %in% colnames(data))) {
    stop("variable '", runVariable, "' not found in data")
  }

  runID <- rep(1:nrow(runData), runData$runLength)

  if(length(runID) != nrow(data)) {
    stop("Computed ", length(runID), " run IDs for a table with ", nrow(data), " rows.")
  }

  data[[newRunDataNames[[1]]]] <- runID
  data[[newRunDataNames[[2]]]] <- runData$startBp[runID]
  data[[newRunDataNames[[3]]]] <- runData$endBp[runID]
  data[[newRunDataNames[[7]]]] <- runData$runLength[runID]

  return(data)
}
