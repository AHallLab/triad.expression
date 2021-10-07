#' Compute runs of a variable in a triad data table.
#'
#' @param data A data frame or tibble that contains triad data.
#'
#' @export
compute_runs <- function(data, var) {
  rr <- rle(data[[var]])
  rrEnds <- cumsum(rr$lengths)
  rrStarts <- rrEnds - rr$lengths + 1
  runData <- data.frame(
    startBp = data$start[rrStarts],
    endBp = data$end[rrEnds],
    value = rr$values,
    runStart = rrStarts,
    runEnd = rrEnds,
    runLength = rr$lengths
  )
  attributes(runData)$runVariable <- var
  return(runData)
}
