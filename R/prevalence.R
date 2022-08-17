#' Prevalence
#'
#' Turns simulated epidemic into a vector of cases per day.
#'
#' @param epidemic The output from epidemic.R.
#' @param stop_time The number of days the epidemic simulation was run for.
#'
#' @return A data frame with 2 columns and N+1 rows. Column 1 is day and column 2 is number of cases.
#' @export
#'
#' @examples
#' prevalence(epidemic = epi, stop_time = 50)
prevalence <- function(epidemic, stop_time) {
  t <- seq(0,stop_time,by=1)
  prev <- list("day" = t, "prev" = 0)
  l <- length(t)
  for (i in 1:l) {
    ti <- t[i]
    inf <- sum(epidemic$inf_time <= ti)
    rem <- sum(epidemic$rem_time <= ti, na.rm = T)
    prev$prev[i] <- inf - rem
  }
  prev <- as.data.frame(prev)
  return(prev)
}
