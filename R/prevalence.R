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
