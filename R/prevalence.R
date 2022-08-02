prevalence <- function(epidemic, stop_time, dt) {
  t <- seq(0,stop_time,by=dt)
  prev <- list("day" = t / dt, "prev" = 0)
  l <- length(t)
  for (i in 1:l) {
    ti <- t[i]
    inf <- sum(epidemic$inf_time <= ti)
    rem <- sum(epi$rem_time <= ti, na.rm = T)
    prev$prev[i] <- inf - rem
  }
  prev <- as.data.frame(prev)
  return(prev)
}
