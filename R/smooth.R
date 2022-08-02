smooth <- function(noisy_prevalence, proportion_obs) {
  n_days <- nrow(noisy_prevalence) - 1
  prev <- noisy_prevalence[,2]

  smooth <- list("day" = 0:n_days, "prev" = 1)

  for (i in 2:n_days) {
    smooth$prev[i] <- round(sum(prev[(i-1):(i+1)]/proportion_obs)/3)
  }
  smooth$prev[n_days+1] <- prev[n_days+1]

  smooth <- as.data.frame(smooth)
  return(smooth)
}
