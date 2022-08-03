smooth <- function(noisy_prevalence, proportion_obs) {
  n_days <- nrow(noisy_prevalence) - 1
  prev <- noisy_prevalence[,2]

  smooth <- list("day" = 0:n_days, "prev" = 1)

  smooth$prev[2] <- floor(prev[2]/proportion_obs)
  for (i in 3:n_days) {
    smooth$prev[i] <- round(sum(floor(prev[(i-1):(i+1)]/proportion_obs)/3))
  }
  smooth$prev[n_days+1] <- floor(prev[n_days+1]/proportion_obs)

  smooth <- as.data.frame(smooth)
  return(smooth)
}
