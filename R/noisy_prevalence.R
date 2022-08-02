noisy_prevalence <- function(prevalence, proportion_obs) {
  n_days <- nrow(prevalence) - 1
  noisy_prev <- list("day" = 0:n_days, "prev" = 1)
  noisy_prev$prev[2:(n_days+1)] <- rbinom(n = n_days, size = prevalence[-1,2], prob = proportion_obs)
  noisy_prev <- as.data.frame(noisy_prev)
  return(noisy_prev)
}
