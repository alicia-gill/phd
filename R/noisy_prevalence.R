#' Noisy Prevalence
#'
#' Generates a partially observed prevalence.
#'
#' @param prevalence data frame of prevalence per day.
#' @param proportion_obs proportion of cases observed.
#'
#' @return data frame of observed prevalence per day.
#' @export
#'
#' @examples
#' noisy_prevalence(prevalence = prev, proportion_obs = 0.2)
noisy_prevalence <- function(prevalence, proportion_obs) {
  n_days <- nrow(prevalence) - 1
  noisy_prev <- list("day" = 0:n_days, "prev" = 1)
  noisy_prev$prev[2:(n_days+1)] <- rbinom(n = n_days, size = prevalence[-1,2], prob = proportion_obs)
  noisy_prev <- as.data.frame(noisy_prev)
  return(noisy_prev)
}
