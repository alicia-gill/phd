#' Smoothing
#'
#' Scales up and smooths the observed prevalence.
#'
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#'
#' @return data frame of estimated prevalence per day.
#' @export
#'
#' @examples
#' smooth(noisy_prevalence = noisy_prev, proportion_obs = 0.2)
smooth <- function(noisy_prevalence, proportion_obs) {
  n_days <- nrow(noisy_prevalence) - 1
  prev <- noisy_prevalence[,2]

  smooth <- list("day" = 0:n_days, "prev" = 1)

  for (i in 2:n_days) {
    smooth$prev[i] <- round(sum(floor(prev[(i-1):(i+1)]/proportion_obs)/3))
  }
  smooth$prev[n_days+1] <- floor(prev[n_days+1]/proportion_obs)

  smooth <- as.data.frame(smooth)
  return(smooth)
}
