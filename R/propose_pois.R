#' Propose Epidemic using Smoothed Poisson
#'
#' Simulates an epidemic by scaling up and smoothing the partially observed epidemic and sampling prevalence from a poisson distribution centred on this scaled and smoothed epidemic.
#'
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param ... allows for the use of lapply to run it multiple times.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' propose_pois(noisy_prevalence = noisy_prev, proportion_obs = 0.2, birth_rate = 0.15, death_rate = 0.1, ptree = sample_tree)
propose_pois <- function(noisy_prevalence, proportion_obs, birth_rate, death_rate, ptree, ...) {
  smooth_prev <- smooth(noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs)
  lambda <- 0.1 + smooth_prev[-1,2]

  n_days <- nrow(noisy_prevalence) - 1

  prev_new <- rep(1, n_days + 1)
  prev_new[-1] <- smooth_prev[-1,2] +
                  extraDistr::rtpois(n_days, lambda = lambda, a = noisy_prevalence[-1,2] + round(lambda) - smooth_prev[-1,2] - 1) -
                  round(lambda)

  prev_new <- data.frame("day" = 0:n_days, "prev" = prev_new)

  f_hat <- skellam_llik(birth_rate = birth_rate, death_rate = death_rate, prevalence = prev_new) +
           genetic_llik(birth_rate = birth_rate, ptree = ptree, prevalence = prev_new, stop_time = n_days) +
           sum(dbinom(x = noisy_prevalence[-1,2], size = prev_new[-1,2], prob = proportion_obs, log = TRUE))
  q <- sum(extraDistr::dtpois(prev_new[-1,2] - smooth_prev[-1,2] + round(lambda), lambda = lambda, a = noisy_prevalence[-1,2] + round(lambda) - smooth_prev[-1,2] - 1, log = TRUE))

  return(list("f_hat"=f_hat, "q"=q))
}
