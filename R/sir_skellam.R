#' Sample Importance Resample - Skellam
#'
#' Implements an adaptive SIR algorithm and returns an approximated log-likelihood. Proposals are Skellam.
#'
#' @param n_particles number of particles used in the sampling.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param genetic_data data frame of day, number of lineages and number of coalescences.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param plot logical; if TRUE, then plots a trajectory according to weight.
#'
#' @return log-likelihood
#' @export
#'
#' @examples
#' sir_skellam(n_particles = 100, birth_rate = 0.2, death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
sir_skellam <- function(n_particles, birth_rate, death_rate, noisy_prevalence, proportion_obs, genetic_data, ess_threshold = n_particles/2, plot=F) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  samples <- matrix(nrow = N + 1, ncol = n_particles)
  weights <- matrix(nrow = N + 1, ncol = n_particles)
  resample <- rep(NA, N)
  particles <- rep(n_particles, N + 1)
  #first day always has 1
  samples[1, ] <- 1
  weights[1, ] <- 1/n_particles

  #initialise smc likelihood approximation
  int_llik <- 0
  #initialise x_resample to be 1
  x_resample <- samples[1, ]
  w <- weights[1, ]

  length <- length(genetic_data)

  for (i in 1:N) {
    #sample
    x_sample <- x_resample + extraDistr::rskellam(n_particles, birth_rate * x_resample, death_rate * x_resample)

    if (length == 0) {
      gen_llik <- 0
    } else {
      gen_llik <- dbinom(genetic_data[i + 1, 3], choose(genetic_data[i + 1, 2], 2), 2 * birth_rate / x_sample, log = T)
    }

    #compute weights
    log_weights <- gen_llik + dbinom(noisy_prevalence[i + 1, 2], x_sample, proportion_obs, log = T)

    log_weights <- ifelse(is.nan(log_weights), -Inf, log_weights)
    log_weights <- ifelse(is.na(log_weights), -Inf, log_weights)

    #if all impossible, then mission abort
    if (max(log_weights) == -Inf) {
      int_llik <- -Inf
      return(int_llik)
    }

    #normalise weights
    lse_weights <- matrixStats::logSumExp(log_weights)
    mean_weights <- matrixStats::logSumExp(log_weights + log(w))
    int_llik <- int_llik + mean_weights
    norm_weights <- exp(log_weights - lse_weights)

    #resampling
    ess <- 1 / sum(norm_weights^2)
    particles[i+1] <- round(ess)
    #if the ess is below threshold, then resample
    if (ess <= ess_threshold) {
      resample[i] <- 1
      w <- rep(1/n_particles, n_particles)
      #if n_particles==1, then sample() is weird
      if (n_particles > 1) {
        x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
      } else {
        x_resample <- x_sample
      }
    } else {
      resample[i] <- 0
      w <- norm_weights
      x_resample <- x_sample
    }

    samples[i + 1, ] <- x_resample
    weights[i + 1, ] <- w
  }

  #if plot==T, sample one trajectory according to final weight
  if (plot == T) {
    j <- sample(1:n_particles, 1, prob=norm_weights)
    sample <- samples[,j]
    lines(sample)
  }

  return(int_llik)
#  return(list("int_llik"=int_llik, "resample"=resample, "particles"=particles, "samples"=samples, "weights"=weights))
}
