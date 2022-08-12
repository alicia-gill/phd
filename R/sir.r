sir <- function(n_particles, birth_rate, death_rate, proportion_obs, noisy_prevalence, ptree) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = N)
  samples <- matrix(nrow = N + 1, ncol = n_particles)
  #first day always has 1
  samples[,1] <- 1

  #initialise smc likelihood approximation
  int_llik <- 0

  #Day 1
  #sample from a truncated poisson
  x_sample <- extraDistr::rtpois(n_particles, 1 + noisy_prevalence[2, 2] / proportion_obs, a = 0)
  log_weights <- extraDistr::dskellam(x_sample - 1, birth_rate, death_rate, log = T) +
                 dbinom(genetic_data[2, 3], choose(genetic_data[2, 2], 2), 2 * birth_rate / x_sample, log = T) +
                 dbinom(noisy_prevalence[2, 2], x_sample, proportion_obs, log = T) -
                 extraDistr::dtpois(x_sample, 1 + noisy_prevalence[2, 2] / proportion_obs, a = 0, log = T)
  lse_weights <- matrixStats::logSumExp(log_weights)
  log_n_particles <- log(n_particles)
  int_llik <- int_llik + lse_weights - log_n_particles
  #normalise weightsb
  norm_weights <- exp(log_weights - lse_weights)
  if (sum(norm_weights, na.rm = T) != 1) {
    norm_weights <- rep(1/n_particles, n_particles)
  }
  if (n_particles > 1) {
    x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
  } else {
    x_resample <- x_sample
  }

  #Day 2:N
  for (i in 2:N) {
    x_sample <- extraDistr::rtpois(n_particles, 1 + noisy_prevalence[i + 1, 2] / proportion_obs, a = 0)
    log_weights <- extraDistr::dskellam(x_sample - x_resample, birth_rate, death_rate, log = T) +
                   dbinom(genetic_data[i + 1, 3], choose(genetic_data[i + 1, 2], 2), 2 * birth_rate / x_sample, log = T) +
                   dbinom(noisy_prevalence[i + 1, 2], x_sample, proportion_obs, log = T) -
                   extraDistr::dtpois(x_sample, 1 + noisy_prevalence[i + 1, 2] / proportion_obs, a = 0, log = T)
    lse_weights <- matrixStats::logSumExp(log_weights)
    int_llik <- int_llik + lse_weights - log_n_particles
    #normalise weights
    norm_weights <- exp(log_weights - lse_weights)
    if (sum(norm_weights, na.rm = T) != 1) {
      norm_weights <- rep(1/n_particles, n_particles)
    }
    if (n_particles > 1) {
      x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
    } else {
      x_resample <- x_sample
    }
  }

  return(int_llik)
}
