sir <- function(n_particles, birth_rate, death_rate, proportion_obs, noisy_prevalence, genetic_data, plot=F) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  samples <- matrix(nrow = N + 1, ncol = n_particles)
  #first day always has 1
  samples[1, ] <- 1

  #initialise smc likelihood approximation
  int_llik <- 0

  # bias <- matrix(nrow=N, ncol=3)

  #Day 1
  #sample from a truncated poisson
  lambda <- max(1, noisy_prevalence[2,2] / proportion_obs)
  a <- max(0, noisy_prevalence[2, 2] - 1)
  x_sample <- extraDistr::rtpois(n_particles, lambda = lambda, a = a)
  log_weights <- extraDistr::dskellam(x_sample - 1, birth_rate, death_rate, log = T) +
                 dbinom(genetic_data[2, 3], choose(genetic_data[2, 2], 2), 2 * birth_rate / x_sample, log = T) +
                 dbinom(noisy_prevalence[2, 2], x_sample, proportion_obs, log = T) -
                 extraDistr::dtpois(x_sample, lambda = lambda, a = a, log = T)

  #if all impossible, then mission abort
  if (max(log_weights) == -Inf) {
    int_llik <- -Inf
    break
  }

  lse_weights <- matrixStats::logSumExp(log_weights)
  log_n_particles <- log(n_particles)
  int_llik <- int_llik + lse_weights - log_n_particles
  #normalise weights
  norm_weights <- exp(log_weights - lse_weights)
  if (n_particles > 1) {
    x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
  } else {
    x_resample <- x_sample
  }
  samples[2, ] <- x_resample

  # bias[1,1] <- mean(x_sample)
  # bias[1,2] <- mean(x_resample)
  # bias[1,3] <- bias[1,1] - bias[1,2]

  #Day 2:N
  for (i in 2:N) {
    lambda <- max(1, noisy_prevalence[i + 1, 2] / proportion_obs)
    a <- max(0, noisy_prevalence[i + 1, 2] - 1)
    x_sample <- extraDistr::rtpois(n_particles, lambda, a = a)
    log_weights <- extraDistr::dskellam(x_sample - x_resample, birth_rate, death_rate, log = T) +
                   dbinom(genetic_data[i + 1, 3], choose(genetic_data[i + 1, 2], 2), 2 * birth_rate / x_sample, log = T) +
                   dbinom(noisy_prevalence[i + 1, 2], x_sample, proportion_obs, log = T) -
                   extraDistr::dtpois(x_sample, lambda, a = a, log = T)

    #if all impossible, then mission abort
    if (max(log_weights) == -Inf) {
      int_llik <- -Inf
      break
    }

    lse_weights <- matrixStats::logSumExp(log_weights)
    int_llik <- int_llik + lse_weights - log_n_particles
    #normalise weights
    norm_weights <- exp(log_weights - lse_weights)
    if (n_particles > 1) {
      x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
    } else {
      x_resample <- x_sample
    }
    samples[i + 1, ] <- x_resample
    # bias[i, 1] <- mean(x_sample)
    # bias[i, 2] <- mean(x_resample)
    # bias[i, 3] <- bias[i, 1] - bias[i, 2]
  }

  if (plot == T) {
    j <- sample(1:n_particles, 1, prob=norm_weights)
    sample <- samples[,j]
    lines(sample)
  }

  return(int_llik)
}
