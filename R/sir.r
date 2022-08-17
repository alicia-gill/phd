sir <- function(n_particles, birth_rate, death_rate, proportion_obs, noisy_prevalence, genetic_data, plot=F) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  samples <- matrix(nrow = N + 1, ncol = n_particles)
  #first day always has 1
  samples[1, ] <- 1

  #initialise smc likelihood approximation
  int_llik <- 0
  #initialise x_resample to be 1
  x_resample <- samples[1,]

  # bias <- matrix(nrow=N, ncol=3)

  log_n_particles <- log(n_particles)

  for (i in 1:N) {
    #sample from poisson centred at noisy_prev/prop_obs, or 1 if noisy_prev is 0
    lambda <- max(1, noisy_prevalence[i + 1, 2] / proportion_obs)
    #truncate at maximum of 1 and noisy_prev
    a <- max(0, noisy_prevalence[i + 1, 2] - 1)
    #sample
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
    #get int_llik by adding log of the mean weight
    int_llik <- int_llik + lse_weights - log_n_particles
    #normalise weights
    norm_weights <- exp(log_weights - lse_weights)
    #if n_particles==1, then sample() is weird
    if (n_particles > 1) {
      x_resample <- sample(x_sample, n_particles, replace = T, prob = norm_weights)
    } else {
      x_resample <- x_sample
    }
    samples[i + 1, ] <- x_resample
  }

  #if plot==T, sample one trajectory according to final
  if (plot == T) {
    j <- sample(1:n_particles, 1, prob=norm_weights)
    sample <- samples[,j]
    lines(sample)
  }

  return(list("int_llik"=int_llik, "samples"=samples))
}
