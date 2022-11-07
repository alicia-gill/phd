#' Sample Importance Resample - Adaptive
#'
#' Implements an adaptive SIR algorithm and returns an approximated log-likelihood. Proposals are Skellam.
#'
#' @param n_particles number of particles used in the sampling.
#' @param birth_rate_max maximum value of the birth rate on day 1 of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param genetic_data data frame of day, number of lineages and number of coalescences.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param plot logical; if TRUE, then plots a prevalence trajectory according to weight.
#'
#' @return log-likelihood
#' @export
#'
#' @examples
#' sir_be(n_particles = 100, death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
sir_be <- function(n_particles, max_birth_rate = 10, death_rate, noisy_prevalence, proportion_obs, genetic_data, ess_threshold = n_particles/2, plot=F) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  anc <- matrix(nrow = N, ncol = n_particles)
  birth_rate <- matrix(nrow = N, ncol = n_particles)
  prevalence <- matrix(nrow = N + 1, ncol = n_particles)
  weights <- matrix(nrow = N, ncol = n_particles)
  resample <- rep(NA, N)
  particles <- rep(NA, N)
  #first day always has 1
  prevalence[1, ] <- 1

  #initialise smc likelihood approximation
  int_llik <- 0
  #initialise x_resample to be 1
  x_resample <- prevalence[1, ]
  w <- 1/n_particles

  for (i in 1:N) {
    #sample b
    if (i == 1) {
      b_sample <- runif(n_particles, min=0, max=max_birth_rate)
    } else {
      b_sample <- rnorm(n_particles, mean = b_resample, sd = 0.1)
      #reflect off 0
      b_sample <- abs(b_sample)
    }

    #sample x
    a <- pmax(1, noisy_prevalence[i+1,2], 2*b_sample)
    x_sample <- rep(NA, n_particles)

    set <- which(is.na(x_sample))
    check <- length(set)
    count <- 0
    while (check > 0) {
      skel_samp <- x_resample[set] + extraDistr::rskellam(check, b_sample[set]*x_resample[set], death_rate*x_resample[set])
      for (k in 1:check) {
        if (skel_samp[k] >= a[set[k]]) {
          x_sample[set[k]] <- skel_samp[k]
        }
      }
      set <- which(is.na(x_sample))
      check <- length(set)
      count <- count + 1
      if (count > 10000000) {
        int_llik <- -Inf
        #return(int_llik)
        return(list("int_llik"=int_llik, "resample_count"=resample, "particles"=particles, "birth_rate"=birth_rate, "prevalence"=prevalence, "weights"=weights))
      }
    }

    prevalence[i + 1, ] <- x_sample
    birth_rate[i, ] <- b_sample

    #compute weights
    log_weights <- dbinom(genetic_data[i + 1, 3], choose(genetic_data[i + 1, 2], 2), 2 * b_sample / x_sample, log = T) +
      dbinom(noisy_prevalence[i + 1, 2], x_sample, proportion_obs, log = T)

    log_weights <- ifelse(is.nan(log_weights), -Inf, log_weights)

    #if all impossible, then mission abort
    if (max(log_weights) == -Inf) {
      int_llik <- -Inf
      #return(int_llik)
      return(list("int_llik"=int_llik, "resample_count"=resample, "particles"=particles, "birth_rate"=birth_rate, "prevalence"=prevalence, "weights"=weights))
    }

    #normalise weights
    lse_weights <- matrixStats::logSumExp(log_weights)
    mean_weights <- matrixStats::logSumExp(log_weights + log(w))
    int_llik <- int_llik + mean_weights
    norm_weights <- exp(log_weights - lse_weights)

    #resampling
    ess <- 1 / sum(norm_weights^2)
    particles[i] <- round(ess)
    #if the ess is below threshold, then resample
    if (ess <= ess_threshold) {
      resample[i] <- 1
      w <- rep(1/n_particles, n_particles)
      #if n_particles==1, then sample() is weird
      index <- sample(1:n_particles, n_particles, replace = T, prob = norm_weights)
      anc[i,] <- index
      x_resample <- x_sample[index]
      b_resample <- b_sample[index]
    } else {
      resample[i] <- 0
      w <- norm_weights
      anc[i,] <- 1:n_particles
      x_resample <- x_sample
      b_resample <- b_sample
    }

    weights[i, ] <- w
  }

  #if plot==T, sample one trajectory according to final weight
  if (plot == T) {
    j <- sample(1:n_particles, 1, prob=norm_weights)
    prevalence <- prevalence[,j]
    lines(prevalence)
  }

  #return(int_llik)
  return(list("int_llik"=int_llik, "resample_count"=resample, "particles"=particles, "birth_rate"=birth_rate, "prevalence"=prevalence, "weights"=weights, "ancestors"=anc))
}
