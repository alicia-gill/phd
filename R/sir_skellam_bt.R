#' Sample Importance Resample - Adaptive
#'
#' Implements an adaptive SIR algorithm and returns an approximated log-likelihood. Proposals are Skellam.
#'
#' @param n_particles number of particles used in the sampling.
#' @param birth_rate vector of birth rates per day of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param genetic_data data frame of day, number of lineages and number of coalescences.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param plot logical; if TRUE, then plots a trajectory according to weight.
#'
#' @return log-likelihood, ancestral matrix, sampled prevalences, weights, ess
#' @export
#'
#' @examples
#' sir_skellam_bt(n_particles = 100, birth_rate = rep(0.2, 50), death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
sir_skellam_bt <- function(n_particles, birth_rate, death_rate, noisy_prevalence, proportion_obs, genetic_data, ess_threshold = n_particles/2, plot=F) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  ancestors <- matrix(nrow = N, ncol = n_particles)
  prevalence <- matrix(nrow = N + 1, ncol = n_particles)
  weights <- matrix(nrow = N + 1, ncol = n_particles)
  ess <- rep(NA, N)
  resample <- rep(NA, N)
  #first day always has 1
  prevalence[1, ] <- 1
  weights[1, ] <- 1/n_particles

  #initialise smc likelihood approximation
  int_llik <- 0
  #initialise x_resample to be 1
  x_resample <- prevalence[1, ]
  w <- weights[1, ]

  for (i in 1:N) {
    #birth rate
    bt <- birth_rate[i]
    #truncation
    a <- max(1, noisy_prevalence[i+1,2], 2*bt)
    x_sample <- rep(NA, n_particles)

    #which indexes haven't been successfully sampled yet
    set <- which(is.na(x_sample))
    #how many indexes is that
    check <- length(set)
    count <- 0
    while (check > 0) {
      #for unsampled, generate samples
      skel_samp <- x_resample[set] + extraDistr::rskellam(check, bt*x_resample[set], death_rate*x_resample[set])
      #if truncation applies, save them
      for (k in 1:check) {
        if (skel_samp[k] >= a) {
          x_sample[set[k]] <- skel_samp[k]
        }
      }
      set <- which(is.na(x_sample))
      check <- length(set)
      count <- count + 1
      #if this loop is run more than 100,000 times, then consider this birth rate trajectory impossible
      if (count > 100000) {
        skel_samp <- x_resample[set] + extraDistr::rskellam(check, bt*x_resample[set], death_rate*x_resample[set])
        for (k in 1:check) {
          if (skel_samp[k] >= 1) {
            x_sample[set[k]] <- skel_samp[k]
          }
        }
        set <- which(is.na(x_sample))
        check <- length(set)
      }
    }

    prevalence[i + 1, ] <- x_sample

    #compute weights
    if (is.null(genetic_data)==1) {
      genetic_llik <- 0
    } else {
      genetic_llik <- dbinom(genetic_data[i + 1, 3], choose(genetic_data[i + 1, 2], 2), 2 * bt / x_sample, log = T)
    }

    log_weights <- genetic_llik + dbinom(noisy_prevalence[i + 1, 2], x_sample, proportion_obs, log = T)

    log_weights <- ifelse(is.nan(log_weights), -Inf, log_weights)

    #if all impossible, then mission abort
    if (max(log_weights) == -Inf) {
      int_llik <- -Inf
      #return(int_llik)
      return(list("int_llik"=int_llik, "ancestors"=ancestors, "prevalence"=prevalence, "weights"=weights, "ess"=ess))
    }

    #normalise weights
    lse_weights <- matrixStats::logSumExp(log_weights)
    mean_weights <- matrixStats::logSumExp(log_weights + log(w))
    int_llik <- int_llik + mean_weights
    norm_weights <- exp(log_weights - lse_weights)

    #resampling
    ess_calc <- 1 / sum(norm_weights^2)
    ess[i] <- ess_calc
    #if the ess is below threshold, then resample
    if (ess_calc <= ess_threshold) {
      resample[i] <- 1
      w <- rep(1/n_particles, n_particles)
      index <- sample(1:n_particles, n_particles, replace = T, prob = norm_weights)
      ancestors[i, ] <- index
      x_resample <- x_sample[index]
    } else {
      resample[i] <- 0
      w <- norm_weights
      ancestors[i, ] <- 1:n_particles
      x_resample <- x_sample
    }

    weights[i + 1, ] <- w
  }

  #if plot==T, sample one trajectory according to final weight
  if (plot == T) {
    j <- sample(1:n_particles, 1, prob=w)
    prev_plot <- prevalence[,j]
    lines(prev_plot)
  }

  #return(int_llik)
  return(list("int_llik"=int_llik, "ancestors"=ancestors, "prevalence"=prevalence, "weights"=weights, "ess"=ess))
}
