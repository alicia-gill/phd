#' Sample Importance Resample - Adaptive
#'
#' Implements an adaptive SIR algorithm and returns an approximated log-likelihood. Prevalence proposals are negative binomial and birth rate proposals are linear Gaussian.
#'
#' @param n_particles number of particles used in the sampling.
#' @param sigma standard deviation of the linear gaussian proposal for birth rates.
#' @param death_rate death rate of the epidemic.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param genetic_data data frame of day, number of lineages and number of coalescences.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param backward_sim logical; if TRUE, uses backward simulation.
#'
#' @return log-likelihood
#' @export
#'
#' @examples
#' sir_mix(n_particles = 100, max_birth_rate = 1, death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
sir_mix <- function(n_particles, sigma, death_rate, noisy_prevalence, proportion_obs, genetic_data, ess_threshold = n_particles/2, resampling_scheme = "systematic", backward_sim = TRUE) {
  #N is number of days of the epidemic
  N <- nrow(noisy_prevalence) - 1
  #ancestor matrix
  anc <- matrix(nrow = N, ncol = n_particles)
  #birth rate matrix
  birth_rate <- matrix(nrow = N, ncol = n_particles)
  #prevalence matrix: (N+1) rows to include 1 case on day 0
  prevalence <- matrix(nrow = N + 1, ncol = n_particles)
  #particle weights
  normweights <- matrix(nrow = N, ncol = n_particles)
  logweights <- matrix(nrow = N, ncol = n_particles)
  #resample indicator
  resample <- rep(NA, N)
  #ess
  particles <- rep(NA, N)
  #zeroth day always has 1
  prevalence[1, ] <- 1

  #initialise smc likelihood approximation
  int_llik <- 0
  #initialise x_resample to be 1
  x_resample <- prevalence[1, ]
  #initial weights are equal
  logw <- rep(0, n_particles)

  for (i in 1:N) {
    #sample b
    if (i == 1) {
      rate <- 1/(death_rate*2)
      b_sample <- rexp(n_particles, rate=rate)
    } else {
      b_sample <- rnorm(n_particles, mean = b_resample, sd = sigma)
      #reflect off 0
      b_sample <- abs(b_sample)
    }

    #sample x
    if (noisy_prevalence[i+1,2] == 0) {
      x_sample <- x_resample + rtskellam(n = n_particles, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate)
      epi_llik <- 0
    } else {
      #q is the sample to propose according to the data
      #if proportion_obs is 0, we want to sample entirely from the prior
      #if proportion_obs is 1, we want to sample entirely from the data
      #from proportion_obs>=0.2, only sample from the data
      q <- min(0.9, proportion_obs/0.2)
      #if index is 0, sample from prior
      #if index is 1, sample from data
      index <- rep(0, n_particles)
      u <- runif(n_particles)
      index[u <= q] <- 1
      x_sample <- rep(NA, n_particles)
      sum_index <- sum(index)
      if (sum_index == 0) {
        x_sample <- x_resample + rtskellam(n = n_particles, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate)
      } else if (sum_index == n_particles) {
        x_sample <- noisy_prevalence[i+1,2] + rnbinom(n = n_particles, size = noisy_prevalence[i + 1, 2], p = proportion_obs)
      } else {
        x_sample[index==0] <- x_resample[index==0] + rtskellam(n = n_particles-sum_index, old_x = x_resample[index==0], birth_rate = b_sample[index==0], death_rate = death_rate)
        x_sample[index==1] <- noisy_prevalence[i+1,2] + rnbinom(n = sum_index, size = noisy_prevalence[i + 1, 2], p = proportion_obs)
      }
      #proposal probability is (1-q)*prior + q*data
      lse_vector <- matrix(data = c(log(1-q) + smc_skellam(new_x = x_sample, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate, log = T), log(q) + dnbinom(x = x_sample - noisy_prevalence[i + 1, 2], size = noisy_prevalence[i + 1, 2], prob = proportion_obs, log = T)), ncol=2)
      epi_llik <- smc_skellam(new_x = x_sample, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate, log = T) - apply(lse_vector, 1, matrixStats::logSumExp)
    }

    # if (noisy_prevalence[i+1,2] == 0) {
    #   x_sample <- x_resample + rtskellam(n = n_particles, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate)
    #   epi_llik <- 0
    # } else {
    #   x_sample <- noisy_prevalence[i + 1, 2] + rnbinom(n = n_particles, size = noisy_prevalence[i + 1, 2], p = proportion_obs)
    #   epi_llik <- smc_skellam(new_x = x_sample, old_x = x_resample, birth_rate = b_sample, death_rate = death_rate, log = T) - dnbinom(x = x_sample - noisy_prevalence[i + 1, 2], size = noisy_prevalence[i + 1, 2], prob = proportion_obs, log = T)
    # }

    prevalence[i + 1, ] <- x_sample
    birth_rate[i, ] <- b_sample

    #compute weights
    if (is.null(genetic_data)==1) {
      genetic_llik <- 0
    } else {
      genetic_llik <- dbinom(x = genetic_data[i + 1, 3], size = choose(genetic_data[i + 1, 2], 2), prob = 1 - exp( - 2 * b_sample / x_sample), log = T)
    }

    log_weights <- logw + epi_llik + genetic_llik + dbinom(x = noisy_prevalence[i + 1, 2], size = x_sample, prob = proportion_obs, log = T)

    #if all impossible, then mission abort
    if (max(log_weights) == -Inf) {
      int_llik <- -Inf
      #return(int_llik)
      return(list("int_llik"=int_llik, "birth_rate"=rep(NA,N), "prevalence"=data.frame("day"=0:N, "prev"=c(1, rep(NA, N)))))
    }

    #normalise weights
    lse_weights <- matrixStats::logSumExp(log_weights)
    mean_weights <- lse_weights - matrixStats::logSumExp(logw)
    int_llik <- int_llik + mean_weights
    norm_weights <- exp(log_weights - lse_weights)
    logweights[i,] <- log_weights
    normweights[i,] <- norm_weights

    #resampling
    ess <- 1 / sum(norm_weights^2)
    particles[i] <- ess
    #if the ess is below threshold, then resample
    if (ess <= ess_threshold) {
      resample[i] <- 1
      logw <- rep(0, n_particles)
      if (resampling_scheme != "multinomial" & resampling_scheme != "systematic") {
        print("resampling scheme must be multinomial or systematic")
        break
      }
      if (resampling_scheme == "multinomial") {
        index <- sample(1:n_particles, n_particles, replace = T, prob = norm_weights)
      }
      if (resampling_scheme == "systematic") {
        index <- systematic_sample_cpp(n_particles, norm_weights)
      }
      anc[i,] <- index
      x_resample <- x_sample[index]
      b_resample <- b_sample[index]
    } else {
      resample[i] <- 0
      logw <- log_weights
      anc[i,] <- 1:n_particles
      x_resample <- x_sample
      b_resample <- b_sample
    }
  }

  b <- rep(NA, N)
  p <- data.frame("day"=0:N, "prev"=rep(1,N+1))
  if (backward_sim == TRUE) {
    jt <- rep(NA, N)
    #day N
    jt[N] <- sample(1:n_particles, 1, prob=normweights[N,])
    for (t in (N-1):1) {
      x <- birth_rate[t+1, jt[t+1]]
      y <- prevalence[t+2, jt[t+1]]
      lin_gau <- dnorm(x, mean = birth_rate[t,], sd = sigma, log = T)
      skel <- smc_skellam(new_x = y, old_x = prevalence[t+1,], birth_rate = x, death_rate = death_rate)
      denom <- matrixStats::logSumExp(logweights[t,] + lin_gau + skel)
      wtT <- logweights[t,] + lin_gau + skel - denom
      jt[t] <- sample(1:n_particles, 1, prob = exp(wtT))
    }
    for (i in 1:N) {
      b[i] <- birth_rate[i,jt[i]]
      p[i+1,2] <- prevalence[i+1, jt[i]]
    }
  } else {
    j <- sample(1:n_particles, 1, prob=normweights[N,])
    a <- anc[,j]
    for (i in 1:N) {
      b[i] <- birth_rate[i,a[i]]
      p[i+1,2] <- prevalence[i+1, a[i]]
    }
  }

  return(list("int_llik"=int_llik, "birth_rate"=b, "prevalence"=p, "ancestor"=anc, "weights"=normweights, "br_matrix"=birth_rate, "prev_matrix"=prevalence))
}
