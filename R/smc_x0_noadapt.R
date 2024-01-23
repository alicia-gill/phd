#' PMMH with Unknown Number of Particles
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion. MH adaptation of shape and scale of rho and sigma. Unknown optimum number of particles.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param max_time maximum number of seconds to run the algorithm for.
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs0 initial value of the proportion of cases observed.
#' @param x0 initial value of the prevalence on day 0.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param day number of days in the past the most recent leaf was sampled.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param max_n_particles upper bound on the number of particles to use in the importance sampling.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param backward_sim logical; if TRUE, uses backward simulation.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, prevalence, proportion observed, maximum birth rate, linear gaussian variance, acceptance rate, run time in seconds, smc log-likelihood and optimal number of particles
#' @export
#'
#' @examples
#' smc_x0_noadapt(iter = 100000, sigma0 = 0.1, x0 = 1, proportion_obs0 = 0.5, death_rate = 0.1, ptree = sample_tree, day = 0, noisy_prevalence = noisy_prev, print = T)
smc_x0_noadapt <- function(iter, max_time = Inf, sigma0, proportion_obs0, x0 = 1, death_rate, ptree, day = 0, genetic_data = NULL, noisy_prevalence, n_particles = NULL, ess_threshold = n_particles/2, max_n_particles = 10000, resampling_scheme = "systematic", backward_sim = TRUE, print = F) {
  sys_time <- as.numeric(Sys.time())

  if (is.null(n_particles)) {
    n_particles <- rep(NA, 3)
    for (i in 1:3) {
      n_particles[i] <- find_nopt_noadapt(sigma0 = sigma0, proportion_obs0 = proportion_obs0, death_rate = death_rate, ptree = ptree, day = day, noisy_prevalence = noisy_prevalence, resampling_scheme = resampling_scheme)
    }
    n_particles <- min(max(n_particles), max_n_particles)
    ess_threshold <- n_particles/2
  }

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  if (is.null(genetic_data)) {
    genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)
  }

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  p_matrix <- matrix(NA, nrow=iter, ncol=stop_time+1)
  sigma <- rep(NA, iter)
  x0_output <- rep(NA, iter)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)
  particles <- matrix(NA, nrow=iter, ncol=stop_time)
  acceptance <- rep(0, iter)

  sigma_old <- sigma0
  x0_old <- x0

  sum_noisy <- sum(noisy_prevalence[-1,2])
  if (sum_noisy == 0) {
    p_obs <- rep(0, iter)
    p_obs_old <- 0
  } else {
    p_obs <- rep(NA, iter)
    p_obs_old <- proportion_obs0
  }

  lambda_sigma <- 1/0.1 #mean of 0.1
  lambda_x0 <- 1 #mean and var of 1
  alpha_pobs <- 1
  beta_pobs <- 3

  #prior on sigma is exponential
  #prior on x0-1 is uniform(1,Inf)
  #prior on p_obs is uniform(0,1) or beta(1,3)
  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dunif(x = p_obs_old, min = 0, max = 1, log = T)
  #  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dbeta(x = p_obs_old, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
  sir <- sir_mix(n_particles = n_particles, sigma = sigma_old, x0 = x0_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_old, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
  f_hat_old <- sir$int_llik
  b_old <- sir$birth_rate
  p_old <- sir$prevalence[,2]

  i <- 1
  run_time <- as.numeric(Sys.time()) - sys_time
  while (i <= iter & run_time <= max_time) {
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1: sample sample sigma_new
    #if proposal is negative, then bounce back
    sigma_new <- sigma_old + rnorm(n=1, mean=0, sd=0.01)
    sigma_new <- abs(sigma_new)
    x0_new <- x0_old + extraDistr::rsign(1)
    if (x0_new <= 0) {
      x0_new <- x0_old + 1
    }
    if (sum_noisy == 0) {
      p_obs_new <- 0
    } else {
      p_obs_new <- p_obs_old + rnorm(n=1, mean=0, sd=0.01)
      while (p_obs_new < 0 | p_obs_new > 1) {
        if (p_obs_new > 1) {
          p_obs_new <- 1 - p_obs_new
        }
        if (p_obs_new < 0) {
          p_obs_new <- -p_obs_new
        }
      }
    }

    #step 2: compute likelihood
    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dunif(x = p_obs_new, min = 0, max = 1, log = T)
    #    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dbeta(x = p_obs_new, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
    sir <- sir_mix(n_particles = n_particles, sigma = sigma_new, x0 = x0_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_new, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
    f_hat_new <- sir$int_llik
    b_new <- sir$birth_rate
    p_new <- sir$prevalence[,2]

    #step 3: compute acceptance probability
    if (x0_old > 1) {
      logq <- 0
    } else {
      logq <- -log(2)
    }
    logr <- prior_new + f_hat_new + logq - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)

    #step 4: accept/reject
    # -log(U) ~ Exp(1) where U ~ Uniform(0,1)
    logu <- -rexp(1)
    if (logu <= loga) {
      acceptance[i] <- 1
      n_accepted <- n_accepted + 1
      b_old <- b_new
      p_old <- p_new
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      x0_old <- x0_new
      p_obs_old <- p_obs_new
    }
    b_matrix[i,] <- b_old
    p_matrix[i,] <- p_old
    sigma[i] <- sigma_old
    x0_output[i] <- x0_old
    smc_llik[i] <- f_hat_old

    if (sum_noisy > 0) {
      p_obs[i] <- p_obs_old
    }

    i <- i + 1
    run_time <- as.numeric(Sys.time()) - sys_time
  }

  output <- list("birth_rate" = b_matrix, "prevalence" = p_matrix, "sigma" = sigma, "x0" = x0_output, "proportion_obs" = p_obs, "acceptance_rate" = n_accepted/(i-1), "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik" = smc_llik, "accept"=acceptance, "n_particles" = n_particles)
  return(output)
}
