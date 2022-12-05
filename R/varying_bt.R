#' Adaptive SMC-MCMC with Partial Epidemic Data with Varying Birth Rates
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param max_time maximum number of seconds to run the algorithm for.
#' @param birth_rate0 initial value of the vector of birth rates.
#' @param max_birth_rate0 initial value of the maximum birth rate.
#' @param sigma0 initial value of the linear Gaussian standard deviation.
#' @param proportion_obs0 initial value of the proportion of cases observed.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, proportion observed, maximum birth rate, linear gaussian variance, acceptance rate, run time in seconds and smc log-likelihood
#' @export
#'
#' @examples
#' varying_bt(iter = 100000, birth_rate_0 = rep(0.1, 50), max_birth_rate0 = 1, sigma0 = 0.1, proportion_obs0 = 0.1, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, n_particles = 100)
varying_bt <- function(iter, max_time=Inf, birth_rate0, max_birth_rate0, sigma0, proportion_obs0, death_rate, ptree, noisy_prevalence, n_particles, ess_threshold = n_particles/2, print=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  max_b <- rep(NA, iter)
  sigma <- rep(NA, iter)
  p_obs <- rep(NA, iter)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)
  particles <- matrix(NA, nrow=iter, ncol=stop_time)

  b_old <- birth_rate0
  max_b_old <- max_birth_rate0
  sigma_old <- sigma0
  p_obs_old <- proportion_obs0

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

  #adaptive mh
  s <- 0
  mu <- rep(0, stop_time+3)
  Sigma <- diag(1, nrow = stop_time+3, ncol = stop_time+3)
  sqrtSigma <- expm::sqrtm(Sigma)

  #prior on b1 is uniform, prior on bt|bt-1 is norm(bt-1, 0.01)
  prior_old <- sum(dnorm(b_old[2:stop_time], mean = b_old[1:(stop_time-1)], sd = sigma_old, log = T))
  sir <- sir_bin_bt(n_particles = n_particles, birth_rate = b_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_old, genetic_data = genetic_data, ess_threshold = ess_threshold)
  f_hat_old <- sir$int_llik

  i <- 1
  run_time <- as.numeric(Sys.time()) - sys_time
  while (i <= iter & run_time <= max_time) {
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1: sample b_new and sample prev_new
    w <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
    new <- c(b_old, max_b_old, sigma_old, p_obs_old) + exp(s) * sqrtSigma %*% w
    new <- abs(new)
    b_new <- new[1:stop_time]
    max_b_new <- new[stop_time+1]
    sigma_new <- new[stop_time+2]
    p_obs_new <- new[stop_time+3]
    while (p_obs_new < 0 | p_obs_new > 1) {
      if (p_obs_new > 1) {
        p_obs_new <- 1 - p_obs_new
      }
      if (p_obs_new < 0) {
        p_obs_new <- -p_obs_new
      }
    }

    #step 2: compute likelihood
    prior_new <- sum(dnorm(b_new[2:stop_time], mean = b_new[1:(stop_time-1)], sd = sigma_new, log = T))
    sir <- sir_bin_bt(n_particles = n_particles, birth_rate = b_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_new, genetic_data = genetic_data, ess_threshold = ess_threshold)
    f_hat_new <- sir$int_llik

    particles[i,] <- sir$ess
    smc_llik[i] <- f_hat_new

    #step 3: compute acceptance probability
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)

    #adaptive mh
    eta <- (i + 10)^(-0.6)
    s <- s + (a - 0.1) * eta

    #step 4: accept/reject
    u <- runif(1,0,1)
    if (u <= a) {
      n_accepted <- n_accepted + 1
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      b_old <- b_new
      max_b_old <- max_b_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }
    b_matrix[i,] <- b_old
    max_b[i] <- max_b_old
    sigma[i] <- sigma_old
    p_obs[i] <- p_obs_old

    i <- i + 1
    run_time <- as.numeric(Sys.time()) - sys_time
  }

  output <- list("birth_rate" = b_matrix, "max_birth_rate" = max_b, "sigma" = sigma, "proportion_obs" = p_obs, "acceptance_rate" = n_accepted/(i-1), "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik"=smc_llik, "particles"=particles)
  return(output)
}
