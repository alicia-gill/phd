#' Adaptive SMC-MCMC with Partial Epidemic Data with Varying Birth Rates
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param birth_rate_0 initial value of the vector of birth rates.
#' @param max_birth_rate maximum value for the birth rate. 100 by default.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#' @param plot logical; if TRUE, then plots a trajectory according to weight.
#'
#' @return list containing: birth rate, acceptance rate, run time in seconds
#' @export
#'
#' @examples
#' varying_bt(iter = 100000, birth_rate_0 = 0.1, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, proportion_obs = 0.2, n_particles = 100)
varying_bt <- function(iter, birth_rate_0, max_birth_rate = 100, death_rate, ptree, noisy_prevalence, proportion_obs, n_particles, ess_threshold = n_particles/2, print=F, plot=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)

  b_old <- birth_rate_0

  if (plot == T) {
    plot(NULL, xlim=c(0, stop_time), ylim=c(0,10000), xlab="Day", ylab="Prevalence")
  }

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

  #adaptive mh
  s <- 0

  mu <- rep(0, stop_time)
  Sigma <- diag(1, nrow = stop_time, ncol = stop_time)

  #prior on b1 is uniform, prior on bt|bt-1 is norm(bt-1, 0.01)
  prior_old <- sum(dnorm(b_old[2:stop_time], mean = b_old[1:(stop_time-1)], sd = 0.1, log = T))
  f_hat_old <- sir_adaptive_bt(n_particles = n_particles, birth_rate = b_old, death_rate = death_rate, proportion_obs = proportion_obs, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, ess_threshold = ess_threshold)

  for (i in 1:iter) {
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1: sample b_new and sample prev_new
    w <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
    sqrtSigma <- expm::sqrtm(Sigma)
    b_new <- b_old + exp(s) * sqrtSigma %*% w

    #if proposal is negative or larger than max_birth_rate, then bounce back
    for (k in 1:stop_time) {
      if (b_new[k] < 0) {
        b_new[k] <- -b_new[k]
      }
    }

    #step 2: compute likelihood
    prior_new <- sum(dnorm(b_new[2:stop_time], mean = b_new[1:(stop_time-1)], sd = 0.1, log = T))
    f_hat_new <- sir_adaptive_bt(n_particles = n_particles, birth_rate = b_new, death_rate = death_rate, proportion_obs = proportion_obs, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, ess_threshold = ess_threshold, plot = plot)
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
      b_old <- b_new
      prior_old <- prior_new
      f_hat_old <- f_hat_new
    }
    b_matrix[i,] <- b_old
  }

  output <- list("birth_rate" = b_matrix, "acceptance_rate" = n_accepted/iter, "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik"=smc_llik)
  return(output)
}
