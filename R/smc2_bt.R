#' Adaptive SMC-MCMC with Partial Epidemic Data with Varying Birth Rates
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param max_time maximum number of seconds to run the algorithm for.
#' @param max_birth_rate0 initial value of the maximum birth rate.
#' @param sigma0 intital value of the linear gaussian standard deviation.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, prevalence, acceptance rate, run time in seconds
#' @export
#'
#' @examples
#' smc2_bt(iter = 100000, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, proportion_obs = 0.2, n_particles = 100)
smc2_bt <- function(iter, max_time=Inf, max_birth_rate0, sigma0, death_rate, ptree, noisy_prevalence, proportion_obs, n_particles, ess_threshold = n_particles/2, print=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  p_matrix <- matrix(NA, nrow=iter, ncol=stop_time+1)
  max_b <- rep(NA, iter)
  sigma <- rep(NA, iter)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

  #prior on max_birth_rate is uniform
  #prior on sigma0 is uniform

  particles <- matrix(NA, nrow=iter, ncol=stop_time)

  max_b_old <- max_birth_rate0
  sigma_old <- sigma0

  sir <- sir_be(n_particles = n_particles, max_birth_rate = max_b_old, sigma = sigma_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold)
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

    #step 1: sample max_b_new and sample sigma_new
    #if proposal is negative, then bounce back
    max_b_new <- rnorm(1, max_b_old, 0.1)
    max_b_new <- abs(max_b_new)
    sigma_new <- rnorm(1, sigma_old, 0.1)
    sigma_new <- abs(sigma_new)

    #step 2: compute likelihood
    sir <- sir_be(n_particles = n_particles, max_birth_rate = max_b_new, sigma = sigma_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold)
    f_hat_new <- sir$int_llik
    b_new <- sir$birth_rate
    p_new <- sir$prevalence[,2]
    smc_llik[i] <- f_hat_new

    #step 3: compute acceptance probability
    logr <- f_hat_new - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)

    #step 4: accept/reject
    u <- runif(1,0,1)
    if (u <= a) {
      n_accepted <- n_accepted + 1
      b_old <- b_new
      p_old <- p_new
      f_hat_old <- f_hat_new
      max_b_old <- max_b_new
      sigma_old <- sigma_new
    }
    b_matrix[i,] <- b_old
    p_matrix[i,] <- p_old
    max_b[i] <- max_b_old
    sigma[i] <- sigma_old

    i <- i + 1
    run_time <- as.numeric(Sys.time()) - sys_time
  }

  output <- list("birth_rate" = b_matrix, "prevalence" = p_matrix, "max_birth_rate" = max_b, "sigma" = sigma, "acceptance_rate" = n_accepted/(i-1), "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik"=smc_llik)
  return(output)
}
