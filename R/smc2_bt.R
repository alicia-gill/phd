#' Adaptive SMC-MCMC with Partial Epidemic Data with Varying Birth Rates
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param max_time maximum number of seconds to run the algorithm for.
#' @param max_birth_rate0 initial value of the maximum birth rate.
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs0 initial value of the proportion of cases observed.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param backward_sim logical; if TRUE, uses backward simulation.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, prevalence, proportion observed, maximum birth rate, linear gaussian variance, acceptance rate, run time in seconds and smc log-likelihood
#' @export
#'
#' @examples
#' smc2_bt(iter = 100000, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, proportion_obs = 0.2, n_particles = 100)
smc2_bt <- function(iter, max_time=Inf, max_birth_rate0, sigma0, proportion_obs0, death_rate, ptree, noisy_prevalence, n_particles, ess_threshold = n_particles/2, resampling_scheme = "multinomial", backward_sim = TRUE, print=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  p_matrix <- matrix(NA, nrow=iter, ncol=stop_time+1)
  max_b <- rep(NA, iter)
  sigma <- rep(NA, iter)
  p_obs <- rep(NA, iter)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)

  scale <- rep(NA, iter+1)
  acceptance <- rep(0, iter)

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

  #prior on max_birth_rate is uniform
  #prior on sigma0 is uniform

  particles <- matrix(NA, nrow=iter, ncol=stop_time)

  max_b_old <- max_birth_rate0
  sigma_old <- sigma0
  p_obs_old <- proportion_obs0

  lambda_mbr <- 1 #mean of 1
  lambda_sigma <- 1/0.1 #mean of 0.1

  #prior on max_birth_rate and sigma are exponential
  #prior on p_obs is uniform(0,1)
  prior_old <- dexp(x = max_b_old, rate = lambda_mbr, log = T) + dexp(x = sigma_old, rate = lambda_sigma, log = T) + dunif(x = p_obs_old, min = 0, max = 1, log = T)
  sir <- sir_be(n_particles = n_particles, max_birth_rate = max_b_old, sigma = sigma_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_old, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
  f_hat_old <- sir$int_llik
  b_old <- sir$birth_rate
  p_old <- sir$prevalence[,2]

  s <- 0
  mu <- rep(0, 3)
  Sigma <- diag(1, nrow = 3, ncol = 3)
  sqrtSigma <- expm::sqrtm(Sigma)

  scale[1] <- s

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
    w <- MASS::mvrnorm(n = 1, mu = mu, Sigma = diag(1, nrow=3, ncol=3))
    new <- c(max_b_old, sigma_old, p_obs_old) + exp(s) * sqrtSigma %*% w
    new <- abs(new)
    max_b_new <- new[1]
    sigma_new <- new[2]
    p_obs_new <- new[3]
    while (p_obs_new < 0 | p_obs_new > 1) {
      if (p_obs_new > 1) {
        p_obs_new <- 1 - p_obs_new
      }
      if (p_obs_new < 0) {
        p_obs_new <- -p_obs_new
      }
    }
    # print(paste0("scale: ", s))
    # print(paste0("max_b_old: ", max_b_old))
    # print(paste0("max_b_new: ", max_b_new))
    # print(paste0("sigma_old: ", sigma_old))
    # print(paste0("sigma_new: ", sigma_new))
    # print(paste0("p_obs_old: ", p_obs_old))
    # print(paste0("p_obs_new: ", p_obs_new))

    #step 2: compute likelihood
    prior_new <- dexp(x = max_b_new, rate = lambda_mbr, log = T) + dexp(x = sigma_new, rate = lambda_sigma, log = T) + dunif(x = p_obs_new, min = 0, max = 1, log = T)
    sir <- sir_be(n_particles = n_particles, max_birth_rate = max_b_new, sigma = sigma_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_new, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
    f_hat_new <- sir$int_llik
    b_new <- sir$birth_rate
    p_new <- sir$prevalence[,2]
    smc_llik[i] <- f_hat_new
    # print(paste0("prior_old: ", prior_old))
    # print(paste0("prior_new: ", prior_new))
    # print(paste0("f_hat_old: ", f_hat_old))
    # print(paste0("f_hat_new: ", f_hat_new))

    #step 3: compute acceptance probability
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    # print(paste0("a: ", a))

    eta <- (i + 10)^(-0.6)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta

    scale[i+1] <- s

    #step 4: accept/reject
    # u <- runif(1,0,1)
    # -log(U) ~ Exp(1) where U ~ Uniform(0,1)
    logu <- -rexp(1)
    if (logu <= loga) {
      acceptance[i] <- 1
      n_accepted <- n_accepted + 1
      b_old <- b_new
      p_old <- p_new
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      max_b_old <- max_b_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }
    b_matrix[i,] <- b_old
    p_matrix[i,] <- p_old
    max_b[i] <- max_b_old
    sigma[i] <- sigma_old
    p_obs[i] <- p_obs_old

    i <- i + 1
    run_time <- as.numeric(Sys.time()) - sys_time
  }

  output <- list("birth_rate" = b_matrix, "prevalence" = p_matrix, "max_birth_rate" = max_b, "sigma" = sigma, "proportion_obs" = p_obs, "acceptance_rate" = n_accepted/(i-1), "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik"=smc_llik, "s"=scale, "accept"=acceptance)
  return(output)
}
