#' PMMH with Unknown Number of Particles
#'
#' Implements an SMC-MCMC with adaptive resampling when the epidemic with birth rate varying over time has been partially observed with known proportion. MH adaptation of shape and scale. Unknown optimum number of particles.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param max_time maximum number of seconds to run the algorithm for.
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs proportion of cases observed.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param day number of days in the past the most recent leaf was sampled.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param backward_sim logical; if TRUE, uses backward simulation.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, prevalence, proportion observed, maximum birth rate, linear gaussian variance, acceptance rate, run time in seconds, smc log-likelihood and optimal number of particles
#' @export
#'
#' @examples
#' smc_nopt_knownp(iter = 100000, sigma0 = 0.1, proportion_obs = 0.5, death_rate = 0.1, ptree = sample_tree, day = 0, noisy_prevalence = noisy_prev, print = T)
smc_nopt_knownp <- function(iter, max_time=Inf, sigma0, proportion_obs, death_rate, ptree, day = 0, noisy_prevalence, n_particles=NULL, ess_threshold = n_particles/2, resampling_scheme = "systematic", backward_sim = TRUE, print=F) {
  sys_time <- as.numeric(Sys.time())

  if (is.null(n_particles)) {
    n_particles <- rep(NA, 3)
    for (i in 1:3) {
      n_particles[i] <- find_nopt_knownp(sigma0 = sigma0, proportion_obs = proportion_obs, death_rate = death_rate, ptree = ptree, day = day, noisy_prevalence = noisy_prevalence, resampling_scheme = resampling_scheme)
    }
    n_particles <- min(3*max(n_particles), 5000)
    ess_threshold <- n_particles/2
  }

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)

  b_matrix <- matrix(NA, nrow=iter, ncol=stop_time)
  p_matrix <- matrix(NA, nrow=iter, ncol=stop_time+1)
  sigma <- rep(NA, iter)
  n_accepted <- 0
  smc_llik <- rep(NA, iter)
  particles <- matrix(NA, nrow=iter, ncol=stop_time)
  scale <- rep(NA, iter+1)
  acceptance <- rep(0, iter)

  sigma_old <- sigma0

  s <- 0
  zeta <- 10000 #can be any constant >= 1
  inv_zeta <- 1/zeta

  scale[1] <- s

  lambda_sigma <- 1/0.1 #mean of 0.1

  mu <- matrix(NA, nrow=iter+1, ncol=1)
  Sigma <- array(NA, c(iter+1,1,1))
  mu_old <- sigma_old
  Sigma_old <- 1
  sqrtSigma_old <- sqrt(Sigma_old)
  mu[1,] <- mu_old
  Sigma[1,,] <- Sigma_old

  #prior on sigma is exponential
  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T)
  sir <- sir_mix(n_particles = n_particles, sigma = sigma_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
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
    w <- rnorm(n = 1, mean = 0, sd = 1)
    sigma_new <- sigma_old + exp(s) * sqrtSigma_old * w
    sigma_new <- abs(sigma_new)

    #step 2: compute likelihood
    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T)
    sir <- sir_mix(n_particles = n_particles, sigma = sigma_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = backward_sim)
    f_hat_new <- sir$int_llik
    b_new <- sir$birth_rate
    p_new <- sir$prevalence[,2]

    #step 3: compute acceptance probability
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)

    eta <- (i + 10)^(-0.6)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta

    scale[i+1] <- s

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
    }
    b_matrix[i,] <- b_old
    p_matrix[i,] <- p_old
    sigma[i] <- sigma_old
    smc_llik[i] <- f_hat_old

    Xn <- sigma_old
    mu_new <- (1 - eta) * mu_old + eta * Xn
    Sigma_new <- (1 - eta) * Sigma_old + eta * (Xn - mu_old) %*% t(Xn - mu_old)
    evalues <- eigen(Sigma_new, symmetric = T)$values
    min_evalue <- min(evalues)*exp(s)
    max_evalue <- max(evalues)*exp(s)
    if (norm(mu_new, type="2") <= zeta & max_evalue <= zeta & min_evalue >= inv_zeta) {
      mu_old <- mu_new
      Sigma_old <- Sigma_new
      sqrtSigma_old <- expm::sqrtm(Sigma_old)
    }
    mu[i+1,] <- mu_old
    Sigma[i+1,,] <- Sigma_old

    i <- i + 1
    run_time <- as.numeric(Sys.time()) - sys_time
  }

  output <- list("birth_rate" = b_matrix, "prevalence" = p_matrix, "sigma" = sigma, "acceptance_rate" = n_accepted/(i-1), "run_time" = as.numeric(Sys.time()) - sys_time, "smc_llik" = smc_llik, "s"=scale, "accept"=acceptance, "mu"=mu, "Sigma"=Sigma, "n_particles" = n_particles)
  return(output)
}
