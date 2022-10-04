#' SMC-MCMC with Partial Epidemic Data
#'
#' Implements an SMC-MCMC when the epidemic has been partially observed with unknown proportion.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param birth_rate_0 initial value of the birth rate.
#' @param max_birth_rate maximum value for the birth rate. 100 by default.
#' @param proportion_obs_0 initial value of the proportion of cases observed.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param n_particles number of particles used in the importance sampling.
#' @param ess_threshold threshold of ESS below which triggers resampling.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#' @param plot logical; if TRUE, then plots a trajectory according to weight.
#'
#' @return list containing: birth rate, proportion observed, acceptance rate, run time in seconds
#' @export
#'
#' @examples
#' smc_adaptive_bp(iter = 100000, birth_rate_0 = 0.1, proportion_obs_0 = 0.2, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, n_particles = 100)
smc_skellam_bp <- function(iter, birth_rate_0, max_birth_rate=100, proportion_obs_0, death_rate, ptree, noisy_prevalence, n_particles, ess_threshold = n_particles/2, print=F, plot=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  b_old <- birth_rate_0
  p_old <- proportion_obs_0

  if (plot == T) {
    plot(NA, xlim=c(0, stop_time), ylim=c(0,2000), xlab="Day", ylab="Prevalence")
  }

  output <- list("birth_rate" = b_old, "proportion_obs" = p_old, "smc_llik" = 0, "acceptance_rate" = 0, "run_time" = 0)
  n_accepted <- 0

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

  #adaptive mh
  s <- 0
  X <- c(b_old, p_old)
  mu_old <- X
  Sigma_old <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

  mu <- c(0, 0)
  Sigma <- matrix(c(0.1, 0.03, 0.03, 0.1), nrow = 2, ncol = 2)

  #prior is uniform
  f_hat_old <- sir_skellam(n_particles = n_particles, birth_rate = b_old, death_rate = death_rate, proportion_obs = p_old, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, ess_threshold = ess_threshold)

  for (i in 1:iter) {
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1: sample b_new, p_new and sample prev_new
    #adaptive mh
    w <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
    sqrtSigma <- expm::sqrtm(Sigma_old)
    Y <- X + exp(s) * sqrtSigma %*% w

    b_new <- Y[1]
    p_new <- Y[2]

    #if proposal is negative or larger than max_birth_rate, then bounce back
    while (b_new < 0 | b_new > max_birth_rate) {
      if (b_new < 0) {
        b_new <- -b_new
      }
      if (b_new > max_birth_rate) {
        b_new <- max_birth_rate - b_new
      }
    }
    while (p_new < 0 | p_new > 1) {
      if (p_new < 0) {
        p_new <- -p_new
      }
      if (p_new > 1) {
        p_new <- 1 - p_new
      }
    }

    #step 2: compute likelihood
    f_hat_new <- sir_skellam(n_particles = n_particles, birth_rate = b_new, death_rate = death_rate, proportion_obs = p_new, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, ess_threshold = ess_threshold, plot = plot)
    output$smc_llik[i] <- f_hat_new

    #step 3: compute acceptance probability
    logr <- f_hat_new - f_hat_old
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
      p_old <- p_new
      f_hat_old <- f_hat_new
    }
    output$birth_rate[i] <- b_old
    output$proportion_obs[i] <- p_old

    #adaptive mh
    X <- c(b_old, p_old)
    mu_new <- (1 - eta) * mu_old + eta * X
    Sigma_new <- (1 - eta) * Sigma_old + eta * (X - mu_old) %*% t(X - mu_old)

    norm <- norm(mu_new, type="2")
    evalues <- eigen(Sigma_new)$values
    if (norm <= 1000 && evalues >= 1/1000 && evalues <= 1000) {
      mu_old <- mu_new
      Sigma_old <- Sigma_new
    }
  }

  output$acceptance_rate <- n_accepted/iter
  output$run_time <- as.numeric(Sys.time()) - sys_time
  return(output)
}

