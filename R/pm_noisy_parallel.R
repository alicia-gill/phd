#' Pseudo-Marginal MCMC with Partial Epidemic Data - Parallelised
#'
#' Parallelised version of pm_noisy().
#'
#' @param iter number of iterations to run the algorithm for.
#' @param birth_rate_0 initial value of the birth rate.
#' @param max_birth_rate maximum value for the birth rate. 100 by default.
#' @param prevalence_0 data frame of initial values for the prevalence per day.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param proportion_obs proportion of cases observed.
#' @param n_particles number of particles used in the importance sampling.
#' @param n_cores number of cores to use in the parallelisation.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, acceptance rate, run time in seconds
#' @export
#'
#' @examples
#' prev_0 <- data.frame("day"=0:50, "prev"=c(1, floor(noisy_prev[-1,2] / proportion_obs)))
#' prev_0[prev_0[,2] == 0, 2] <- 1
#' pm_noisy_parallel(iter = 100000, birth_rate_0 = 0.1, prevalence_0 = prev_0, death_rate = 0.1, ptree = sample_tree, noisy_prevalence = noisy_prev, proportion_obs = 0.2, n_particles = 100, n_cores = 8)
pm_noisy_parallel <- function(iter, birth_rate_0, max_birth_rate=100, prevalence_0, death_rate, ptree, noisy_prevalence, proportion_obs, n_particles, n_cores, print=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  output <- list("birth_rate" = birth_rate_0, "acceptance_rate" = 0, "run_time" = 0)
  n_accepted <- 0

  b_old <- birth_rate_0
  prev_old <- prevalence_0

  #prior is uniform
  f_hat_old <- skellam_llik(birth_rate = b_old, death_rate = death_rate, prevalence = prev_old) +
    genetic_llik(birth_rate = b_old, ptree = ptree, prevalence = prev_old, stop_time = stop_time) +
    sum(dbinom(x = noisy_prevalence[-1,2], size = prev_old[-1,2], prob = proportion_obs, log = TRUE))

  smooth_prev <- smooth(noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs)
  lambda <- pmax(1, smooth_prev[-1,2])
  a <- noisy_prevalence[-1,2] + round(lambda) - smooth_prev[-1,2]
  q_old <- sum(extraDistr::dtpois(prev_old[-1,2] - smooth_prev[-1,2] + round(lambda), lambda = lambda, a = a, log = TRUE))

  llik_old <- f_hat_old - q_old

  for (i in 1:iter) {
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1: sample b_new and sample prev_new
    eps <- rnorm(1,0,0.01)
    b_new <- b_old + eps
    #if proposal is negative or larger than max_birth_rate, then bounce back
    while (b_new < 0 | b_new > max_birth_rate) {
      if (b_new < 0) {
        b_new <- -b_new
      }
      if (b_new > max_birth_rate) {
        b_new <- max_birth_rate - b_new
      }
    }

    #step 2: compute likelihood
    llik <- unlist(parallel::mclapply(1:n_particles, propose_pois2, birth_rate = b_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, ptree = ptree, mc.cores = n_cores))
    llik_new <- matrixStats::logSumExp(llik) - log(n_particles)

    #step 3: compute acceptance probability
    logr <- llik_new - llik_old
    loga <- min(0,logr)
    a <- exp(loga)

    #step 4: accept/reject
    u <- runif(1,0,1)
    if (u <= a) {
      n_accepted <- n_accepted + 1
      b_old <- b_new
      llik_old <- llik_new
    }
    output$birth_rate[i] <- b_old
  }

  output$acceptance_rate <- n_accepted/iter
  output$run_time <- as.numeric(Sys.time()) - sys_time
  return(output)
}

