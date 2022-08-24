#' Metropolis--Hastings using Epidemic and Genetic data
#'
#' Runs a Metropolis--Hastings algorithm to find the birth rate of the epidemic.
#'
#' @param iter number of iterations to run the algorithm for.
#' @param birth_rate_0 initial value of the birth rate.
#' @param max_birth_rate maximum value for the birth rate. 100 by default.
#' @param death_rate death rate of the epidemic.
#' @param prevalence data frame of prevalence per day.
#' @param ptree object of class phylo.
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return list containing: birth rate, acceptance rate, run time in seconds
#' @export
#'
#' @examples
#' mh_epigen(iter = 100000, birth_rate_0 = 0.1, death_rate = 0.1, prevalence = prev, ptree = sample_tree)
mh_epigen <- function(iter, birth_rate_0, max_birth_rate=100, death_rate, prevalence, ptree, print=F) {
  sys_time <- as.numeric(Sys.time()) #time at beginning of run

  #initialise output
  output <- list("birth_rate"=0, "acceptance_rate"=0, "run_time"=0)

  #number of accepted moves
  n_accepted <- 0

  #initialise s (adaptive MCMC)
  s <- 0

  stop_time <- (nrow(prevalence) - 1)

  b_old <- birth_rate_0
  #uniform prior, so always the same
  epi_llik_old <- skellam_llik(birth_rate = b_old, death_rate = death_rate, prevalence = prevalence)
  gen_llik_old <- genetic_llik(birth_rate = b_old, ptree = ptree, prevalence = prevalence, stop_time = stop_time)

  for (i in 1:iter) {
    #prints how far through the chain we are
    if (print == T) {
      j <- 100*i/iter
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #step 1 - generate new value of beta (adaptive MCMC)
    w <- rnorm(1,0,1)
    b_new <- b_old + exp(s)*w
    #if proposal is negative, then bounce back

    while (b_new < 0 | b_new > max_birth_rate) {
      if (b_new < 0) {
        b_new <- -b_new
      }
      if (b_new > max_birth_rate) {
        b_new <- max_birth_rate - b_new
      }
    }

    #step 2 - calculate acceptance probability
    #prior for beta is uniform
    #r is ratio of new posterior to old posterior
    epi_llik_new <- skellam_llik(birth_rate = b_new, death_rate = death_rate, prevalence = prev)
    gen_llik_new <- genetic_llik(birth_rate = b_new, ptree = ptree, prevalence = prev, stop_time = stop_time)

    logr <- epi_llik_new + gen_llik_new - epi_llik_old - gen_llik_old

    #acceptance probability
    loga <- min(0,logr)
    a <- exp(loga)
    s <- s + (a - 0.5) / i #update s (adaptive MCMC)

    #step 3 - accept/reject update
    u <- runif(1,0,1)
    #accept
    if (u <= a) {
      b_old <- b_new
      epi_llik_old <- epi_llik_new
      gen_llik_old <- gen_llik_new
      n_accepted <- n_accepted + 1
    }
    output$birth_rate[i] <- b_old
  }

  output$acceptance_rate <- n_accepted/iter
  output$run_time <- as.numeric(Sys.time()) - sys_time

  return(output)
}
