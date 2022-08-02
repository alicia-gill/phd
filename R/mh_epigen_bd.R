mh_epigen_bd <- function(iter, birth_rate_0, max_birth_rate=100, death_rate_0, max_death_rate=100, prevalence, ptree) {
  sys_time <- as.numeric(Sys.time()) #time at beginning of run

  #initialise output
  output <- list("birth_rate"=0, "death_rate"=0, "acceptance_rate"=0, "run_time"=0)

  #number of accepted moves
  n_accepted <- 0

  #initialise s (adaptive MCMC)
  s <- 0

  stop_time <- (nrow(prevalence) - 1)

  b_old <- birth_rate_0
  d_old <- death_rate_0
  #uniform prior, so always the same
  epi_llik_old <- skellam_llik(birth_rate = b_old, death_rate = d_old, prevalence = prevalence)
  gen_llik_old <- genetic_llik(birth_rate = b_old, ptree = ptree, prevalence = prevalence, stop_time = stop_time)

  for (i in 1:iter) {
    #prints how far through the chain we are
    j <- 100*i/iter
    if (j %% 1 == 0) {
      print(paste0(j,"%"))
    }

    #step 1 - generate new value of beta (adaptive MCMC)
    w1 <- rnorm(1,0,1)
    b_new <- b_old + exp(s)*w1
    #if proposal is negative, then bounce back
    while (b_new < 0 | b_new > max_birth_rate) {
      if (b_new < 0) {
        b_new <- -b_new
      }
      if (b_new > max_birth_rate) {
        b_new <- max_birth_rate - b_new
      }
    }

    w2 <- rnorm(1,0,1)
    d_new <- d_old + exp(s)*w2
    while (d_new < 0 | d_new > max_death_rate) {
      if (d_new < 0) {
        d_new <- -d_new
      }
      if (d_new > max_death_rate) {
        d_new <- max_death_rate - d_new
      }
    }

    #step 2 - calculate acceptance probability
    #prior for beta is uniform
    #r is ratio of new posterior to old posterior
    epi_llik_new <- skellam_llik(birth_rate = b_new, death_rate = d_new, prevalence = prev)
    gen_llik_new <- genetic_llik(birth_rate = b_new, ptree = ptree, prevalence = prev, stop_time = stop_time)

    logr <- epi_llik_new + gen_llik_new - epi_llik_old - gen_llik_old

    #acceptance probability
    loga <- min(0,logr)
    a <- exp(loga)
    s <- s + (a - 0.234) / i #update s (adaptive MCMC)

    #step 3 - accept/reject update
    u <- runif(1,0,1)
    #accept
    if (u <= a) {
      b_old <- b_new
      d_old <- d_new
      epi_llik_old <- epi_llik_new
      gen_llik_old <- gen_llik_new
      n_accepted <- n_accepted + 1
    }
    output$birth_rate[i] <- b_old
    output$death_rate[i] <- d_old
  }

  output$acceptance_rate <- n_accepted/iter
  output$run_time <- as.numeric(Sys.time()) - sys_time

  return(output)
}
