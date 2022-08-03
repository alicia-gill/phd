#no epidemic
pm_noepi_parallel <- function(iter, birth_rate_0, max_birth_rate=100, prevalence_0, death_rate, ptree, n_particles, n_cores) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(prevalence_0)
  stop_time <- n - 1

  output <- list("birth_rate" = birth_rate_0, "acceptance_rate" = 0, "run_time" = 0)
  n_accepted <- 0

  b_old <- birth_rate_0
  prev_old <- prevalence_0

  #prior is uniform
  f_hat_old <- genetic_llik(birth_rate = b_old, ptree = ptree, prevalence = prev_old, stop_time = stop_time)
  #proposal is skellam, so cancels

  for (i in 1:iter) {
    j <- 100*i/iter
    if (j %% 1 == 0) {
      print(paste0(j,"%"))
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
    f_hat_log <- unlist(parallel::mclapply(1:n_particles, propose_epi, birth_rate = b_new, death_rate = death_rate, stop_time = stop_time, ptree = ptree, mc.cores = n_cores))
    f_hat_new <- matrixStats::logSumExp(f_hat_log) - log(n_particles)

    #step 3: compute acceptance probability
    logr <- f_hat_new - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)

    #step 4: accept/reject
    u <- runif(1,0,1)
    if (u <= a) {
      n_accepted <- n_accepted + 1
      b_old <- b_new
      f_hat_old <- f_hat_new
    }
    output$birth_rate[i] <- b_old
  }

  output$acceptance_rate <- n_accepted/iter
  output$run_time <- as.numeric(Sys.time()) - sys_time
  return(output)
}

