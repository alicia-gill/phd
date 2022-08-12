#noisy epidemic
smc <- function(iter, birth_rate_0, max_birth_rate=100, prevalence_0, death_rate, ptree, noisy_prevalence, proportion_obs, n_particles, plot=F) {
  sys_time <- as.numeric(Sys.time())

  n <- nrow(noisy_prevalence)
  stop_time <- n - 1

  output <- list("birth_rate" = birth_rate_0, "acceptance_rate" = 0, "run_time" = 0)
  n_accepted <- 0

  b_old <- birth_rate_0
  prev_old <- prevalence_0

  if (plot == T) {
    plot(prevalence_0, type="l")
  }

  #prior is uniform
  f_hat_old <- skellam_llik(birth_rate = b_old, death_rate = death_rate, prevalence = prev_old) +
    genetic_llik(birth_rate = b_old, ptree = ptree, prevalence = prev_old, stop_time = stop_time) +
    sum(dbinom(x = noisy_prevalence[-1,2], size = prev_old[-1,2], prob = proportion_obs, log = TRUE))

  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time)

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
    if (plot == T) {
      f_hat_new <- sir(n_particles = n_particles, birth_rate = b_new, death_rate = death_rate, proportion_obs = proportion_obs, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, plot = T)
    } else {
      f_hat_new <- sir(n_particles = n_particles, birth_rate = b_new, death_rate = death_rate, proportion_obs = proportion_obs, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data)
    }

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

