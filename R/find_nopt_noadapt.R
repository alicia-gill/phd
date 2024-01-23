#' Find optimal number of particles
#'
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs0 initial value of the proportion of cases observed.
#' @param x0 prevalence on day 0.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param day number of days in the past the most recent leaf was sampled.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return optimal number of particles
#' @export
#'
#' @examples
#' find_nopt_noadapt(sigma0 = 0.1, proportion_obs0 = 0.5, death_rate = 0.1, ptree = sample_tree, day = 1, noisy_prevalence = noisy_prev, print = T)
find_nopt_noadapt <- function(sigma0, proportion_obs0, x0 = 1, death_rate, ptree, day = 0, noisy_prevalence, resampling_scheme = "systematic", print = F) {
  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)

  sigma_old <- sigma0

  lambda_sigma <- 1/0.1 #mean of 0.1
  #alpha_pobs <- 1
  #beta_pobs <- 3

  sum_noisy <- sum(noisy_prevalence[-1,2])

  if (sum_noisy == 0) {
    p_obs_old <- 0
  } else {
    p_obs_old <- proportion_obs0
  }

  n_particles <- 10000
  ess_threshold <- n_particles/2

  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dunif(x = p_obs_old, min = 0, max = 1, log = T)
#  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dbeta(x = p_obs_old, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
  f_hat_old <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_old, proportion_obs = p_obs_old, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik

  #need to find a good estimate of the posterior mean for sigma and p_obs
  #first 500 for burn-in
  #second 500 to estimate means
  for (i in 1:500) {
    if (print == T) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }
    #if no data observed, then p_obs is always 0
    if (sum_noisy == 0) {
      sigma_new <- sigma_old + rnorm(n=1, mean=0, sd=0.01)
      sigma_new <- abs(sigma_new)
      p_obs_new <- 0
    } else {
      sigma_new <- sigma_old + rnorm(n=1, mean=0, sd=0.01)
      sigma_new <- abs(sigma_new)
      p_obs_new <- p_obs_old + rnorm(n=1, mean=0, sd=0.01)
      while (p_obs_new < 0 | p_obs_new > 1) {
        if (p_obs_new > 1) {
          p_obs_new <- 1 - p_obs_new
        }
        if (p_obs_new < 0) {
          p_obs_new <- -p_obs_new
        }
      }
    }

    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dunif(x = p_obs_new, min = 0, max = 1, log = T)
#    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dbeta(x = p_obs_new, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
    f_hat_new <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_new, proportion_obs = p_obs_new, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }
  }

  sigma_mean <- 0
  p_obs_mean <- 0
  for (i in 501:1000) {
    if (print == T) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #if no data observed, then p_obs is always 0
    if (sum_noisy == 0) {
      sigma_new <- sigma_old + rnorm(n=1, mean=0, sd=0.01)
      sigma_new <- abs(sigma_new)
      p_obs_new <- 0
    } else {
      sigma_new <- sigma_old + rnorm(n=1, mean=0, sd=0.01)
      sigma_new <- abs(sigma_new)
      p_obs_new <- p_obs_old + rnorm(n=1, mean=0, sd=0.01)
      while (p_obs_new < 0 | p_obs_new > 1) {
        if (p_obs_new > 1) {
          p_obs_new <- 1 - p_obs_new
        }
        if (p_obs_new < 0) {
          p_obs_new <- -p_obs_new
        }
      }
    }

    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dunif(x = p_obs_new, min = 0, max = 1, log = T)
    #    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T) + dbeta(x = p_obs_new, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
    f_hat_new <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_new, proportion_obs = p_obs_new, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }

    sigma_mean <- ((i-501)*sigma_mean + sigma_old)/(i-500)
    p_obs_mean <- ((i-501)*p_obs_mean + p_obs_old)/(i-500)
  }

  R <- 100
  Ns <- 10000
  nopt_llik <- rep(NA, R)
  for (r in 1:R) {
    if (print == T) {
      print(r)
    }
    nopt_llik[r] <- sir_mix(n_particles = Ns, sigma = sigma_mean, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_mean, genetic_data = genetic_data, ess_threshold = Ns/2, resampling_scheme = resampling_scheme, backward_sim = FALSE)$int_llik
  }
  var <- sum(nopt_llik^2)/R - mean(nopt_llik)^2
  Nopt <- ceiling(Ns * (var) / (0.92^2))

  n_particles <- max(1, Nopt)

  return(n_particles)
}
