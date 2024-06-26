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
#' find_nopt(sigma0 = 0.1, proportion_obs0 = 0.5, death_rate = 0.1, ptree = sample_tree, day = 1, noisy_prevalence = noisy_prev, print = T)
find_nopt <- function(sigma0, proportion_obs0, x0 = 1, death_rate, ptree, day = 0, noisy_prevalence, ess_threshold_prop = 0.5, resampling_scheme = "systematic", backward_sim = F, print = F) {
  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)

  sigma_old <- sigma0

  s <- 0
  zeta <- 10000 #can be any constant >= 1
  inv_zeta <- 1/zeta

  lambda_sigma <- 1/0.1 #mean of 0.1
  #alpha_pobs <- 1
  #beta_pobs <- 3

  sum_noisy <- sum(noisy_prevalence[-1,2])

  if (sum_noisy == 0) {
    p_obs_old <- 0
    mu_old <- sigma_old
    Sigma_old <- 1
    sqrtSigma_old <- sqrt(Sigma_old)
  } else {
    p_obs_old <- proportion_obs0
    mu_old <- c(sigma_old, p_obs_old)
    Sigma_old <- diag(1, nrow=2, ncol=2)
    sqrtSigma_old <- expm::sqrtm(Sigma_old)
  }

  n_particles <- 1000
  ess_threshold <- n_particles*ess_threshold_prop

  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dunif(x = p_obs_old, min = 0, max = 1, log = T)
#  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T) + dbeta(x = p_obs_old, shape1 = alpha_pobs, shape2 = beta_pobs, log = T)
  f_hat_old <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_old, proportion_obs = p_obs_old, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = backward_sim)$int_llik

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
      w <- rnorm(n = 1, mean = 0, sd = 1)
      new <- sigma_old + exp(s) * sqrtSigma_old %*% w
      new <- abs(new)
      sigma_new <- new
      p_obs_new <- 0
    } else {
      w <- rnorm(n=2, mean=0, sd=0.01)
      #      w <- MASS::mvrnorm(n = 1, mu = rep(0, 2), Sigma = diag(0.01, nrow=2))
      new <- c(sigma_old, p_obs_old) + exp(s) * sqrtSigma_old %*% w
      new <- abs(new)
      sigma_new <- new[1]
      p_obs_new <- new[2]
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
    f_hat_new <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_new, proportion_obs = p_obs_new, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = backward_sim)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 100)^(-0.9)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }

    if (sum_noisy == 0) {
      Xn <- sigma_old
    } else {
      Xn <- c(sigma_old, p_obs_old)
    }
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
      w <- rnorm(n = 1, mean = 0, sd = 1)
      new <- sigma_old + exp(s) * sqrtSigma_old %*% w
      new <- abs(new)
      sigma_new <- new
      p_obs_new <- 0
    } else {
      w <- rnorm(n=2, mean=0, sd=0.01)
#      w <- MASS::mvrnorm(n = 1, mu = rep(0, 2), Sigma = diag(0.01, nrow=2))
      new <- c(sigma_old, p_obs_old) + exp(s) * sqrtSigma_old %*% w
      new <- abs(new)
      sigma_new <- new[1]
      p_obs_new <- new[2]
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
    f_hat_new <- sir_mix(n_particles = n_particles, ess_threshold = ess_threshold, x0 = x0, death_rate = death_rate, sigma = sigma_new, proportion_obs = p_obs_new, noisy_prevalence = noisy_prevalence, genetic_data = genetic_data, resampling_scheme = resampling_scheme, backward_sim = backward_sim)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 100)^(-0.9)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
    }
    sigma_mean <- ((i-501)*sigma_mean + sigma_old)/(i-500)
    p_obs_mean <- ((i-501)*p_obs_mean + p_obs_old)/(i-500)

    if (sum_noisy == 0) {
      Xn <- sigma_old
    } else {
      Xn <- c(sigma_old, p_obs_old)
    }
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
  }

  R <- 100
  Ns <- 1000
  nopt_llik <- rep(NA, R)
  for (r in 1:R) {
    if (print == T) {
      print(r)
    }
    nopt_llik[r] <- sir_mix(n_particles = Ns, sigma = sigma_mean, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_mean, genetic_data = genetic_data, ess_threshold = Ns*ess_threshold_prop, resampling_scheme = resampling_scheme, backward_sim = backward_sim)$int_llik
  }
  var <- sum(nopt_llik^2)/R - mean(nopt_llik)^2
  Nopt <- ceiling(Ns * (var) / (0.92^2))

  n_particles <- max(1, Nopt)

  return(n_particles)
}
