#' Find optimal number of particles
#'
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs0 initial value of the proportion of cases observed.
#' @param x0 prevalence on day 0.
#' @param death_rate death rate of the epidemic.
#' @param ptree object of class phylo.
#' @param day number of days in the past the most recent leaf was sampled.
#' @param noisy_prevalence data frame of observed prevalence per day.
#' @param sigma_mean exponential prior mean of sigma.
#' @param pobs_prior "uniform" or "beta"; prior distribution on reporting probability.
#' @param pobs_min minimum for uniform reporting probability prior.
#' @param pobs_max maximum for uniform reporting probability prior.
#' @param pobs_alpha shape1 for beta reporting probability prior.
#' @param pobs_beta shape2 for beta reporting probability prior.
#' @param x0_prior "uniform" or "nbinom"; prior distribution on day 0 prevalence.
#' @param x0_min minimum for uniform day 0 prevalence prior.
#' @param x0_max maximum for uniform day 0 prevalence prior.
#' @param x0_mean mean for negative binomial day 0 prevalence prior.
#' @param x0_var variance for negative binomial day 0 prevalence prior.
#' @param resampling_scheme "multinomial" or "systematic".
#' @param print logical; if TRUE, prints percentage of the way through the chain.
#'
#' @return optimal number of particles
#' @export
#'
#' @examples
#' find_nopt(sigma0 = 0.1, proportion_obs0 = 0.5, death_rate = 0.1, ptree = sample_tree, day = 1, noisy_prevalence = noisy_prev, print = T)
find_nopt <- function(sigma0, proportion_obs0, x0 = 1, death_rate, ptree, day = 0, noisy_prevalence, sigma_mean = 0.1, pobs_prior = "uniform", pobs_min = 0, pobs_max = 1, pobs_alpha = 1, pobs_beta = 1, x0_prior = "uniform", x0_min=1, x0_max=Inf, x0_mean=10, x0_var=100, ess_threshold_prop = 0.5, resampling_scheme = "systematic", print = F) {
  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)

  sigma_old <- sigma0
  x0_old <- x0

  #setup for adaptation
  s <- 0
  zeta <- Inf #can be any constant >= 1
  inv_zeta <- 1/zeta

  sum_noisy <- sum(noisy_prevalence[-1,2])

  #if sum_noisy=0, then pobs=0, so there are only 2 hyperparameters
  if (sum_noisy == 0) {
    p_obs_old <- 0
    mu_old <- c(sigma_old, x0_old)
    Sigma_old <- diag(c(0.01, max(1, round(x0_old/10,0))))
    # sqrtSigma_old <- expm::sqrtm(Sigma_old)
    trace <- Sigma_old[1,1] + Sigma_old[2,2]
    det <- Sigma_old[1,1]*Sigma_old[2,2] - Sigma_old[1,2]*Sigma_old[2,1]
    sqrtdet <- sqrt(det)
    tt <- sqrt(trace + 2*sqrtdet)
    sqrtSigma_old <- (Sigma_old + sqrtdet*diag(1,nrow=2,ncol=2))/tt
  } else {
    p_obs_old <- proportion_obs0
    mu_old <- c(sigma_old, p_obs_old, x0_old)
    Sigma_old <- diag(c(0.01, 0.01, max(1, round(x0_old/10,0))))
    sqrtSigma_old <- expm::sqrtm(Sigma_old)
  }

  n_particles <- 10000
  ess_threshold <- n_particles*ess_threshold_prop

  #prior on sigma is exponential
  #prior on x0-1 is uniform(1,Inf)
  #prior on p_obs is uniform(0,1) or beta(1,3)
  sigma_prior_old <- dexp(x = sigma_old, rate = sigma_mean, log = T)
  if (pobs_prior == "uniform") {
    pobs_prior_old <- dunif(x = p_obs_old, min = pobs_min, max = pobs_max, log = T)
  } else if (pobs_prior == "beta") {
    pobs_prior_old <- dbeta(x = p_obd_ols, shape1 = pobs_alpha, shape2 = pobs_beta, log = T)
  } else {
    stop("Reporting probability prior should be 'uniform' or 'beta'")
  }
  if (x0_prior == "uniform") {
    if (x0_max <= 1e300) {
      x0_prior_old <- dunif(x = x0_old, min = x0_min, max = x0_max, log = T)
    } else {
      x0_prior_old <- 0
    }
  } else if (x0_prior == "nbinom") {
    if (x0_var <= x0_mean) {
      stop("X0 variance must be larger than X0 mean")
    }
    size_nb <- x0_mean * x0_mean / (x0_var - x0_mean)
    x0_prior_old <- dnbinom(x = x0_old, size = size_nb, mu = x0_mean, log = T)
  } else {
    stop("Day 0 prevalence prior should be 'uniform' or 'nbinom'")
  }
  prior_old <- sigma_prior_old + pobs_prior_old + x0_prior_old
  f_hat_old <- sir_mix(n_particles = n_particles, sigma = sigma_old, x0 = x0_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_old, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik

  #need to find a good estimate of the posterior mean for sigma and p_obs
  #first 500 for burn-in
  #second 500 to estimate means
  for (i in 1:500) {
#    print(i)
    if (print) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #if no data observed, then p_obs is always 0
    if (sum_noisy == 0) {
      w <- rnorm(n = 2, mean = 0, sd = 1)
      move <- exp(s) * sqrtSigma_old %*% w
      sigma_new <- abs(sigma_old + move[1])
      x0_move <- round(move[2],0)
      x0_sign <- sign(x0_move)
      x0_move <- x0_sign*max(1, abs(x0_move))
      x0_new <- abs(x0_old + x0_move - 1) + 1
      # x0_new <- abs(x0_old + round(move[2], 0) - 1) + 1
      p_obs_new <- 0
    } else {
      w <- rnorm(n = 3, mean = 0, sd = 1)
      move <- exp(s) * sqrtSigma_old %*% w
      sigma_new <- abs(sigma_old + move[1])
      p_obs_new <- abs(p_obs_old + move[2])
      x0_move <- round(move[3],0)
      x0_sign <- sign(x0_move)
      x0_move <- x0_sign*max(1, abs(x0_move))
      x0_new <- abs(x0_old + x0_move - 1) + 1
      # x0_new <- abs(x0_old + round(move[3], 0) - 1) + 1
      while (p_obs_new < 0 | p_obs_new > 1) {
        if (p_obs_new > 1) {
          p_obs_new <- 1 - p_obs_new
        }
        if (p_obs_new < 0) {
          p_obs_new <- -p_obs_new
        }
      }
    }

    sigma_prior_new <- dexp(x = sigma_new, rate = sigma_mean, log = T)
    if (pobs_prior == "uniform") {
      pobs_prior_new <- dunif(x = p_obs_new, min = pobs_min, max = pobs_max, log = T)
    } else if (pobs_prior == "beta") {
      pobs_prior_new <- dbeta(x = p_obd_ols, shape1 = pobs_alpha, shape2 = pobs_beta, log = T)
    } else {
      stop("Reporting probability prior should be 'uniform' or 'beta'")
    }
    if (x0_prior == "uniform") {
      if (x0_max <= 1e300) {
        x0_prior_new <- dunif(x = x0_new, min = x0_min, max = x0_max, log = T)
      } else {
        x0_prior_new <- 0
      }
    } else if (x0_prior == "nbinom") {
      x0_prior_new <- dnbinom(x = x0_new, size = size_nb, mu = x0_mean, log = T)
    } else {
      stop("Day 0 prevalence prior should be 'uniform' or 'nbinom'")
    }
    prior_new <- sigma_prior_new + pobs_prior_new + x0_prior_new
    f_hat_new <- sir_mix(n_particles = n_particles, sigma = sigma_new, x0 = x0_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_new, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 100)^(-0.8)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
      x0_old <- x0_new
    }

    if (sum_noisy == 0) {
      Xn <- c(sigma_old, x0_old)
    } else {
      Xn <- c(sigma_old, p_obs_old, x0_old)
    }
    mu_new <- ((1 - eta) * mu_old) + (eta * Xn)
    Sigma_new <- ((1 - eta) * Sigma_old) + (eta * (Xn - mu_old) %*% t(Xn - mu_old))
    evalues <- eigen(Sigma_new, symmetric = T)$values
    min_evalue <- min(evalues)
    max_evalue <- max(evalues)
    if (norm(mu_new, type="2") <= zeta & max_evalue <= zeta & min_evalue >= inv_zeta) {
      mu_old <- mu_new
      Sigma_old <- Sigma_new
      # sqrtSigma_old <- expm::sqrtm(Sigma_old)
      if (sum_noisy == 0) {
        trace <- Sigma_old[1,1] + Sigma_old[2,2]
        det <- Sigma_old[1,1]*Sigma_old[2,2] - Sigma_old[1,2]*Sigma_old[2,1]
        sqrtdet <- sqrt(det)
        tt <- sqrt(trace + 2*sqrtdet)
        sqrtSigma_old <- (Sigma_old + sqrtdet*diag(1,nrow=2,ncol=2))/tt
      } else {
        sqrtSigma_old <- expm::sqrtm(Sigma_old)
      }
    }
  }

  sigma_mean <- 0
  p_obs_mean <- 0
  x0_mean <- 0
  for (i in 501:1000) {
    if (print == T) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    #if no data observed, then p_obs is always 0
    if (sum_noisy == 0) {
      w <- rnorm(n = 2, mean = 0, sd = 1)
      move <- exp(s) * sqrtSigma_old %*% w
      sigma_new <- abs(sigma_old + move[1])
      x0_move <- round(move[2],0)
      x0_sign <- sign(x0_move)
      x0_move <- x0_sign*max(1, abs(x0_move))
      x0_new <- abs(x0_old + x0_move - 1) + 1
      # x0_new <- abs(x0_old + round(move[2], 0) - 1) + 1
      p_obs_new <- 0
    } else {
      w <- rnorm(n = 3, mean = 0, sd = 1)
      move <- exp(s) * sqrtSigma_old %*% w
      sigma_new <- abs(sigma_old + move[1])
      p_obs_new <- abs(p_obs_old + move[2])
      x0_move <- round(move[3],0)
      x0_sign <- sign(x0_move)
      x0_move <- x0_sign*max(1, abs(x0_move))
      x0_new <- abs(x0_old + x0_move - 1) + 1
      # x0_new <- abs(x0_old + round(move[3], 0) - 1) + 1
      while (p_obs_new < 0 | p_obs_new > 1) {
        if (p_obs_new > 1) {
          p_obs_new <- 1 - p_obs_new
        }
        if (p_obs_new < 0) {
          p_obs_new <- -p_obs_new
        }
      }
    }

    sigma_prior_new <- dexp(x = sigma_new, rate = sigma_mean, log = T)
    if (pobs_prior == "uniform") {
      pobs_prior_new <- dunif(x = p_obs_new, min = pobs_min, max = pobs_max, log = T)
    } else if (pobs_prior == "beta") {
      pobs_prior_new <- dbeta(x = p_obd_ols, shape1 = pobs_alpha, shape2 = pobs_beta, log = T)
    } else {
      stop("Reporting probability prior should be 'uniform' or 'beta'")
    }
    if (x0_prior == "uniform") {
      if (x0_max <= 1e300) {
        x0_prior_new <- dunif(x = x0_new, min = x0_min, max = x0_max, log = T)
      } else {
        x0_prior_new <- 0
      }
    } else if (x0_prior == "nbinom") {
      p_nb <- x0_mean / x0_var
      n_nb <- x0_mean * p_nb / (1-p_nb)
      x0_prior_new <- dnbinom(x = x0_new, size = n_nb, prob = p_nb, log = T)
    } else {
      stop("Day 0 prevalence prior should be 'uniform' or 'nbinom'")
    }
    prior_new <- sigma_prior_new + pobs_prior_new + x0_prior_new
    f_hat_new <- sir_mix(n_particles = n_particles, sigma = sigma_new, x0 = x0_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_new, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 100)^(-0.8)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
      p_obs_old <- p_obs_new
      x0_old <- x0_new
    }

    sigma_mean <- ((i-501)*sigma_mean + sigma_old)/(i-500)
    p_obs_mean <- ((i-501)*p_obs_mean + p_obs_old)/(i-500)
    x0_mean <- ((i-501)*x0_mean + x0_old)/(i-500)

    if (sum_noisy == 0) {
      Xn <- c(sigma_old, x0_old)
    } else {
      Xn <- c(sigma_old, p_obs_old, x0_old)
    }
    mu_new <- ((1 - eta) * mu_old) + (eta * Xn)
    Sigma_new <- ((1 - eta) * Sigma_old) + (eta * (Xn - mu_old) %*% t(Xn - mu_old))
    evalues <- eigen(Sigma_new, symmetric = T)$values
    min_evalue <- min(evalues)
    max_evalue <- max(evalues)
    if (norm(mu_new, type="2") <= zeta & max_evalue <= zeta & min_evalue >= inv_zeta) {
      mu_old <- mu_new
      Sigma_old <- Sigma_new
      # sqrtSigma_old <- expm::sqrtm(Sigma_old)
      if (sum_noisy == 0) {
        trace <- Sigma_old[1,1] + Sigma_old[2,2]
        det <- Sigma_old[1,1]*Sigma_old[2,2] - Sigma_old[1,2]*Sigma_old[2,1]
        sqrtdet <- sqrt(det)
        tt <- sqrt(trace + 2*sqrtdet)
        sqrtSigma_old <- (Sigma_old + sqrtdet*diag(1,nrow=2,ncol=2))/tt
      } else {
        sqrtSigma_old <- expm::sqrtm(Sigma_old)
      }
    }
  }

  R <- 100
  Ns <- 1000
  nopt_llik <- rep(NA, R)
  for (r in 1:R) {
    if (print == T) {
      print(r)
    }
    nopt_llik[r] <- sir_mix(n_particles = Ns, x0 = min(1, round(x0_mean,0)), sigma = sigma_mean, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = p_obs_mean, genetic_data = genetic_data, ess_threshold = Ns*ess_threshold_prop, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
  }
  var <- sum(nopt_llik^2)/R - mean(nopt_llik)^2
  Nopt <- ceiling(Ns * (var) / (0.92^2))

  n_particles <- max(1, Nopt)

  return(n_particles)
}
