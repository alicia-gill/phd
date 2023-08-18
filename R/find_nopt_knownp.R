#' Find optimal number of particles with known proportion observed
#'
#' @param sigma0 initial value of the linear gaussian standard deviation.
#' @param proportion_obs proportion of cases observed.
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
#' find_nopt_knownp(sigma0 = 0.1, proportion_obs = 0.5, death_rate = 0.1, ptree = sample_tree, day = 1, noisy_prevalence = noisy_prev, print = T)
find_nopt_knownp <- function(sigma0, proportion_obs, death_rate, ptree, day = 0, noisy_prevalence, resampling_scheme = "systematic", print=F) {
  n <- nrow(noisy_prevalence)
  stop_time <- n - 1
  genetic_data <- genetic_data(ptree = ptree, stop_time = stop_time, day = day)

  sigma_old <- sigma0

  s <- 0
  zeta <- 10000 #can be any constant >= 1
  inv_zeta <- 1/zeta

  lambda_sigma <- 1/0.1 #mean of 0.1

  sum_noisy <- sum(noisy_prevalence[-1,2])

  mu_old <- sigma_old
  Sigma_old <- 1
  sqrtSigma_old <- sqrt(Sigma_old)

  n_particles <- 5000
  ess_threshold <- n_particles/2

  prior_old <- dexp(x = sigma_old, rate = lambda_sigma, log = T)
  f_hat_old <- sir_mix(n_particles = n_particles, sigma = sigma_old, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik

  #need to find a good estimate of the posterior mean for sigma
  #first 500 for burn-in
  #second 500 to estimate means
  for (i in 1:500) {
    if (print == T) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }
    w <- rnorm(n = 1, mean = 0, sd = 1)
    new <- sigma_old + exp(s) * sqrtSigma_old %*% w
    new <- abs(new)
    sigma_new <- new

    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T)
    f_hat_new <- sir_mix(n_particles = n_particles, sigma = sigma_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 10)^(-0.6)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
    }

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
  }

  sigma_mean <- 0
  for (i in 501:1000) {
    if (print == T) {
      j <- i/10
      if (j %% 1 == 0) {
        print(paste0(j,"%"))
      }
    }

    w <- rnorm(n = 1, mean = 0, sd = 1)
    new <- sigma_old + exp(s) * sqrtSigma_old %*% w
    new <- abs(new)
    sigma_new <- new

    prior_new <- dexp(x = sigma_new, rate = lambda_sigma, log = T)
    f_hat_new <- sir_mix(n_particles = n_particles, sigma = sigma_new, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = ess_threshold, resampling_scheme = resampling_scheme, backward_sim = F)$int_llik
    logr <- prior_new + f_hat_new - prior_old - f_hat_old
    loga <- min(0,logr)
    a <- exp(loga)
    eta <- (i + 10)^(-0.6)
    #targeting 10% acceptance
    s <- s + (a - 0.1) * eta
    logu <- -rexp(1)
    if (logu <= loga) {
      prior_old <- prior_new
      f_hat_old <- f_hat_new
      sigma_old <- sigma_new
    }
    sigma_mean <- ((i-501)*sigma_mean + sigma_old)/(i-500)

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
  }

  R <- 100
  Ns <- 5000
  nopt_llik <- rep(NA, R)
  for (r in 1:R) {
    nopt_llik[r] <- sir_mix(n_particles = Ns, sigma = sigma_mean, death_rate = death_rate, noisy_prevalence = noisy_prevalence, proportion_obs = proportion_obs, genetic_data = genetic_data, ess_threshold = Ns/2, resampling_scheme = resampling_scheme, backward_sim = FALSE)$int_llik
  }
  var <- sum(nopt_llik^2)/R - mean(nopt_llik)^2
  Nopt <- ceiling(Ns * (var) / (0.92^2))

  n_particles <- max(1, min(Nopt, 5000))

  return(n_particles)
}
