#' Skellam Log-Likelihood
#'
#' Approximates the probability of a change in prevalence given the birth and death rate.
#'
#' @param new_x new prevalence.
#' @param old_x old prevalence.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param log logical; if TRUE, log-likelihood is given.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' smc_skellam(new_x = 10, old_x = 7, birth_rate = 0.2, death_rate = 0.1)
smc_skellam <- function(new_x, old_x, birth_rate, death_rate, log=T) {

  nu <- abs(new_x - old_x)
  length <- length(nu)
  logI <- rep(NA, length)
  cutoff <- 0
  set0 <- (1:length)[nu <= cutoff]
  set1 <- (1:length)[nu > cutoff]

  if (length(birth_rate) > 1) {
    #if nu=0, then use normal besselI
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate[set0] * death_rate), nu = nu[set0], expon.scaled = T))
    #if nu>0, then use asymptotic besselI
    logI[set1] <- Bessel::besselI.nuAsym(x = 2 * old_x[set1] * sqrt(birth_rate[set1] * death_rate), nu = nu[set1], k.max = 5, expon.scaled = T, log = T)
  } else {
    #if nu=0, then use normal besselI
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate * death_rate), nu = nu[set0], expon.scaled = T))
    #if nu>0, then use asymptotic besselI
    logI[set1] <- Bessel::besselI.nuAsym(x = 2 * old_x[set1] * sqrt(birth_rate * death_rate), nu = nu[set1], k.max = 5, expon.scaled = T, log = T)
  }

  logl <- ((new_x - old_x)/2 * (log(birth_rate) - log(death_rate))) -
          (old_x * (birth_rate + death_rate)) +
          logI +
          (2 * old_x * sqrt(birth_rate * death_rate))

  if (log==T) {
    return(logl)
  } else {
    return(exp(logl))
  }
}
