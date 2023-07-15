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
  x <- 2 * old_x * sqrt(birth_rate * death_rate)
  length <- length(nu)
  logI <- rep(NA, length)
  nu_cutoff <- 10
  x_cutoff <- 10000
  set0 <- (1:length)[nu <= nu_cutoff & x <= x_cutoff] #use log(besselI)
  set1 <- (1:length)[nu <= nu_cutoff & x > x_cutoff] #use Bessel::besselIasym for small/moderate nu and large x
  set2 <- (1:length)[nu > nu_cutoff] #use Bessel::besselI.nuAsym for large nu (and x)

  if (length(birth_rate) > 1) {
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate[set0] * death_rate), nu = nu[set0], expon.scaled = T))
    logI[set1] <- Bessel::besselIasym(x = 2 * old_x[set1] * sqrt(birth_rate[set1] * death_rate), nu = nu[set1], k.max = 20, expon.scaled = T, log = T)
    logI[set2] <- Bessel::besselI.nuAsym(x = 2 * old_x[set2] * sqrt(birth_rate[set2] * death_rate), nu = nu[set1], k.max = 5, expon.scaled = T, log = T)
  } else {
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate * death_rate), nu = nu[set0], expon.scaled = T))
    logI[set1] <- Bessel::besselIasym(x = 2 * old_x[set1] * sqrt(birth_rate * death_rate), nu = nu[set1], k.max = 20, expon.scaled = T, log = T)
    logI[set2] <- Bessel::besselI.nuAsym(x = 2 * old_x[set2] * sqrt(birth_rate * death_rate), nu = nu[set1], k.max = 5, expon.scaled = T, log = T)
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
