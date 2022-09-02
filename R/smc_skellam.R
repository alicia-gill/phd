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

  logl <- ((new_x - old_x)/2 * (log(birth_rate) - log(death_rate))) -
          (old_x * (birth_rate + death_rate)) +
          log(besselI(x = 2 * old_x * sqrt(birth_rate * death_rate), nu = abs(new_x - old_x), expon.scaled = T)) +
          (2 * old_x * sqrt(birth_rate * death_rate))

  if (log==T) {
    return(logl)
  } else {
    return(exp(logl))
  }
}
