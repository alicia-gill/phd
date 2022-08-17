#' Propose Epidemic
#'
#' Simulates an epidemic with a given birth rate, death rate and stop time and returns the log-likelihood corresponding to that epidemic.
#'
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param stop_time number of days the epidemic simulation was run for.
#' @param ptree object of class phylo.
#' @param ... allows for the use of lapply to run it multiple times.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' propose_epi(birth_rate = 0.15, death_rate = 0.1, stop_time = 50, ptree = sample_tree)
propose_epi <- function(birth_rate, death_rate, stop_time, ptree, ...) {
  epi <- epidemic(birth_rate = birth_rate, death_rate = death_rate, stop_time = stop_time)
  prev <- prevalence(epidemic = epi, stop_time = stop_time)
  fhat <- genetic_llik(birth_rate = birth_rate, ptree = ptree, prevalence = prev, stop_time = stop_time)
  return(fhat)
}
