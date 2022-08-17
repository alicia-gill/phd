#' Genetic Log-Likelihood
#'
#' Probability of the phylogenetic tree given the birth rate and the prevalence.
#'
#' @param birth_rate birth rate of the epidemic.
#' @param ptree object of class phylo.
#' @param prevalence data frame of prevalences by day.
#' @param stop_time number of days the epidemic simulation was run for.
#' @param log logical; if TRUE, log-likelihood is given.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' genetic_llik(birth_rate = 0.2, ptree = sample_tree, prevalence = prev, stop_time = 50)
genetic_llik <- function(birth_rate, ptree, prevalence, stop_time, log=T) {
  gen_data <- genetic_data(ptree, stop_time)
  data <- list("prob" = 2 * birth_rate / prevalence[-1,2],
               "n_pairs" = choose(gen_data[-1,2], 2),
               "n_coal" = gen_data[-1,3])

  if (min(data$prob)<0 | max(data$prob)>1) {
    if (log==T) {
      return(-Inf)
    } else {
      return(0)
    }
  }

  llik <- sum(dbinom(data$n_coal,data$n_pairs,data$prob,log=T))

  if (log==T) {
    return(llik)
  } else {
    return(exp(llik))
  }
}
