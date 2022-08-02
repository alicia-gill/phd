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
