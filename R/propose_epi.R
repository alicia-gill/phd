propose_epi <- function(birth_rate, death_rate, stop_time, ptree, ...) {
  epi <- epidemic(birth_rate = birth_rate, death_rate = death_rate, stop_time = stop_time)
  prev <- prevalence(epidemic = epi, stop_time = stop_time)
  fhat <- genetic_llik(birth_rate = birth_rate, ptree = ptree, prevalence = prev, stop_time = stop_time)
  return(fhat)
}
