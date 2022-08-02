skellam_llik <- function(birth_rate, death_rate, prevalence, log=T) {
  b <- birth_rate
  d <- death_rate
  prev <- prevalence[,2]

  n <- length(prev)
  I <- 0
  root <- sqrt(birth_rate*death_rate)
  for (i in 2:n) {
    previ <- prev[i-1]
    I <- I + log(besselI(x=2*previ*root,nu=abs(prev[i]-previ),expon.scaled=T)) + 2*previ*root
  }

  logl <- -(birth_rate+d)*(sum(prev)-prev[n]) + ((prev[n]-prev[1])/2)*(log(birth_rate)-log(death_rate)) + I

  if (log==T) {
    return(logl)
  } else {
    return(exp(logl))
  }
}
