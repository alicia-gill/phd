skellam_llik <- function(birth_rate, death_rate, dt, prevalence, log=T) {
  b <- birth_rate*dt
  d <- death_rate*dt
  prev <- prevalence[,2]

  n <- length(prev)
  I <- 0
  root <- sqrt(b*d)
  for (i in 2:n) {
    previ <- prev[i-1]
    I <- I + log(besselI(x=2*previ*root,nu=abs(prev[i]-previ),expon.scaled=T)) + 2*previ*root
  }

  logl <- -(b+d)*(sum(prev)-prev[n]) + ((prev[n]-prev[1])/2)*(log(b)-log(d)) + I

  if (log==T) {
    return(logl)
  } else {
    return(exp(logl))
  }
}
