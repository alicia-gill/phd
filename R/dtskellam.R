#' Density function for the truncated Skellam distribution
#'
#' Density function for the truncated Skellam distribution.
#'
#' @param x vector of quantities.
#' @param mu1,mu2 positive valued parameters.
#' @param a,b lower and upper truncation points (a < x <= b).
#' @param log logical; if TRUE, probabilities are given as log(p).
#'
#' @return density.
#' @export
#'
#' @examples
#' dtskellam(x = 10, a=0, mu1 = 20, mu2 = 10)
dtskellam <- function(x, a=-1e150, b=1e300, mu1, mu2, log=T) {
  N <- length(x)
  numer <- rep(NA, N)
  for (i in 1:N) {
    if (x[i] <= a | x[i] > b) {
      numer[i] <- -Inf
    } else {
      numer[i] <- skellam::dskellam(x=x[i], lambda1=mu1, lambda2=mu2, log=T)
    }
  }
  denom <- skellam::pskellam(q=b, lambda1=mu1, lambda2=mu2) - skellam::pskellam(q=a, lambda1=mu1, lambda2=mu2)
  output <- numer - log(denom)
  if (log) {
    return(output)
  } else {
    return(exp(output))
  }
}
