#' Density function for the truncated negative binomial distribution
#'
#' Density function for the truncated negative binomial distribution.
#'
#' @param x vector of quantities.
#' @param size dispersion parameter.
#' @param mu mean.
#' @param a,b lower and upper truncation points (a <= x <= b).
#' @param log logical; if TRUE, probabilities are given as log(p).
#'
#' @return density.
#' @export
#'
#' @examples
#' dtnbinom(x = 1, size = 10, mu = 2, a = 1, b = 100)
dtnbinom <- function(x, size, mu, a, b, log = TRUE) {
  log_num <- dnbinom(x = x, size = size, mu = mu, log = TRUE)
  denominator <- 1 - pnbinom(a-1, size = size, mu = mu) - pnbinom(b, size = size, mu = mu, lower.tail = FALSE)
  log_denom <- log(denominator)

  if (x < a | x > b) {
    output <- -Inf
  } else {
    output <- log_num - log_denom
  }

  if (log) {
    return(output)
  } else {
    return(exp(output))
  }
}
