#' Truncated Skellam distribution
#'
#' Random generation truncated at 0 for the Skellam distribution
#'
#' @param n number of observations.
#' @param old_x old prevalence.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' rtskellam(n = 100, old_x = x_resample, birth_rate = b_sample, death_rate = 0.1)
rtskellam <- function(n, old_x, birth_rate, death_rate) {
  N <- length(old_x)
  B <- length(birth_rate)
  output <- rep(NA, n)
  for (i in 1:n) {
    count <- 0
    j <- i %% N
    if (j == 0) {
      j <- N
    }
    mu1 <- (old_x*birth_rate)[j]
    mu2 <- (old_x*death_rate)[j]
    #truncate proposals st epidemic cannot die out
    t <- -old_x[j]
    while (is.na(output[i]) == TRUE & count < 1000) {
      b <- rpois(1, mu1)
      d <- rpois(1, mu2)
      x <- b-d
      if (x>t) {
        output[i] <- x
      }
      count <- count + 1
    }
    if (is.na(output[i]) == TRUE) {
      w <- smc_skellam((t+1):(t+100), rep(old_x[j], 100), birth_rate[j], death_rate, log=T)
      w <- w - matrixStats::logSumExp(w)
      output[i] <- sample((t+1):(t+100), 1, prob=exp(w))
    }
  }
  return(output)
}
