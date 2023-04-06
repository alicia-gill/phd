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
  D <- length(death_rate)
  output <- rep(NA, n)
  for (i in 1:n) {
    j <- i %% N
    if (j == 0) {
      j <- N
    }
    k <- i %% B
    if (k == 0) {
      k <- B
    }
    l <- i %% D
    if (l == 0) {
      l <- D
    }
    oldx <- old_x[j]
    mu1 <- oldx*birth_rate[k]
    mu2 <- oldx*death_rate[l]
    mu <- c(mu1, mu2)
    #truncate proposals st epidemic cannot die out
    t <- -oldx
    for (count in 1:1000) {
      bd <- rpois(2, mu)
      x <- bd[1] - bd[2]
      if (x>t) {
        output[i] <- x
        break
      }
    }
    if (count == 1000) {
      w <- smc_skellam((t+1):(t+100), rep(oldx, 100), birth_rate[k], death_rate[l], log=T)
      w <- w - matrixStats::logSumExp(w)
      output[i] <- sample((t+1):(t+100), 1, prob=exp(w))
    }
  }
  return(output)
}
