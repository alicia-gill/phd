#' Truncated Skellam distribution
#'
#' Random generation truncated at specified location for the Skellam distribution
#'
#' @param n number of observations.
#' @param old_x old prevalence.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param trunc truncation point; -old_x by default to prevent the epidemic dying out.
#'
#' @return log-likelihood.
#' @export
#'
#' @examples
#' rtskellam_smc(n = 100, old_x = x_resample, birth_rate = b_sample, death_rate = 0.1)
rtskellam_smc <- function(n, old_x, birth_rate, death_rate, trunc = -old_x) {
  N <- length(old_x)
  B <- length(birth_rate)
  D <- length(death_rate)
  if (N != n) {
    stop("length of x_resample not the same as n_particles")
  }
  if (B != n) {
    stop("length of b_sample not the same as n_particles")
  }
  if (D != 1) {
    stop("length of death rate should be 1")
  }

  mu1 <- old_x * birth_rate
  mu2 <- old_x * death_rate
  mu <- c(mu1, mu2)

  output <- rep(NA, n)

  index <- (1:n)[is.na(output)]
  m <- length(index)
  count <- 0
  while(m > 0 & count <= 100) {
    bd <- rpois(2*m, mu[c(index, n+index)])
    x <- bd[1:m] - bd[(m+1):(2*m)]
    success <- index[x > trunc[index]]
    output[success] <- x[x > trunc[index]]

    index <- (1:n)[is.na(output)]
    m <- length(index)
    count <- count + 1
  }
  if (m > 0) {
    for (i in 1:m) {
      j <- index[i]
      w <- smc_skellam(new_x = (trunc[j]+1):(trunc[j]+100), rep(old_x[j], 100), birth_rate[j], death_rate, log=T)
      w <- w - matrixStats::logSumExp(w)
      output[j] <- sample((trunc[j]+1):(trunc[j]+100), 1, prob=exp(w))
    }
  }

  return(output)
}
