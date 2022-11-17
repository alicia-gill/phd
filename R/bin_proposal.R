#' Propose according to binomial distribution
#'
#' Generates proposals for unknown number of trials n given x successes and success probability p.
#'
#' @param n number of particles.
#' @param x number of successes.
#' @param p success probability.
#'
#' @return list containing: number of trials, corresponding log-likelihood.
#' @export
#'
#' @examples
#' bin_proposal(n=100, x=10, p=0.1)
bin_proposal <- function(n, x, p) {
  out_n <- rep(NA, n)
  out_p <- rep(NA, n)
  logp <- log(p)
  logq <- log(1-p)
  for (i in 1:n) {
    success <- 0
    total <- 0
    while (success <= x) {
      u <- runif(1,0,1)
      if (u <= p) {
        success <- success + 1
      }
      total <- total + 1
    }
    out_n[i] <- total - 1
    out_p[i] <- lchoose(out_n[i], x) + (x+1)*logp + (out_n[i]-x)*logq
  }
  return(list("n"=out_n, "loglik"=out_p))
}
