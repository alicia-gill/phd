#' Plotting functions
#'
#' @name plot_fns
#' @title A collection of functions to plot inference from MCMC chains
#'
#' @param chain MCMC chain
#' @param death_rate death rate of the epidemic
#' @param q coverage of credible intervals; default 95%
#' @param burn_in burn in
#'
#' @return plot
NULL

#' @rdname plot_fns
#' @export
plot_birthrates <- function(chain, q=0.95, burn_in=0) {
  if (burn_in > 0) {
    birth_rates <- chain$birth_rate[-(1:burn_in),]
  } else {
    birth_rates <- chain$birth_rate
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- min(matrixStats::colQuantiles(birth_rates, probs=lq))
  max <- max(matrixStats::colQuantiles(birth_rates, probs=uq))
  par(mar=c(4,4,1.5,1.5))
  plot(apply(birth_rates, 2, mean), type="l", ylim=c(min,max), xlab="Day", ylab="B(t)")
  lines(matrixStats::colQuantiles(birth_rates, probs=lq), col="blue")
  lines(matrixStats::colQuantiles(birth_rates, probs=uq), col="blue")
}

#' @rdname plot_fns
#' @export
plot_rt <- function(chain, death_rate, q=0.95, burn_in=0) {
  if (burn_in > 0) {
    rt <- chain$birth_rate[-(1:burn_in),]/death_rate
  } else {
    rt <- chain$birth_rate/death_rate
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- min(matrixStats::colQuantiles(rt, probs=lq))
  max <- max(matrixStats::colQuantiles(rt, probs=uq))
  par(mar=c(4,4,1.5,1.5))
  plot(apply(rt, 2, mean), type="l", ylim=c(min,max), xlab="Day", ylab="R(t)")
  lines(matrixStats::colQuantiles(rt, probs=lq), col="blue")
  lines(matrixStats::colQuantiles(rt, probs=uq), col="blue")
}

#' @rdname plot_fns
#' @export
plot_prevalence <- function(chain, q=0.95, burn_in=0) {
  if (burn_in > 0) {
    prevalence <- chain$prevalence[-(1:burn_in),]
  } else {
    prevalence <- chain$prevalence
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- 0
  max <- max(matrixStats::colQuantiles(prevalence, probs=uq))
  n <- ncol(prevalence)-1
  par(mar=c(4,4,1.5,1.5))
  plot(0:n, apply(prevalence, 2, mean), type="l", ylim=c(min,max), xlab="Day", ylab="Prevalence")
  lines(0:n, matrixStats::colQuantiles(prevalence, probs=lq), col="blue")
  lines(0:n, matrixStats::colQuantiles(prevalence, probs=uq), col="blue")
}

#' @rdname plot_fns
#' @export
plot_pobs <- function(chain, burn_in=0) {
  pobs <- chain$proportion_obs
  if (burn_in > 0) {
    min <- min(pobs[-(1:burn_in)])
    max <- max(pobs[-(1:burn_in)])
  } else {
    min <- min(pobs)
    max <- max(pobs)
  }
  iter <- length(pobs)
  par(mar=c(4,4,1.5,1.5))
  plot(pobs, type="l", ylim=c(min,max), xlab="Iteration", ylab="Reporting probability", xaxt="n")
  at <- seq(0,iter,length.out=6)
  labels <- as.character(format(at, scientific=F))
  axis(1, at=at, labels=labels)
}

#' @rdname plot_fns
#' @export
plot_sigma <- function(chain, burn_in=0) {
  sigma <- chain$sigma
  if (burn_in > 0) {
    min <- min(sigma[-(1:burn_in)])
    max <- max(sigma[-(1:burn_in)])
  } else {
    min <- min(sigma)
    max <- max(sigma)
  }
  iter <- length(sigma)
  par(mar=c(4,4,1.5,1.5))
  plot(sigma, type="l", ylim=c(min,max), xlab="Iteration", ylab="Sigma", xaxt="n")
  at <- seq(0,iter,length.out=6)
  labels <- as.character(format(at, scientific=F))
  axis(1, at=at, labels=labels)
}
