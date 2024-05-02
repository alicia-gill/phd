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
plot_birthrates <- function(chain, q=0.95, burn_in=0, start_time=NA, end_time=NA) {
  if (burn_in > 0) {
    birth_rates <- chain$birth_rate[-(1:burn_in),]
  } else {
    birth_rates <- chain$birth_rate
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- min(matrixStats::colQuantiles(birth_rates, probs=lq, na.rm=T))
  max <- max(matrixStats::colQuantiles(birth_rates, probs=uq, na.rm=T))
  ncol <- ncol(birth_rates)
  if (is.na(start_time) & is.na(end_time)) {
    warning("Specify start or end time")
    break
  }
  if (is.na(end_time)) {
    end_time <- start_time + ncol - 1
  }
  if (is.na(start_time)) {
    start_time <- end_time - ncol + 1
  }
  par(mar=c(4,4,1.5,1.5))
  plot(apply(birth_rates, 2, mean, na.rm=T), type="l", xlim=c(1, ncol), ylim=c(min,max), xlab="Time", ylab="B(t)", xaxt="n")
  axis(1, at=seq(0, ncol, by=10), labels=seq(start_time-1, end_time, by=10))
  lines(matrixStats::colQuantiles(birth_rates, probs=lq, na.rm=T), col="blue")
  lines(matrixStats::colQuantiles(birth_rates, probs=uq, na.rm=T), col="blue")
}

#' @rdname plot_fns
#' @export
plot_rt <- function(chain, death_rate, q=0.95, burn_in=0, start_time=NA, end_time=NA) {
  if (burn_in > 0) {
    rt <- chain$birth_rate[-(1:burn_in),]/death_rate
  } else {
    rt <- chain$birth_rate/death_rate
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- min(matrixStats::colQuantiles(rt, probs=lq, na.rm=T))
  max <- max(matrixStats::colQuantiles(rt, probs=uq, na.rm=T))
  ncol <- ncol(rt)
  if (is.na(start_time) & is.na(end_time)) {
    warning("Specify start or end time")
    break
  }
  if (is.na(end_time)) {
    end_time <- start_time + ncol - 1
  }
  if (is.na(start_time)) {
    start_time <- end_time - ncol + 1
  }
  par(mar=c(4,4,1.5,1.5))
  plot(apply(rt, 2, mean, na.rm=T), type="l", xlim=c(1, ncol), ylim=c(min,max), xlab="Time", ylab="R(t)", xaxt="n")
  axis(1, at=seq(0, ncol, by=10), labels=seq(start_time-1, end_time, by=10))
  lines(matrixStats::colQuantiles(rt, probs=lq, na.rm=T), col="blue")
  lines(matrixStats::colQuantiles(rt, probs=uq, na.rm=T), col="blue")
}

#' @rdname plot_fns
#' @export
plot_prevalence <- function(chain, q=0.95, burn_in=0, start_time=NA, end_time=NA) {
  if (burn_in > 0) {
    prevalence <- chain$prevalence[-(1:burn_in),]
  } else {
    prevalence <- chain$prevalence
  }
  lq <- (1-q)/2
  uq <- 1-lq
  min <- 0
  max <- max(matrixStats::colQuantiles(prevalence, probs=uq, na.rm=T))
  ncol <- ncol(prevalence)
  n <- ncol-1
  if (is.na(start_time) & is.na(end_time)) {
    warning("Specify start or end time")
    break
  }
  if (is.na(end_time)) {
    end_time <- start_time + ncol
  }
  if (is.na(start_time)) {
    start_time <- end_time - ncol + 1
  }
  par(mar=c(4,4,1.5,1.5))
  plot(0:n, apply(prevalence, 2, mean, na.rm=T), type="l", xlim=c(0, ncol-1), ylim=c(min,max), xlab="Time", ylab="Prevalence", xaxt="n")
  axis(1, at=seq(0, ncol, by=10), labels=seq(start_time, end_time, by=10))
  lines(0:n, matrixStats::colQuantiles(prevalence, probs=lq, na.rm=T), col="blue")
  lines(0:n, matrixStats::colQuantiles(prevalence, probs=uq, na.rm=T), col="blue")
}

#' @rdname plot_fns
#' @export
plot_pobs <- function(chain, burn_in=0, iter=NA) {
  pobs <- chain$proportion_obs
  if (burn_in > 0) {
    min <- min(pobs[-(1:burn_in)], na.rm=T)
    max <- max(pobs[-(1:burn_in)], na.rm=T)
  } else {
    min <- min(pobs, na.rm=T)
    max <- max(pobs, na.rm=T)
  }
  if (is.na(iter)) {
    which_max <- which.max(is.na(chain$sigma))-1
    iter <- ifelse(which_max==0, length(pobs), which_max)
    iter <- 5*ceiling(iter/5)
  }
  par(mar=c(4,4,1.5,1.5))
  plot(pobs, type="l", xlim=c(0,iter), ylim=c(min,max), xlab="Iteration", ylab="Reporting probability", xaxt="n")
  at <- seq(0,iter,length.out=6)
  labels <- as.character(format(at, scientific=F))
  axis(1, at=at, labels=labels)
}

#' @rdname plot_fns
#' @export
plot_sigma <- function(chain, burn_in=0, iter=NA) {
  sigma <- chain$sigma
  if (burn_in > 0) {
    min <- min(sigma[-(1:burn_in)], na.rm=T)
    max <- max(sigma[-(1:burn_in)], na.rm=T)
  } else {
    min <- min(sigma, na.rm=T)
    max <- max(sigma, na.rm=T)
  }
  if (is.na(iter)) {
    which_max <- which.max(is.na(sigma))-1
    iter <- ifelse(which_max==0, length(sigma), which_max)
    iter <- 5*ceiling(iter/5)
  }
  par(mar=c(4,4,1.5,1.5))
  plot(sigma, type="l", xlim=c(0,iter), ylim=c(min,max), xlab="Iteration", ylab="Sigma", xaxt="n")
  at <- seq(0,iter,length.out=6)
  labels <- as.character(format(at, scientific=F))
  axis(1, at=at, labels=labels)
}

#' @rdname plot_fns
#' @export
plot_x0 <- function(chain, burn_in=0, iter=NA) {
  x0 <- chain$x0
  if (burn_in > 0) {
    min <- min(x0[-(1:burn_in)], na.rm=T)
    max <- max(x0[-(1:burn_in)], na.rm=T)
  } else {
    min <- min(x0, na.rm=T)
    max <- max(x0, na.rm=T)
  }
  if (is.na(iter)) {
    which_max <- which.max(is.na(x0))-1
    iter <- ifelse(which_max==0, length(x0), which_max)
    iter <- 5*ceiling(iter/5)
  }
  par(mar=c(4,4,1.5,1.5))
  plot(x0, type="l", xlim=c(0,iter), ylim=c(min,max), xlab="Iteration", ylab="X0", xaxt="n")
  at <- seq(0,iter,length.out=6)
  labels <- as.character(format(at, scientific=F))
  axis(1, at=at, labels=labels)
}
