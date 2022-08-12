t <- 10

VV <- 0.5
vv <- rnorm(t, 0, VV)
WW <- 1
ww <- rnorm(t, 0, WW)

#FF <- 0.95
#GG <- 1

x0 <- rnorm(1, 0, 1)

# x[i] <- x[i-1] + ww
# y[i] <- x[i] + vv
x <- c()
y <- c()
for (i in 1:t) {
  if (i==1) {
    x[i] <- x0 + ww[i]
    y[i] <- x[i] + vv[i]
  } else {
    x[i] <- x[i-1] + ww[i]
    y[i] <- x[i] + vv[i]
  }
}

# plot(1:t, x, type="o", pch=20, ylim=c(min(x,y),max(x,y)))
# lines(1:t, y, type="o", pch=20, col="red")
# legend("topright", legend=c("True X", "Observed Y", "KF", "SMC"), col=c("Black", "Red", "Blue", "Green"), lty=1)

library(dlm)
mod <- dlm(FF=1, V=VV, GG=1, W=WW, m0=0, C0=1)
filter <- dlmFilter(y, mod=mod)
# lines(1:t, filter$m[-1], type="o", pch=20, col="blue")

smc <- function(n_particles, y, ww, vv) {
  t <- length(y)

  x <- matrix(nrow=n_particles, ncol=t)
  x[,1] <- rnorm(n_particles, y[1], vv)

  logw <- rep(0,n_particles)
  int_llik <- 0
  for (j in 1:n_particles) {
    #likelihood: y1|x1 ~ N(x1,vv), prior: x1 ~ N(x0,1), proposal: q(x1) ~ N(y1,1)
    logw[j] <- dnorm(y[1], x[j,1], vv, log=T) + dnorm(x[j,1], x0, 1, log=T) - dnorm(x[j,1], y[1], vv, log=T)
  }
  int_llik <- int_llik + logSumExp(logw) - log(n_particles)
  w.norm <- exp(logw - logSumExp(logw))
  x_resample <- matrix(nrow=n_particles, ncol=t)
  if (n_particles==1) {
    x_resample[,1] <- x[,1]
  } else {
    x_resample[,1] <- sample(x[,1], size=n_particles, replace=T, prob=w.norm)
  }
  for (i in 2:t) {
    x[,i] <- rnorm(n_particles, y[i], vv)
    for (j in 1:n_particles) {
      #likelihood: yi|xi~N(xi,vv), xi|xi-1~N(xi-1,ww), q(xi)~N(yi,1)
      w[j] <- dnorm(y[i], x[j,i], vv, log=T) + dnorm(x[j,i], x_resample[j,i-1], ww, log=T) - dnorm(x[j,i], y[i], vv, log=T)
    }
    int_llik <- int_llik + logSumExp(logw) - log(n_particles)
    w.norm <- exp(logw - logSumExp(logw))
    if (n_particles==1) {
      x_resample[,i] <- x[,i]
    } else {
      x_resample[,i] <- sample(x[,i], size=n_particles, replace=T, prob=w.norm)
    }
  }
  return(int_llik)
}


# lines(1:t, apply(x_resample,2,mean),type="o",pch=20,col="green")

#find smc package and compare likelihoods
#implement into mh
#infer one parameter (FF,GG=1, infer variance of unobserved process)

library(nimble)
code <- nimble::nimbleCode({
  x0 ~ dnorm(0, var = 1)
  x[1] ~ dnorm(x0, var = 1)
  y[1] ~ dnorm(x[1], var = 0.5)
  for (t in 2:10) {
    x[t] ~ dnorm(x[t-1], var = 1)
    y[t] ~ dnorm(x[t], var = 0.5)
  }
})
model <- nimble::nimbleModel(code=code, data=list(y=y), inits=list(x0=x0, x=x))

library(nimbleSMC)
BF <- nimbleSMC::buildBootstrapFilter(model=model, nodes="x", control=list(thresh=1, saveAll=TRUE, smoothing=FALSE))

# N <- 100
# BF$run(N)
# smc(N,y=y,ww=WW,vv=VV)
#
# bf_mean <- 0
# smc_mean <- 0
# time <- Sys.time()
# for (i in 1:10) {
#   bf_mean <- bf_mean + BF$run(N)
#   smc_mean <- smc_mean + smc(N,y=y,ww=WW,vv=VV)
# }
# bf_mean <- bf_mean/10
# smc_mean <- smc_mean/10
# time <- Sys.time() - time
# c(bf_mean, smc_mean, time)

mh_smc <- function(iter, n_particles, y, ww0, vv) {
  time <- as.numeric(Sys.time())
  output <- list("var_obs"=0, "acceptance"=0, "runtime"=0)

  count <- 0

  s <- 0

  ww_old <- ww0
  prior_old <- invgamma::dinvgamma(ww_old, shape=2, scale=1,log=T)
  llik_old <- smc(n_particles=n_particles, y=y, ww=ww_old, vv=vv)

  for (i in 1:iter) {
    w <- rnorm(1,0,1)
    ww_new <- ww_old + exp(s)*w
    if (ww_new < 0) {
      ww_new <- -ww_new
    }

    prior_new <- invgamma::dinvgamma(ww_new, shape=2, scale=1, log=T)
    llik_new <- smc(n_particles=n_particles, y=y, ww=ww_new, vv=vv)

    logr <- prior_new + llik_new - prior_old - llik_old
    loga <- min(0, logr)
    a <- exp(loga)
    s <- s + (a - 0.5) / i #update s (adaptive MCMC)

    u <- runif(1, 0, 1)
    if (u <= a) {
      ww_old <- ww_new
      prior_old <- prior_new
      llik_old <- llik_new
      count <- count + 1
    }
    output$var_obs[i] <- ww_old
  }

  output$acceptance <- count/iter
  output$runtime <- as.numeric(Sys.time()) - time
  return(output)
}

# set.seed(1)
# test <- mh_smc(iter=100000, n_particles=1000, y=y, ww0=1, vv=0.5)
# test$acceptance
# test$runtime
# plot(test$var_obs, type="l")
# plot(dplyr::cummean(test$var_obs),type="l")
