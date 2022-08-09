t <- 10

V <- 0.5
v <- rnorm(t, 0, V)
W <- 1
w <- rnorm(t, 0, W)

#FF <- 0.95
#GG <- 1

x0 <- rnorm(1, 0, 1)

# x[i] <- x[i-1] + w
# y[i] <- x[i] + v
x <- c()
y <- c()
for (i in 1:t) {
  if (i==1) {
    x[i] <- x0 + w[i]
    y[i] <- x[i] + v[i]
  } else {
    x[i] <- x[i-1] + w[i]
    y[i] <- x[i] + v[i]
  }
}

plot(1:t, x, type="o", pch=20, ylim=c(min(x,y),max(x,y)))
lines(1:t, y, type="o", pch=20, col="red")
legend("topleft", legend=c("True X", "Observed Y", "KF", "SMC"), col=c("Black", "Red", "Blue", "Green"), lty=1)

library(dlm)
mod <- dlm(FF=1, V=V, GG=1, W=W, m0=0, C0=1)
filter <- dlmFilter(y, mod=mod)
lines(1:t, filter$m[-1], type="o", pch=20, col="blue")

smc <- function(n_particles, y) {
  t <- length(y)
  x <- matrix(nrow=n_particles, ncol=t)
  x[,1] <- rnorm(n_particles, y[1], 1)
  w <- rep(0,n_particles)
  int_llik <- 0
  for (j in 1:n_particles) {
    #likelihood: y1|x1 ~ N(FF*x1,1), prior: x1 ~ N(GG*x0,1), proposal: q(x1) ~ N(GG*y1,1)
    w[j] <- dnorm(y[1], x[j,1], 1) * dnorm(x[j,1], x0, 1) / dnorm(x[j,1], y[1], 1)
  }
  int_llik <- int_llik + log(mean(w))
  w.norm <- w/sum(w)
  x_resample <- matrix(nrow=n_particles, ncol=t)
  x_resample[,1] <- sample(x[,1], size=n_particles, replace=T, prob=w.norm)
  for (i in 2:t) {
    x[,i] <- rnorm(n_particles, y[i], 1)
    for (j in 1:n_particles) {
      #likelihood: yi|xi~N(xi,1), xi|xi-1~N(xi-1,1), q(xi)~N(GG*yi,1)
      w[j] <- dnorm(y[i], x[j,i], 1) * dnorm(x[j,i], x_resample[j,i-1], 1) / dnorm(x[j,i], y[i], 1)
    }
    int_llik <- int_llik + log(mean(w))
    w.norm <- w/sum(w)
    x_resample[,i] <- sample(x[,i], size=n_particles, replace=T, prob=w.norm)
  }
  return(int_llik)
}

lines(1:t, apply(X_resample,2,mean),type="o",pch=20,col="green")

#find smc package and compare likelihoods
#implement into mh
#infer one parameter (FF,GG=1, infer variance of unobserved process)

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

BF <- nimbleSMC::buildBootstrapFilter(model=model, nodes="x", control=list(thresh=1, saveAll=TRUE, smoothing=FALSE))
BF$run(N)
smc(N,y=y)

bf_mean <- 0
smc_mean <- 0
time <- Sys.time()
for (i in 1:10) {
  bf_mean <- bf_mean + BF$run(N)
  smc_mean <- smc_mean + smc(N,y=y)
}
bf_mean <- bf_mean/10
smc_mean <- smc_mean/10
time <- Sys.time() - time
c(bf_mean, smc_mean, time)

