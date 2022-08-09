t <- 100

V <- 1
v <- rnorm(t, 0, V)
W <- 1
w <- rnorm(t, 0, W)

FF <- 0.95
GG <- 1

m0 <- 0
C0 <- 1
x0 <- rnorm(1, m0, C0)

x <- c()
y <- c()
for (i in 1:t) {
  if (i==1) {
    x[i] <- GG*x0 + w[i]
    y[i] <- FF*x[i] + v[i]
  } else {
    x[i] <- GG*x[i-1] + w[i]
    y[i] <- FF*x[i] + v[i]
  }
}

plot(1:t, x, type="o", pch=20, ylim=c(min(x,y),max(x,y)))
lines(1:t, y, type="o", pch=20, col="red")
legend("topleft", legend=c("True X", "Observed Y", "KF", "SMC"), col=c("Black", "Red", "Blue", "Green"), lty=1)

library(dlm)
mod <- dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)
filter <- dlmFilter(y, mod=mod)
lines(1:t, filter$m[-1], type="o", pch=20, col="blue")

N <- 1000
X <- matrix(nrow=N, ncol=t)
X[,1] <- rnorm(N, FF*y[1], 1)
w <- rep(0,N)
int_llik <- 0
for (j in 1:N) {
  #likelihood: y1|x1 ~ N(FF*x1,1), prior: x1 ~ N(GG*m0,1), proposal: q(x1) ~ N(GG*y1,1)
  w[j] <- dnorm(y[1], FF*X[j,1], 1) * dnorm(X[j,1], GG*m0, 1) / dnorm(X[j,1], GG*y[1], 1)
}
int_llik <- int_llik + log(mean(w))
w.norm <- w/sum(w)
X_resample <- matrix(nrow=N, ncol=t)
X_resample[,1] <- sample(X[,1], size=N, replace=T, prob=w.norm)
for (i in 2:t) {
  X[,i] <- rnorm(N, FF*y[i], 1)
  for (j in 1:N) {
    #likelihood: yi|xi~N(FF*xi,1), xi|xi-1~N(GG*xi-1,1), q(xi)~N(GG*yi,1)
    w[j] <- dnorm(y[i], FF*X[j,i], 1) * dnorm(X[j,i], GG*X_resample[j,i-1], 1) / dnorm(X[j,i], GG*y[i], 1)
  }
  int_llik <- int_llik + log(mean(w))
  w.norm <- w/sum(w)
  X_resample[,i] <- sample(X[,i], size=N, replace=T, prob=w.norm)
}
print(int_llik)

lines(1:t, apply(X_resample,2,mean),type="o",pch=20,col="green")

#find smc package and compare likelihoods
#implement into mh
#infer one parameter (FF,GG=1, infer variance of unobserved process)
