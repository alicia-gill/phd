systematic_sample <- function(n_particles, norm_weights) {
  u <- rep(NA, n_particles)
  u[1] <- runif(1, 0, 1/n_particles)
  for (i in 2:n_particles) {
    u[i] <- u[1] + (i-1)/n_particles
  }
  cumulative_sum <- cumsum(norm_weights)
  resample <- rep(NA, n_particles)
  i <- 1
  j <- 1
  while (i <= n_particles) {
    if (u[i] < cumulative_sum[j]) {
      resample[i] <- j
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  return(resample)
}
