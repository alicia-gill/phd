#' Systematic Resampling
#'
#' Implements systematic resampling.
#'
#' @param n_particles number of particles used in the sampling.
#' @param norm_weights normalised weights of the particles.
#'
#' @return resampled particles
#' @export
#'
#' @examples
#' systematic_sample <- function(n_particles = 1000, norm_weights = rep(1/1000, 1000))
systematic_sample <- function(n_particles, norm_weights) {
  u1 <- runif(1, 0, 1/n_particles)
  u <- (0:(n_particles-1))/n_particles + u1
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
