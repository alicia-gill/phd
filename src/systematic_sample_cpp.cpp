#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]
//' Systematic Resampling
//'
//' Implements systematic resampling.
//'
//' @param n_particles number of particles used in the sampling.
//' @param norm_weights normalised weights of the particles.
//'
//' @return resampled particles
//' @export
//'
//' @examples
//' systematic_sample_cpp <- function(n_particles = 1000, norm_weights = rep(1/1000, 1000))
//'
//' @useDynLib phd, .registration = TRUE
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector systematic_sample_cpp(int n_particles, NumericVector norm_weights) {
  double u1 = runif(1)[0] / n_particles;
  NumericVector u(n_particles);
  for (int k = 0; k < n_particles; ++k) {
    u[k] = u1 + (double)k / (double)n_particles;
  }
  NumericVector cumulative_sum = cumsum(norm_weights);
  NumericVector resample(n_particles);
  int i = 0;
  int j = 0;
  while (i < n_particles) {
    if (u[i] < cumulative_sum[j]) {
      resample[i] = j + 1;
      i = i + 1;
    } else {
      j = j + 1;
    }
  }
  return resample;
}
