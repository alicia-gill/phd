% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sir.r
\name{sir}
\alias{sir}
\title{Sample Importance Resample}
\usage{
sir(
  n_particles,
  birth_rate,
  death_rate,
  noisy_prevalence,
  proportion_obs,
  genetic_data,
  plot = F
)
}
\arguments{
\item{n_particles}{number of particles used in the sampling.}

\item{birth_rate}{birth rate of the epidemic.}

\item{death_rate}{death rate of the epidemic.}

\item{noisy_prevalence}{data frame of observed prevalence per day.}

\item{proportion_obs}{proportion of cases observed.}

\item{genetic_data}{data frame of day, number of lineages and number of coalescences.}

\item{plot}{logical; if TRUE, then plots a trajectory according to weight.}
}
\value{
log-likelihood
}
\description{
Implements an SIR algorithm and returns an approximated log-likelihood.
}
\examples{
sir(n_particles = 100, birth_rate = 0.2, death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
}
