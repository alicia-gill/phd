% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sir_adaptive.R
\name{sir_adaptive}
\alias{sir_adaptive}
\title{Sample Importance Resample - Adaptive}
\usage{
sir_adaptive(
  n_particles,
  birth_rate,
  death_rate,
  noisy_prevalence,
  proportion_obs,
  genetic_data,
  ess_threshold = n_particles/2,
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

\item{ess_threshold}{threshold of ESS below which triggers resampling.}

\item{plot}{logical; if TRUE, then plots a trajectory according to weight.}
}
\value{
log-likelihood
}
\description{
Implements an adaptive SIR algorithm and returns an approximated log-likelihood. Proposals are centred on the smoothed prevalence.
}
\examples{
sir_adaptive(n_particles = 100, birth_rate = 0.2, death_rate = 0.1, noisy_prevalence = noisy_prev, proportion_obs = 0.2, genetic_data = gen_data)
}
