% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noisy_prevalence.R
\name{noisy_prevalence}
\alias{noisy_prevalence}
\title{Noisy Prevalence}
\usage{
noisy_prevalence(prevalence, proportion_obs)
}
\arguments{
\item{prevalence}{data frame of prevalence per day.}

\item{proportion_obs}{proportion of cases observed.}
}
\value{
data frame of observed prevalence per day.
}
\description{
Generates a partially observed prevalence.
}
\examples{
noisy_prevalence(prevalence = prev, proportion_obs = 0.2)
}
