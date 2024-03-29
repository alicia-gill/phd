% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtskellam_smc.R
\name{rtskellam_smc}
\alias{rtskellam_smc}
\title{Truncated Skellam distribution}
\usage{
rtskellam_smc(n, old_x, birth_rate, death_rate, trunc = -old_x)
}
\arguments{
\item{n}{number of observations.}

\item{old_x}{old prevalence.}

\item{birth_rate}{birth rate of the epidemic.}

\item{death_rate}{death rate of the epidemic.}

\item{trunc}{truncation point; -old_x by default to prevent the epidemic dying out.}
}
\value{
log-likelihood.
}
\description{
Random generation truncated at specified location for the Skellam distribution
}
\examples{
rtskellam_smc(n = 100, old_x = x_resample, birth_rate = b_sample, death_rate = 0.1)
}
