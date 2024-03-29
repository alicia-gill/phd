% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propose_epi.R
\name{propose_epi}
\alias{propose_epi}
\title{Propose Epidemic}
\usage{
propose_epi(birth_rate, death_rate, stop_time, ptree, ...)
}
\arguments{
\item{birth_rate}{birth rate of the epidemic.}

\item{death_rate}{death rate of the epidemic.}

\item{stop_time}{number of days the epidemic simulation was run for.}

\item{ptree}{object of class phylo.}

\item{...}{allows for the use of lapply to run it multiple times.}
}
\value{
log-likelihood.
}
\description{
Simulates an epidemic with a given birth rate, death rate and stop time and returns the log-likelihood corresponding to that epidemic.
}
\examples{
propose_epi(birth_rate = 0.15, death_rate = 0.1, stop_time = 50, ptree = sample_tree)
}
