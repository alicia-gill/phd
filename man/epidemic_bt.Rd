% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epidemic_bt.R
\name{epidemic_bt}
\alias{epidemic_bt}
\title{Simulate epidemics with birth rate varying over time.}
\usage{
epidemic_bt(birth_rate_fn, min_r = 1, max_r = 5, death_rate, stop_time)
}
\arguments{
\item{birth_rate_fn}{function of the birth rate over time of the epidemic.}

\item{min_r}{minimum reproduction number}

\item{max_r}{maximum reproduction number}

\item{death_rate}{death rate of the epidemic.}

\item{stop_time}{number of days to run the epidemic simulation for.}
}
\value{
data frame with: Index, Parent Index, Infection Status, Infection Time, Removal Time
}
\description{
Simulates a birth-death epidemic with varying birth rate and constant death rates and 1 infected on day 0.
}
\examples{
epidemic_bt(birth_rate_fn = linear, death_rate = 0.1, stop_time = 50)
}
