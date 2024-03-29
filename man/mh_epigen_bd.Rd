% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mh_epigen_bd.R
\name{mh_epigen_bd}
\alias{mh_epigen_bd}
\title{Metropolis--Hastings using Epidemic and Genetic data}
\usage{
mh_epigen_bd(
  iter,
  birth_rate_0,
  max_birth_rate = 100,
  death_rate_0,
  max_death_rate = 100,
  prevalence,
  ptree,
  print = F
)
}
\arguments{
\item{iter}{number of iterations to run the algorithm for.}

\item{birth_rate_0}{initial value of the birth rate.}

\item{max_birth_rate}{maximum value for the birth rate. 100 by default.}

\item{death_rate_0}{initial value of the death rate.}

\item{max_death_rate}{maximum value for the death rate. 100 by default.}

\item{prevalence}{data frame of prevalence per day.}

\item{ptree}{object of class phylo.}

\item{print}{logical; if TRUE, prints percentage of the way through the chain.}
}
\value{
list containing: birth rate, death rate, acceptance rate, run time in seconds
}
\description{
Runs a Metropolis--Hastings algorithm to find the birth and death rate of the epidemic.
}
\examples{
mh_epigen(iter = 100000, birth_rate_0 = 0.1, death_rate = 0.1, prevalence = prev, ptree = sample_tree)
}
