% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtskellam.R
\name{dtskellam}
\alias{dtskellam}
\title{Density function for the truncated Skellam distribution}
\usage{
dtskellam(x, a = -1e+150, b = 1e+300, mu1, mu2, log = T)
}
\arguments{
\item{x}{vector of quantities.}

\item{a, b}{lower and upper truncation points (a < x <= b).}

\item{mu1, mu2}{positive valued parameters.}

\item{log}{logical; if TRUE, probabilities are given as log(p).}
}
\value{
density.
}
\description{
Density function for the truncated Skellam distribution.
}
\examples{
dtskellam(x = 10, a=0, mu1 = 20, mu2 = 10)
}
