% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fns.R
\name{plot_fns}
\alias{plot_fns}
\alias{plot_birthrates}
\alias{plot_rt}
\alias{plot_prevalence}
\alias{plot_pobs}
\alias{plot_sigma}
\alias{plot_x0}
\title{A collection of functions to plot inference from MCMC chains}
\usage{
plot_birthrates(chain, q = 0.95, burn_in = 0, start_time = NA, end_time = NA)

plot_rt(
  chain,
  death_rate,
  q = 0.95,
  burn_in = 0,
  start_time = NA,
  end_time = NA
)

plot_prevalence(chain, q = 0.95, burn_in = 0, start_time = NA, end_time = NA)

plot_pobs(chain, burn_in = 0, iter = NA)

plot_sigma(chain, burn_in = 0, iter = NA)

plot_x0(chain, burn_in = 0, iter = NA)
}
\arguments{
\item{chain}{MCMC chain}

\item{q}{coverage of credible intervals; default 95\%}

\item{burn_in}{burn in}

\item{death_rate}{death rate of the epidemic}
}
\value{
plot
}
\description{
Plotting functions
}
