% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFF_Demux.R
\name{ParameterScan}
\alias{ParameterScan}
\title{Perform a parameter scan with BFF parameters}
\usage{
ParameterScan(lognormedcounts)
}
\arguments{
\item{lognormedcounts}{Count matrix after BQN and log transformation}
}
\description{
Prints plots that demonstrate the impact of BFF parameters.
Parameters alpha_c, beta_c, and delta_c can take values between 0 and 1.
Reasonable parameter values are less than 0.25.  The scan is over parameter
values (0.05, 0.1, 0.15, 0.2, 0.25, 0.5) for each of the parameters.
}
