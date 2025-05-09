% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preprocessing.R
\name{PlotLibrarySaturation}
\alias{PlotLibrarySaturation}
\title{Plot Library Saturation}
\usage{
PlotLibrarySaturation(citeseqCountDir, metricsFile = NULL)
}
\arguments{
\item{citeseqCountDir}{The root of the Cite-seq-Count output folder, which should contain umi_count and read_count folders.}

\item{metricsFile}{If provided, summary metrics will be written to this file.}
}
\value{
The overall saturation for this library
}
\description{
Create a plot of the library saturation per cell
}
