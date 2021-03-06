% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_spline.R
\name{clean_spline}
\alias{clean_spline}
\title{Smooths NDVI time series using splines}
\usage{
clean_spline(orig_df)
}
\arguments{
\item{orig_df}{A data frame created from the Beck's dataset.}

\item{ID_list}{A list of all unique IDs within the Beck's dataset.}
}
\value{
A data frame of Beck's data w/o divergent data points.
}
\description{
Smooths NDVI time series using splines
}
