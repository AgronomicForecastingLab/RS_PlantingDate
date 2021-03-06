% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_dates.R
\name{aggregate_dates}
\alias{aggregate_dates}
\title{Applies a uniform 3x3 smoothing window across the input matrix
(excluding pixels equal to NA or -999).}
\usage{
aggregate_dates(x, y, n)
}
\arguments{
\item{x}{x-dimension (width) of the NDVI images.}

\item{y}{y-dimension(length) of the NDVI images.}

\item{n}{number of observed dates.}
}
\value{
A x-y-n matrix of NDVI data.
}
\description{
A smaller smoothing window (e.g., 2x2, 2x3) is applied to edge and corner pixels.
}
