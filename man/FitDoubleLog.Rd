% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_double_logistic.R
\name{FitDoubleLog}
\alias{FitDoubleLog}
\title{Function to fit a double logistic function that has been adjusted from FitDoubleLogBeck function
That function did not accurately represent the equation as given in Beck et al. 
In particular, there is a subtract of 1 from the exponential part that does not exist in the original function. This resulted 
in the parameters not having meaningful values in the fit.
You can look at the arguments for the original function to understand the parameters of this function.}
\usage{
FitDoubleLog(
  x,
  t = 1:length(x),
  tout = t,
  weighting = TRUE,
  hessian = FALSE,
  plot = FALSE,
  ninit = 30,
  ...
)
}
\arguments{
\item{dlog.data}{An empty data frame for per-site predictions.}
}
\value{
A data frame of the fitted double logistic function for each site.
}
\description{
Function to fit a double logistic function that has been adjusted from FitDoubleLogBeck function
That function did not accurately represent the equation as given in Beck et al. 
In particular, there is a subtract of 1 from the exponential part that does not exist in the original function. This resulted 
in the parameters not having meaningful values in the fit.
You can look at the arguments for the original function to understand the parameters of this function.
}
