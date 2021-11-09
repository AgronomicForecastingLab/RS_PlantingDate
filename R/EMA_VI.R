#' An Exponential Moving Average (EMA) is a moving average filter that applies
#' exponentially decreasing weighting factors to time-series datum within the
#' moving window.
#'
#' @param series NDVI time-series on which to run the EMA
#' @param n Length of EMA moving window (# of time-series datum in the window)
#' @return An EMA time-series
#' @export
EMA_VI <- function(series, n) {
  ema <- c()
  ema[1:(n - 1)] <- NA
  ema[n] <- mean(series[1:n], na.rm=TRUE)
  beta <- 2.0/(n + 1)
  for (i in (n + 1):length(series)){
    ema[i] <- (beta * series[i]) + ((1-beta) * ema[i - 1])
  }
  return(ema)
}