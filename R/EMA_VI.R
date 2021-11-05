#' An Exponential Moving Average (EMA) is a moving average filter that applies
#' exponentially decreasing weighting factors to time-series datum within the
#' moving window.
#'
#' @param ndvi NDVI time-series on which to run the EMA
#' @param n Length of EMA moving window (# of time-series datum in the window)
#' @return An EMA time-series
#' @export
EMA_VI <- function(ndvi, n) {
  ema <- c()
  ema[1:(n-1)] <- NA
  ema[n]<- mean(ndvi[1:n])
  beta <- 2/(n+1)
  for (i in (n+1):length(ndvi)){
    ema[i]<-beta * ndvi[i] + 
      (1-beta) * ema[i-1]
  }
  return(ema)
}