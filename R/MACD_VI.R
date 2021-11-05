#' The Moving Average Convergence/Divergence (MACD) is an indicator composed of
#' 3 time series (the MACD time series, the "signal"/"average" series, and the
#' "divergence" series)
#' MACD(t) = EMA(v(t),a) - EMA(v(t),b)
#'
#' @param ndvi NDVI time-series on which to run the EMA
#' @param S Length of the short ("fast") EMA moving window
#' @param L Length of the long ("slow") EMA moving window
#' @param K Length of the "signal"/"average" EMA moving window 
#' @return A MACD time-series
#' @export
MACD_VI <- function(ndvi, S, L, K) {
  MACD <- EMA_VI(ndvi,S) - EMA_VI(ndvi,L)
  signal <- EMA_VI(MACD,K)
  output <- cbind(MACD,signal)
  colnames(output) <- c("MACD","signal")
  return(output)
}