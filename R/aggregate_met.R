#' Aggregates ERA5 met. data from an hourly to daily time step for the variable
#' specified. 
#'
#' @param x x-dimension (width) of the NDVI images.
#' @param y y-dimension(length) of the NDVI images.
#' @param n number of observed dates.
#' @return A x-y-n matrix of NDVI data.
#' @export
aggregate_met <- function(series, df) {
  start <- 1
  end <- 24
  days <- length(series) # 8784
  dont_remove <- c()
  
  while (end <= days) {
    agg_var <- sum(series[,start:end]) / 24
    dont_remove <- c(dont_remove, start)
    
    df[start,] = agg_var
    
    start <- start + 24
    end <- end + 24 
  }
  
  df$rownumber = 1:nrow(df)
  df <- df[df$rownumber %in% dont_remove,]

  return(df)
}






