#' Aggregates ERA5 met. data from an hourly to daily time step for the variable
#' specified. 
#'
#' @param series
#' @param df 
#' @return A data frame with a column from the aggregated time series data.
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
  df <- subset(df, select=-c(rownumber))

  return(df)
}






