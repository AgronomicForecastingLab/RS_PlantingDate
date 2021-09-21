#' Filters data using a weighted least-squares linear regression model and a 
#' moving window (aims to minimize squared error).
#'
#' @param orig_df A data frame created from the Beck's dataset.
#' @param ID_list A list of all unique IDs within the Beck's dataset.
#' @return A data frame of Beck's data w/o divergent data points.
#' @export
clean_WLS <- function(orig_df) {
  
  # Add original row numbers as a column.
  orig_df = orig_df %>% 
    rownames_to_column('orig_row') %>%
    mutate(orig_row = as.numeric(orig_row))
  
  ID_list = unique(orig_df$ID)
  pb = txtProgressBar(0, 1, style=3)
  
  # Iterate through all sites for which data was collected (data with same ID).
  for (i in 1:length(ID_list)) {
    setTxtProgressBar(pb, i/length(ID_list))
    
    # Data frame with data for ID `ID_list[i]`.
    this = orig_df %>% filter(ID == ID_list[i]) 
    
    # Remove any values before April with NDVI values greater than 0.3.
    # Here I'm using the range given for bare soil in Mzid et al. (2021)
    bad <- this %>% filter(Date < as.Date('2017-04-01') & NDVI > 0.3)
    if (nrow(bad) > 0){
      this <- this %>% filter(!(orig_row %in% bad$orig_row))
    }
    
    # Convert from formatted date to Julian days. 
    this$Date <- difftime(this$Date, as.Date("2016-12-31","%Y-%m-%d"), units="days")
    
    # predictor variable: date
    # explanatory variable: NDVI
    model <- lm(NDVI ~ Date, data = this)
    summary(model)
    
    # Plot:
    plot(fitted(model), resid(model), xlab='Fitted Values', ylab='Residuals')
    abline(0,0)
    
    # Calculate weight
    wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
    wls_model <- lm(NDVI ~ Date, data = this, weights=wt)
    
    summary(wls_model)
    plot(fitted(wls_model), resid(wls_model), xlab='Fitted Values', ylab='Residuals')
    abline(0,0)
    
    residuals <- residuals(wls_model)
    
    # What I don't understand here is how the "moving window" idea works ^^
    # I.e., is the WLS Model recalculated at every data point?
    
    # Iterate through all data points.
    for (j in first:last) { 
      # Remove points with residuals greater than some value (0.3? 0.4?)
      resid <- residuals[[j]]
      
    }
    
    p = ggplot() + 
      geom_point(data = this, aes(x = Date, y = NDVI), col = 'black') +
      geom_line(data = this, aes(x = Date, y = NDVI), col = 'black')
    
   # orig_df = orig_df %>% filter(!orig_row %in% rem_rows)
  }
  
}
  