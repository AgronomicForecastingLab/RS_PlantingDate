#' Smooths NDVI time series using splines 
#'
#' @param orig_df A data frame created from the Beck's dataset.
#' @param ID_list A list of all unique IDs within the Beck's dataset.
#' @return A data frame of Beck's data w/o divergent data points.
#' @export
clean_spline <- function(orig_df) {
  
  require(splines)
  
  # Add original row numbers as a column
  orig_df = orig_df %>% 
    rownames_to_column('orig_row') %>%
    mutate(orig_row = as.numeric(orig_row))
  
  new_df = orig_df[1,]
  
  ID_list = unique(orig_df$ID)
  pb = txtProgressBar(0, 1, style=3)
  
  # Iterate through all sites for which data was collected (data with same ID).
  for (i in 1:length(ID_list)) {
    setTxtProgressBar(pb, i/length(ID_list))
    
    # Data frame with data for ID `ID_list[i]`.
    this = orig_df %>% dplyr::filter(ID == ID_list[i]) 
    
    if (nrow(this) < 4) next
    
    # Convert from formatted date to Julian days. 
    this$DOY <- as.numeric(difftime(this$Date, as.Date("2016-12-31","%Y-%m-%d"), units="days"))
    
    # Perform cubic spline interpolation
    fit1<-smooth.spline(this$DOY,this$NDVI)

    # Adjust for peak NDVI value 
    max.ind = which.max(this$NDVI)
    max = this$NDVI[max.ind]
    this$NDVI = fit1$y
    this$NDVI[max.ind] = max
    
    new_df = rbind(new_df, 
                   this %>% dplyr::select(-DOY))
  }
  
  new_df = new_df[-1,]
  return(new_df %>% dplyr::select(-orig_row))

}
