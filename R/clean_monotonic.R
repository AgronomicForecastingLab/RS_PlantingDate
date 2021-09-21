#' Removes all non-monotonically increasing or decreasing data points for each
#' site (data with the same ID) in the Beck's Hybrids dataset.
#'
#' **Will be further generalized to handle a range of datasets (with some
#' formatting requirements)
#'
#' @param orig_df A data frame created from the Beck's dataset.
#' @param ID_list A list of all unique IDs within the Beck's dataset.
#' @return A data frame of Beck's data w/o divergent data points.
#' @export
clean_monotonic <- function(orig_df) {
  
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
    
    # Find the highest NDVI value for this point (as the median data point).
    maxNDVI = max(this$NDVI)
    first = min(this$orig_row)
    last = max(this$orig_row)
    max_index = this$orig_row[which(this$NDVI == maxNDVI)]
    
    if (length(max_index) == 0) {
      next
    }
    
    if (length(max_index) > 1) {
      print(i)
      break
    }
    
    # Start a vector of original df rows to remove 
    rem_rows <- c()
    if (nrow(bad) > 0) {
      rem_rows <- c(bad$orig_row)
    }
    
    # Iterate through all data points.
    for (j in first:last) {
      
      # Check to see if we've already removed this point.
      if (!(j %in% this$orig_row)) next
      
      # If this point is the max NDVI point, skip. 
      if (j == max_index) next
      
      # If less than max, we check for increasing for all points up through max.
      if (j < max_index) {
        inds = (this %>% filter(orig_row > j, orig_row <= max_index))$orig_row
        diffs <- (this %>% filter(orig_row %in% inds))$NDVI - (this %>% filter(orig_row == j))$NDVI
        if (any(diffs < 0)){
          ## Only remove rows with a decrease of greater than 0.1
          rem_rows <- c(rem_rows,inds[which(diffs < -0.1)])
          this = this %>% filter(!(orig_row %in% rem_rows))
        }
      } else {
        # If this point is after max_index, we check for monotonic decreasing.
        inds = (this %>% filter(orig_row > j))$orig_row
        diffs <- (this %>% filter(orig_row %in% inds))$NDVI - (this %>% filter(orig_row == j))$NDVI
        if (any(diffs > 0)){
          rem_rows <- c(rem_rows,inds[which(diffs > 0.1)])
          this = this %>% filter(!(orig_row %in% rem_rows))
        }
      }
    }
    
    if (i %% 30 == 0){
       before = orig_df %>% filter(ID == ID_list[i]) 
       p = ggplot() + 
         geom_point(data = this, aes(x = Date, y = NDVI), col = 'black') +
         geom_point(data = before, aes(x = Date, y = NDVI), col = 'red') +
         geom_line(data = this, aes(x = Date, y = NDVI), col = 'black') +
         geom_line(data = before, aes(x = Date, y = NDVI), col = 'red') 
       print(p)
       Sys.sleep(5)
    }
    
    # Remove the rem_rows indices from the original dataset.
    orig_df = orig_df %>% filter(!orig_row %in% rem_rows)
  }
  
  return(orig_df %>% dplyr::select(-orig_row))
}
