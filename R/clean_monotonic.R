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
  
  orig_df = orig_df %>% 
    rownames_to_column('orig_row') %>%
    mutate(orig_row = as.numeric(orig_row))
  ID_list = unique(orig_df$ID)
  pb = txtProgressBar(0, 1, style=3)
  
  # Iterate through all sites for which data was collected (data with same ID).
  for (i in 1:length(ID_list)) {
    setTxtProgressBar(pb, i/length(ID_list))
    
    this = orig_df %>% filter(ID == ID_list[i]) 
    
    # Are there any values before April that are not in the bare soil range?
    # Here I'm using the range given for bare soil in Mzid et al. (2021)
    bad <- this %>% filter(Date < as.Date('2017-04-01'), NDVI > 0.3)
    if (nrow(bad) > 0){
      reduced <- this %>% filter(!(orig_row %in% bad$orig_row))
    }else{
      reduced <- this
    }
    
    # Find the highest NDVI value for this point (as the median data point).
    maxNDVI = max(reduced$NDVI)
    first = min(reduced$orig_row)
    last = max(reduced$orig_row)
    mn_index = reduced$orig_row[which(reduced$NDVI == maxNDVI)]
    
    if (length(mn_index) == 0) {
      next
    }
    if (length(mn_index) > 1){
      print(i)
      break
    }
    
    # Start a vector of original df rows to remove 
    if (nrow(bad) > 0){
      rem_rows <- c(bad$orig_row)
    }else{ 
      rem_rows <- c()
    }
    
    # Iterate through all the points 
    for (j in first:last){
      
      # check to see if we have removed this point 
      if (!(j %in% reduced$orig_row)) next
      
      # if this is max, we can skip
      if (j == mn_index) next
      
      # if less than max, we check for increasing with all points up through max
      if (j < mn_index){
        inds = (reduced %>% filter(orig_row > j, orig_row <= mn_index))$orig_row
        diffs <- (reduced %>% filter(orig_row %in% inds))$NDVI - (reduced %>% filter(orig_row == j))$NDVI
        if (any(diffs < 0)){
          rem_rows <- c(rem_rows,inds[which(diffs < -0.1)])
          reduced = reduced %>% filter(!(orig_row %in% rem_rows))
        }
      }else{
        # after mn_index we check for decreasing 
        inds = (reduced %>% filter(orig_row > j))$orig_row
        diffs <- (reduced %>% filter(orig_row %in% inds))$NDVI - (reduced %>% filter(orig_row == j))$NDVI
        if (any(diffs > 0)){
          rem_rows <- c(rem_rows,inds[which(diffs > 0.1)])
          reduced = reduced %>% filter(!(orig_row %in% rem_rows))
        }
      }
    }
    
    # if (i %% 30 == 0){
    #   p = ggplot() + 
    #     geom_point(data = this, aes(x = Date, y = NDVI), col = 'black') +
    #     geom_point(data = reduced, aes(x = Date, y = NDVI), col = 'red') +
    #     geom_line(data = this, aes(x = Date, y = NDVI), col = 'black') +
    #     geom_line(data = reduced, aes(x = Date, y = NDVI), col = 'red') 
    #   print(p)
    #   Sys.sleep(5)
    # }
    
    # Now let's remove the rem_rows indices from the original dataset
    orig_df = orig_df %>% filter(!orig_row %in% rem_rows)
    
  }
  
  return(orig_df %>% dplyr::select(-orig_row))
}
