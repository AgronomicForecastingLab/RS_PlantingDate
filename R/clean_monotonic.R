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
clean_monotonic <- function(orig_df, ID_list) {
  # Iterate through all sites for which data was collected (data with same ID).
  for (i in 1:length(ID_list)) {
    this = orig_df %>% rownames_to_column('orig_row') %>% filter(ID == ID_list[i])
    
    # Find the highest NDVI value for this point (as the median data point).
    maxNDVI = max(this$NDVI)
    mn_index = which(this$NDVI == maxNDVI)
    
    if (length(mn_index) == 0) {
      next
    }
    
    incr = this$NDVI[1:mn_index]
    decr = this$NDVI[mn_index:nrow(this)]
    
    rem_rows <- c()
    # Check that from planting to peak, data is monotonically increasing.
    if (!all(diff(incr) >= 0)) {
      # Add divergent data points to list of rows to be removed from `temp`.
      rem_rows <- which(diff(incr) < 0)
    }
    
    # Check that from peak to maturity, data is monotonically decreasing.
    if (!all(diff(decr) <= 0)) {
      # Add divergent data points to list of rows to be removed from `temp`.
      rem_rows <- c(rem_rows, which(diff(decr) > 0))
    }
    
    # Convert rem_rows list from `this` rownames to `temp` rownames.
    rem_exes <- c()
    for (i in 1:length(rem_rows)) {
      rem_exes <- c(rem_exes, this[rem_rows[i], ]$X)
    }
    
    # Filter divergent data points from `temp`.
    orig_df = orig_df %>% filter(!X %in% rem_exes)
  }
  return(orig_df)
}
