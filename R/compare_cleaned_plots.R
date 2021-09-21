#' Returns a ggplot comparing uncleaned data plots to cleaned data plots.
#' 
#' @param IDs IDs of data within the data frame to compare.
#' @param n Size of random sample of data to include in comparison. 
#' @param before The data frame before data cleaning.
#' @param after The data frame after data cleaning. 
#' @return A ggplot comparing uncleaned and cleaned data plots. 
#' @export
compare_cleaned_plots <- function(IDs, n, before, after) {
  ids = sample(unique(before$ID), n, replace=F)
  before = before %>% filter(ID %in% ids)
  before$check = 'before'
  
  after = after %>% filter(ID %in% ids)
  after$check = 'after'
  check = rbind(before,after) %>% mutate(Date = as.Date(Date))
  
  comp_plot <- ggplot(check)+
    geom_line(aes(x = Date, y = NDVI, color = check), alpha = 0.5) + 
    geom_point(aes(x = Date, y = NDVI, color = check)) + 
    facet_wrap(~ID, nrow = 4)
  
  return(comp_plot)
}