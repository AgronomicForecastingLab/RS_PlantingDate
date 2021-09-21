#' Returns a ggplot comparing uncleaned data plots to cleaned data plots.
#' 
#' @param IDs IDs of data within the data frame to compare.
#' @param n Size of random sample of data to include in comparison. 
#' @param before The data frame before data cleaning.
#' @param after The data frame after data cleaning. 
#' @return A ggplot comparing uncleaned and cleaned data plots. 
#' @export
compare_cleaned_plots <- function(ids, before, after) {
  before = before %>% filter(ID %in% ids)
  before$check = 'before'
  
  after = after %>% filter(ID %in% ids)
  after$check = 'after'
  check = rbind(before,after) %>% mutate(Date = as.Date(Date))
  
  comp_plot <- ggplot(check)+
    geom_point(aes(x = Date, y = NDVI, color = check)) + 
    geom_line(data = check %>% filter(check == 'after'), aes(x = Date, y = NDVI), linetype = 'dashed') +
    scale_color_manual(values = c('before' = 'red', 'after' = 'black'))+
    facet_wrap(~ID, nrow = 4)
  
  return(comp_plot)
}
