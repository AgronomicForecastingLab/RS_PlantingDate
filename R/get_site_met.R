#' Extracts ERA5 met. data for the site with "site_id" in the provided 
#' "sites_df" dataframe from the met. data dataframe "met_df". 
#'
#' @param sites_df
#' @param site_id
#' @param mx2t_brick
#' @param mn2t_brick
#' @return A data frame with a column from the aggregated time series data.
#' @export
get_site_met <- function(sites_df, site_id, mx2t_brick, mn2t_brick) {
  # Establish a data frame for met. data at the site corresponding to `site_id`:
  site_df <- data.frame(matrix(NA, nrow = 366, ncol = 1))
  site_row <- sites_df[sites_df$ID == site_id,][1,] # Find site with ID `site_id`

  # Retrieve site lon., lat. coordinates (rounded to 2 digits).
  t_lat <- round(site_row$Latitude, digits=2)
  t_lon <- round(site_row$Longitude, digits=2)
  
  # Extract and aggregate "mx2t" data for the site.
  site_series <- extract(mx2t_brick, SpatialPoints(cbind(t_lon, t_lat)), method='simple')
  t_mx2t_df <- data.frame(mx2t=t(site_series))
  agg_mx2t_df <- aggregate_met(site_series, t_mx2t_df)
  site_df$mx2t <- agg_mx2t_df
  site_df$mx2t <- site_df$mx2t - 273.15 # Convert units from to K to C
  
  # Extract and aggregate "mn2t" data for the site.
  site_series <- extract(mn2t_brick, SpatialPoints(cbind(t_lon, t_lat)), method='simple')
  t_mn2t_df <- data.frame(mn2t=t(site_series))
  agg_mn2t_df <- aggregate_met(site_series, t_mn2t_df)
  site_df$mn2t <- agg_mn2t_df
  site_df$mn2t <- site_df$mn2t - 273.15 # Convert units from to K to C
  
  site_df <- site_df[,-c(1)]
  return(site_df)
}