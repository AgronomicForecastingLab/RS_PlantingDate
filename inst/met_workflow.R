require(ncdf4) # package for netcdf manipulation
require(rgdal) # package for geospatial analysis
require(ncmeta)
require(stars)
require(units)
require(ecmwfr)
require(ggplot2)
require(RNetCDF)
require(raster)
require(rasterVis)
require(maptools)
require(maps)

# Open the NetCDF so we can extract met data by variable:
met_nc_data <- nc_open('inst/data/met_data.nc')

met_nc_fname <- 'inst/data/met_data.nc'
# "mx2t": max. air temperature
met_nc_mx2t.b <- brick(met_nc_fname, varname="mx2t")
# "mn2t": min. air temperature
met_nc_mn2t.b <- brick(met_nc_fname, varname="mn2t")

# This is called `siteLoc` for whatever reason:
load('inst/data/finalSiteLocations.Rdata')

n <- nrow(siteLoc) # 280

site_df <- data.frame(matrix(NA, nrow = 366, ncol = 1))

# Extract site-specific timeseries data from the raster brick 
for (i in 1:n) {
  t_lat <- round(siteLoc[i,]$Latitude, digits=2)
  t_lon <- round(siteLoc[i,]$Longitude, digits=2)

  # Extract "mx2t" data for the site
  site_series <- extract(met_nc_mx2t.b, SpatialPoints(cbind(t_lon, t_lat)), method='simple')
  t_mx2t_df <- data.frame(mx2t=t(site_series))
  agg_mx2t_df <- aggregate_met(site_series, t_df)
  
  # Call aggregation function, store returned list ("mx2t")
  site_df$mx2t <- agg_mx2t_df
  
  # Convert units from K to C
  site_df$mx2t <- site_df$mx2t - 273.15
  
  # Extract "mn2t" data for the site
  t_mn2t_df <- data.frame(mn2t=t(site_series))
  agg_mn2t_df <- aggregate_met(site_series, t_df)
  
  # Call aggregation function, store returned list ("mn2t")
  site_df$mn2t <- agg_mn2t_df
  
  # Convert units from K to C
  site_df$mn2t <- site_df$mn2t - 273.15
}

# Remove the NA column from `site_df`
site_df <- site_df[,-c(1)]

# Now that we have the avg. max. and min. daily temps, we can find GDD
# GDD = ((Daily Max Temp °C + Daily Min Temp °C)/2) - 10 °C.
# GDD = (mx2t + mn2t)/2 – 10.

gdd_list <- c()
for (i in 1:nrow(site_df)) {
  max_temp <- site_df[i,]$mx2t
  min_temp <- site_df[i,]$mn2t
  
  gdd <- ((max_temp + min_temp)/2) - 10
  gdd_list <- c(gdd_list, gdd)
}
site_df$GDD <- gdd_list
  
  

  