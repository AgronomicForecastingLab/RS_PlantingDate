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

# There are 3 dimensions in the met .nc (lat., lon., time).
lon <- ncvar_get(met_nc_data, "longitude")
lat <- ncvar_get(met_nc_data, "latitude")
t <- ncvar_get(met_nc_data, "time")
# Convert the time series units from NetCDF time to POSIXct R dates. 
t.POSIXct <- utcal.nc("hours since 1900-01-01 00:00:00 +00:00", t, type="c")

# note that you may have to play around with the transpose (the t() function)
# and flip() before the data are oriented correctly. In this example, the netcdf
# file recorded latitude on the X and longitude on the Y, so both a transpose 
# and a flip in the y direction were required.
# r <- raster(t(met.mx2t), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# r <- flip(r, direction='y')
# plot(r)

met_nc_fname <- 'inst/data/met_data.nc'
# "mx2t": max. air temperature
met_nc_mx2t.b <- brick(met_nc_fname, varname="mx2t")
# "mn2t": min. air temperature
met_nc_mn2t.b <- brick(met_nc_fname, varname="mn2t")

# Extract timeseries data from the raster brick
t_lon <- -89.06
t_lat <- 40.15

# Extract "mx2t" data for the site
site_series <- extract(met_nc_mx2t.b, SpatialPoints(cbind(t_lon, t_lat)), method='simple')
t_df <- data.frame(mx2t=t(site_series))

# Aggregate "mx2t" data (hourly -> daily timestep)
start <- 1
end <- 24
days <- length(site_series) # 8784
dont_remove <- c()

while (end <= days) {
  print(end)
  agg_var <- sum(site_series[,start:end]) / 24
  dont_remove <- c(dont_remove, start)
  
  t_df[start,] = agg_var
    
  start <- start + 24
  end <- end + 24 
}

# Keep only the aggregated rows: 
site_df <- t_df 
site_df$rownumber = 1:nrow(t_df)
site_df <- site_df[site_df$rownumber %in% dont_remove,]

# Convert units from K to C
site_df$mx2t <- site_df$mx2t - 273.15

# Extract "mn2t" data for the site
site_series <- extract(met_nc_mn2t.b, SpatialPoints(cbind(t_lon, t_lat)), method='simple')
t_df <- data.frame(mn2t=t(site_series))

# Aggregate "mn2t" data (hourly -> daily timestep)
start <- 1
end <- 24
days <- length(site_series) # 8784
dont_remove <- c()

while (end <= days) {
  print(end)
  agg_var <- sum(site_series[,start:end]) / 24
  dont_remove <- c(dont_remove, start)
  
  t_df[start,] = agg_var
  
  start <- start + 24
  end <- end + 24 
}

# Keep only the aggregated rows: 
t_df$rownumber = 1:nrow(t_df)
t_df <- t_df[t_df$rownumber %in% dont_remove,]
site_df$mn2t <- t_df$mn2t

# Convert units from K to C
site_df$mn2t <- site_df$mn2t - 273.15

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
  
  

  