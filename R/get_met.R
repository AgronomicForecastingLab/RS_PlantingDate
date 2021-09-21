require(tidyverse)

#' Given a plot location (lat/lon) and year, calculates monthly precipitation
#' and cumulative GDD (Growing Degree Days).
#' 
#' Retrieves precipitation, temperature data from the meterological database
#' 'ERA5 hourly data on single levels from 1979 to present'. 
#'
#' @param lat plot latitude
#' @param lon plot longitude
#' @param year year that we need data for
#' @return A data frame containing cumulative GDD and precipitation for each month of the year
#' @export
get_met <- function(lat, lon, year) {
  load('inst/data/test_train_data.Rdata')
  all = rbind(test, train)
  tile <- c()
  
  # Determine tile location of the given lat/lon.
  # Load directory of available pSIMS tiles.
  pSIMS_extents <- read.csv2("inst/pSIMS_extents.csv", sep = ',')
  
  # Verify that given lat/lon falls within a pSIMS tile.
  for (i in 1:nrow(all)) {
    print(i)
    lat = all$Latitude[i]
    lon = all$Longitude[i]
    
    current_tile <- pSIMS_extents %>%
      dplyr::filter(xmax >= lon,
                    xmin <= lon,
                    ymin <= lat,
                    ymax >= lat)
    if (nrow(current_tile) < 1)
      stop("The lat/lon does not fall into any pSIMS tile")
    tile = c(tile, current_tile$name)
  }
  unique(tile)
  
  # Extract tile lat and long pSIMS indices.
  idxs <- (current_tile %>%
             pull(name) %>%
             strsplit('/'))[[1]] %>% as.numeric()
  
  tlatidx <- idxs[1]
  tlonidx <- idxs[2]
  
  # Tile deltas ____ (?)
  tlatdelta <- 120
  tlondelta <- 120
  
  # what tile we want to simulate
  
  ## Store lat/lon degrees associated with given tile.
  lats <- seq.int(current_tile$ymin, current_tile$ymax, by=1)
  lons <- seq.int(current_tile$xmin, current_tile$xmax, by=1)
 
}