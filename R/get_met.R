#' Takes a given plot location and year and calculates monthly precipitation and cumulative GDD 

#'
#' @param lat plot latitude
#' @param lon plot longitude
#' @param year year that we need data for
#' @return A data frame containing cumulative GDD and precipitation for each month of the year
#' @export
# get_met <- function(lat, lon, year) {

require(tidyverse)

load('inst/data/test_train_data.Rdata')
all = rbind(test, train)
tile <- c()

#Find out in which tile does this lat/lon fall
pSIMS_ext_path <- system.file("Utils", 'pSIMS_extents.csv', package = "pSIMSSiteMaker")
PSIMS_extents <- read.csv2(pSIMS_ext_path, sep=',')

for (i in 1:nrow(all)){
  print(i)
  lat = all$Latitude[i]
  lon = all$Longitude[i]
  
  current_tile <- PSIMS_extents %>%
    dplyr::filter(
      xmax>= lon,
      xmin<=lon,
      ymin <= lat,
      ymax>=lat
    )
  if(nrow(current_tile)<1) stop("The lat/lon does not fall into any pSIMS tile")
  tile = c(tile, current_tile$name)
}
unique(tile)

idxs <- (current_tile%>%
           pull(name) %>%
           strsplit('/'))[[1]] %>% as.numeric()

tlatidx <- idxs[1]
tlonidx <- idxs[2]
# Tile deltas and indices
tlatdelta <- 120
tlondelta <- 120
# what tile we want to simulate

## USER EDIT: what are the associated lat and long degree boundaries of these tiles?
lats  <-  c(current_tile$ymin, current_tile$ymax)
lons  <-  c(current_tile$xmin, current_tile$xmax)
split       <-  1
slatidx     <-  1
slonidx     <-  1
tslatdelta  <-  tlatdelta / split
tslondelta  <-  tlondelta / split
tslatidx    <-  split * (tlatidx - 1) + slatidx
tslonidx    <-  split * (tlonidx - 1) + slonidx
latlines  <-  rev(seq(lats[1], lats[2], length.out = tslatdelta / latdelta + 1))
lonlines  <- seq(lons[1], lons[2], length.out = tslondelta / londelta + 1)
# what coordinates do you need?

# find i and j
latdiffs  <- latlines - latneed
min  <- min(latdiffs[which(latdiffs > 0)])
i  <- which(latdiffs == min)
londiffs  <- lonlines - lonneed
min  <- max(londiffs[which(londiffs < 0)])
j  <- which(londiffs == min)
# Subtile deltas and indices
latidx  <- ceiling((tlatdelta * (tlatidx - 1) + tslatdelta * (slatidx - 1) + latdelta * i) / latdelta)
lonidx  <- ceiling((tlondelta * (tlonidx - 1) + tslondelta * (slonidx - 1) + londelta * j) / londelta)