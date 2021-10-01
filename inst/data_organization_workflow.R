require(tidyverse)
require(RSPlantingDate)

# Organize the Beck's data and the extracted NDVI values

rm(list=ls())
set.seed(102396)
setwd('/Volumes/UIUC/libs/RS_PlantingDate')

# 1: Compile data from CSV files 
csvs = list.files('inst/data/NDVI/', full.names = T)
full = read.csv(csvs[1])
for (c in csvs[-1]){
  temp = read.csv(c)
  full = rbind(full, temp)
}
full = full %>% 
  mutate(Date = as.Date(substr(time, 1,10))) %>% 
  dplyr::select(ID, NDVI, Date)
rm(csvs,c,temp)

# 2: Add Beck's site information
IDs = unique(full$ID)
siteInfo = read.csv('inst/data/Becks_site.csv') %>% 
  dplyr::filter(TestPlotId %in% IDs) %>% 
  dplyr::mutate(ID = TestPlotId) %>%
  dplyr::select(ID, Latitude, Longitude, Family) %>%
  left_join(read.csv('inst/data/PSIMS_Mng.csv') %>% 
              mutate(HARVESTED = as.Date(HARVESTED),
                     PLANTED = as.Date(PLANTED)) %>%
              dplyr::select(Id,PLANTED, HARVESTED) %>%
              dplyr::filter(Id %in% IDs), 
            by = c('ID' = 'Id'))   
  
# 3: Add missing date values for planting date => for these values, Alex and I manually had to look at each PDF and extract the needed information
mine = read.csv('inst/data/missingDates.csv') %>% 
  mutate(PLANTED = as.Date(PLANTED, format = '%m/%d/%y'),
         HARVESTED = as.Date(HARVESTED,format = '%m/%d/%y'))
missing = c()
for (i in 1:nrow(siteInfo)){
  if (is.na(siteInfo$PLANTED[i]) | is.na(siteInfo$HARVESTED[i])){
    X = which(mine$ID == siteInfo$ID[i])
    
    if (length(X) == 0){
      missing = c(missing, i)
      next
    }
    
    siteInfo$PLANTED[i] = mine$PLANTED[X]
    siteInfo$HARVESTED[i] = mine$HARVESTED[X]
  }
  
  string = toString(siteInfo$PLANTED[i])
  string2 = toString(siteInfo$HARVESTED[i])
  siteInfo$PLANTED[i] = as.Date(paste0('2017-',substr(string, nchar(string)-4, nchar(string))))
  siteInfo$HARVESTED[i] = as.Date(paste0('2017-',substr(string2, nchar(string2)-4, nchar(string2))))
}

if (length(missing) > 0){
  print('There are some site IDs that do not have planting dates.')
  break
}   

# 4: Update the latitude and longitude coordinates for this site 
load('inst/data/updated_coordinates.Rdata')
temp = left_join(siteInfo %>% dplyr::select(-Latitude, -Longitude),
                 updated4 %>% dplyr::select(ID, Latitude, Longitude), by = c('ID'))

# 5: Combine information into one data frame
temp = left_join(full, temp, by = 'ID') %>% 
  filter(!is.na(PLANTED)) %>% 
  mutate(das = as.integer(Date - PLANTED), 
         dos = as.integer(PLANTED - as.Date('2016-12-31')), 
         doh = as.integer(HARVESTED - as.Date('2016-12-31'))) 

write.csv(temp, file = 'inst/data/organized/cleanedData.csv')
