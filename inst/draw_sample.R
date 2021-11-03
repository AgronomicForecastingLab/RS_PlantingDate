require(tidyverse)
require(maps)
require(mapdata)

# Draw a random sample of 300 sites for planting date analysis 

rm(list=ls())
set.seed(102396)
setwd('/Volumes/UIUC/libs/RS_PlantingDate')

# 1: Get full site information
siteInfo = read.csv('inst/data/Becks_site.csv') %>% 
  dplyr::mutate(ID = TestPlotId) %>%
  dplyr::select(ID, Latitude, Longitude, Family) %>%
  left_join(read.csv('inst/data/PSIMS_Mng.csv') %>% 
              mutate(HARVESTED = as.Date(HARVESTED),
                     PLANTED = as.Date(PLANTED)) %>%
              dplyr::select(Id,PLANTED, HARVESTED), 
            by = c('ID' = 'Id')) %>% 
  mutate(Year = lubridate::year(PLANTED)) %>%
  filter(Family == 1,
         PLANTED < HARVESTED, 
         Year %in% c(2017:2019))

# There are some plots in 2017 that we know are bad. 
bads = read.table('/Volumes/UIUC/libs/RS_PlantingDate/inst/data/unclearPlots.txt') 
siteInfo = siteInfo %>% filter(!(ID %in% bads))

# We have already fixed locations for 2017 so let's add those adjusted coordinates in prior to plotting
load('inst/data/updated_coordinates.Rdata')
updated = updated4 %>% distinct(ID, .keep_all = TRUE) %>% 
  mutate(newLat = Latitude, newLon = Longitude) %>%
  dplyr::select(ID, newLat, newLon)
rm(updated4)
temp = left_join(siteInfo %>% filter(Year == 2017),
                 updated, 
                 by = 'ID') 
# all IDs are accounted for here so that's a good sign! 

# let's change the coordinates in the original siteInfo data frame for each of these IDs
for (i in 1:nrow(temp)){
  this = temp[i,]
  ind = which(siteInfo$ID == temp$ID[i])
  siteInfo$Latitude[i] = temp$newLat[i]
  siteInfo$Longitude[i] = temp$newLon[i]
}

# Plot points on a map
sts = c('minnesota','wisconsin','north dakota','south dakota',
        'nebraska','kansas','ohio','michigan','indiana','illinois',
        'iowa', 'missouri', 'kentucky', 'tennessee', 'west virginia',
        'pennsylvania','arkansas', 'oklahoma')
usa = map_data('state') %>% filter(region %in% sts)
ggplot(data = usa) + 
  geom_polygon(aes(x = long, y = lat, group = group), alpha = 0.2, color = 'black',
               fill = 'white') + 
  geom_point(data = siteInfo, aes(x = Longitude, y = Latitude), color = 'blue', alpha = 0.5) +
  labs(x = NULL, y = NULL, title = 'Full Data Coverage (2017-2019)') +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84'))


# Random stratified sample by year and location 
summary(siteInfo$Latitude)
summary(siteInfo$Longitude)
latGrid = seq(35.25,46.41,length.out = 3)
lonGrid = seq(-99.05,-79.68,length.out=3)
grid = expand.grid(latGrid, lonGrid)
ggplot(data = usa) + 
  geom_polygon(aes(x = long, y = lat, group = group), color = 'black',
               fill = 'white') + 
  geom_point(data = siteInfo, aes(x = Longitude, y = Latitude, color = as.factor(Year)), alpha = 0.5) +
  labs(x = NULL, y = NULL, title = 'Stratified Sample (2017-2019)') + 
  geom_vline(data = grid, aes(xintercept = Var2), color = 'red') + 
  geom_hline(data = grid, aes(yintercept = Var1), color = 'red') +
  coord_cartesian(xlim = c(-99.15, -79.78), ylim = c(37.15, 46.51))

# Now classify based on grid cell
siteInfo$class = NA
for (i in 1:nrow(siteInfo)){
  lat = siteInfo$Latitude[i]
  lon = siteInfo$Longitude[i]
  if (lat > 40.83){
    if (lon < -89.365){
      siteInfo$class[i] = 1
    }else{
      siteInfo$class[i] = 2
    }
  }else{
    if (lon < -89.365){
      siteInfo$class[i] = 4
    }else{
      siteInfo$class[i] = 3
    }
  }
}

props = siteInfo %>% group_by(Year) %>%
  summarize(total = n()) %>% 
  left_join(siteInfo %>% group_by(Year, class) %>%
              summarize(totalclass = n()) %>%
              ungroup(),
            by = c('Year')) %>%
  mutate(prop = totalclass/total)
  #ungroup() %>% 
  #ggplot(aes(x = Year, y = totalclass, fill = as.factor(class))) + 
  #geom_col()


# Extracting stratified random sample won't work evenly for some because of limited
# data... so let's extract from each grid cell based on the relative proportion of data
# in each per year. 
samp = c()
for (y in 2017:2019){
  this.year = siteInfo %>% filter(Year == y)
  for (cl in 1:4){
    this.class = (this.year %>% filter(class == cl))$ID
    this.n = round((props %>% filter(Year == y, class == cl))$prop * 100)
    samp = c(samp,sample(this.class, size = this.n, replace = FALSE))
  }
}
sampleInfo = siteInfo %>% filter(ID %in% samp)
ggplot(data = usa) + 
  geom_polygon(aes(x = long, y = lat, group = group), alpha = 0.2, color = 'black',
               fill = 'white') + 
  geom_point(data = sampleInfo, aes(x = Longitude, y = Latitude), color = 'purple', alpha = 0.6) +
  labs(x = NULL, y = NULL, title = 'Sample Coverage (2017-2019)', color = NULL) +
  theme_bw() +
  theme(text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) 


write.csv(sampleInfo, file = 'inst/data/sampleInfo.csv')

