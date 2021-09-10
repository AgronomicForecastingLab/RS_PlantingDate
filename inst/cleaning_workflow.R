require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)
require(RSPlantingDate)

# Data cleaning workflow 

rm(list=ls())

# The Beck's Hybrids dataset should be named 'cleanedData.csv'.
# We need to remove duplicates of data points for the correct function of our workflow 
temp <- read.csv('inst/data/cleanedData.csv') %>% 
  mutate(Date = as.Date(Date)) %>%
  distinct(ID, Date,.keep_all = TRUE)

# Remove any sites with late DOS values or with DOS values after harvest (dataset errors).
temp <- temp %>% filter(dos < 182, dos < doh)

# Let's check to see how the cleaning function works
#ids = sample(unique(temp$ID), 12, replace=F)
#before = temp %>% filter(ID %in% ids)
#before$check = 'before'

# Remove any data points that are not monotonically increasing/decreasing 
IDs = unique(temp$ID)

if (!file.exists('inst/data/monotonic_corn_data.Rdata')){
  temp <- clean_monotonic(orig_df = temp)
  
  # Remove any sites with fewer than 5 data points.
  finalIDs = c()
  for (i in 1:length(IDs)) {
    this <- temp %>% filter(ID == IDs[i])
    if (nrow(this) >= 5)
      finalIDs = c(finalIDs, IDs[i])
  }
  temp <- temp %>% filter(ID %in% finalIDs)
  
  # Remove any non-corn sites (corn is 'Family 1').
  temp = temp %>% filter(Family == 1)
  
  save(temp, file = 'inst/data/monotonic_corn_data.Rdata')
}else{
  load('inst/data/monotonic_corn_data.Rdata')
}

# Plot!
#after = temp %>% filter(ID %in% ids)
#after$check = 'after'
#check = rbind(before,after) %>% mutate(Date = as.Date(Date))
#ggplot(check)+
#  geom_line(aes(x = Date, y = NDVI, color = check), alpha = 0.5) + 
#  geom_point(aes(x = Date, y = NDVI, color = check)) + 
#  facet_wrap(~ID, nrow = 4)

cornIDs = unique(temp$ID)

# Randomly select training and test data.
set.seed(102396)
trainIDs = sample(cornIDs, length(cornIDs) / 2, replace = FALSE)
train = temp %>% filter(ID %in% trainIDs)
test = temp %>% filter(!(ID %in% trainIDs))

train = train %>% mutate(
  #  das2 = das ^ 2,
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date('2016-12-31'))
)

test = test %>% mutate(
  #  das2 = das ^ 2,
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date('2016-12-31'))
)

save(test, train, file = 'inst/data/test_train_data.Rdata')
