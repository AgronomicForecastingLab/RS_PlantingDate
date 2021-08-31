# Double-logistic equation prediction workflow:

devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies=TRUE)

require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)

overwrite = F
# The Beck's Hybrids dataset should be named 'cleanedData.csv'.
temp <- read.csv('cleanedData.csv')

# Remove any sites with late DOS values. 
temp <- temp %>% filter(dos < 182)
# Remove any sites with DOS values after harvest (dataset error).
temp <- temp %>% filter(dos > doh) 

# Remove any data points that are not monotonically increasing/decreasing.
IDs = unique(temp$ID)
temp <- clean_monotonic(orig_df = temp, ID_list = IDs)

# Remove any sites with fewer than 5 data points. 
finalIDs = c()
for (i in 1:length(IDs)){
  this <- temp %>% filter(ID == IDs[i])
  if (nrow(this) >= 5) finalIDs = c(finalIDs, IDs[i])
}
temp <- temp %>% filter(ID %in% finalIDs)

# Remove any non-corn sites (corn is 'Family 1').
temp = temp %>% filter(Family == 1)
cornIDs = unique(temp$ID)

# Randomly select training and test data.
set.seed(102396)
trainIDs = sample(cornIDs, length(cornIDs)/2, replace = FALSE)
train = temp %>% filter(ID %in% trainIDs)
test = temp %>% filter(!(ID %in% trainIDs))

train = train %>% mutate(das2 = das^2,
                         Date = as.Date(Date),
                         DOY = as.integer(Date - as.Date('2016-12-31')))

rm(i, overwrite, cornIDs, this)



