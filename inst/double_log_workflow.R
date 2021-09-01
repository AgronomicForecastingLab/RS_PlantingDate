# Double-logistic equation prediction workflow:

devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies =
                           TRUE)

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
for (i in 1:length(IDs)) {
  this <- temp %>% filter(ID == IDs[i])
  if (nrow(this) >= 5)
    finalIDs = c(finalIDs, IDs[i])
}
temp <- temp %>% filter(ID %in% finalIDs)

# Remove any non-corn sites (corn is 'Family 1').
temp = temp %>% filter(Family == 1)
cornIDs = unique(temp$ID)

# Randomly select training and test data.
set.seed(102396)
trainIDs = sample(cornIDs, length(cornIDs) / 2, replace = FALSE)
train = temp %>% filter(ID %in% trainIDs)
test = temp %>% filter(!(ID %in% trainIDs))

train = train %>% mutate(
  das2 = das ^ 2,
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date('2016-12-31'))
)

rm(i, overwrite, cornIDs, this)

# Fit the double-logistic function to each site.
dlog.data = data.frame(
  ID = trainIDs,
  mn = rep(NA, length(trainIDs)),
  # winter NDVI (minimum)
  mx = rep(NA, length(trainIDs)),
  # maximum NDVI
  sos = rep(NA, length(trainIDs)),
  # start of season
  rsp = rep(NA, length(trainIDs)),
  # initial slope
  eos = rep(NA, length(trainIDs)),
  # end of season
  rau = rep(NA, length(trainIDs)),
  # ending slope
  dos = rep(NA, length(trainIDs)),
  # day of sowing
  lat = rep(NA, length(trainIDs)),
  # latitude
  rmse = rep(NA, length(trainIDs)),
  # RMSE value for fitted model
  maxDOY = rep(NA, length(trainIDs)),
  # day of year where NDVI is at maximum following model
  emerg = rep(NA, length(trainIDs))
) # "emergence" = first day where function derivative > 0.0001

# Iterate through all training plots, fits a double logistic, and records the
# fitted parameters.
dlog.data <- fit_double_logistic(dlog.data)

# Consider the relationship between variables and DOS to check for potential use
# in predicting planting date.
# Currently, we remove any observations for which the RMSE value is too high.
ggpairs(dlog.data %>% filter(rmse < 0.2), columns = c(2:7, 11, 12, 9, 8))
mod1 = lm(dos ~ eos + lat , data = dlog.data %>% filter(rmse < 0.15))
summary(mod1)
# RMSE: 13.5

# Optimize the dl function for each test plot with corn.
test.corn = test %>%
  mutate(dos = as.integer(as.Date(PLANTED) - as.Date('2016-12-31')),
         DOY = as.integer(as.Date(Date) - as.Date('2016-12-31')))
test.cornids = unique(test.corn$ID)
try = data.frame(ID = NULL, obs = NULL, pred = NULL)
for (i in 1:length(test.cornids)) {
  print(i)
  this = test.corn %>% filter(ID == test.cornids[i]) %>% distinct(DOY, .keep_all = TRUE)
  
  # fit model based on Beck's equation
  opt = try(FitDoubleLog(
    x = this$NDVI,
    t = this$DOY,
    weighting = TRUE,
    plot = FALSE
  ),
  silent = T)
  
  # if we are not able to fit the model, we should try once more and then skip it if still bad
  if (is.na(opt[2])) {
    opt = try(FitDoubleLog(
      x = this$NDVI,
      t = this$DOY,
      weighting = TRUE,
      plot = FALSE
    ),
    silent = T)
    if (is.na(opt[2])) {
      try = rbind(try, data.frame(
        ID = this$ID[1],
        obs = this$dos[1],
        pred = NA
      ))
      next
    }
  }
  
  this.pred = predict.lm(mod1,
                         data.frame(
                           eos = opt$params[5],
                           lat = this$Latitude[1],
                           dos = NA
                         ))
  try = rbind(try,
              data.frame(
                ID = this$ID[1],
                obs = this$dos[1],
                pred = this.pred
              ))
}
