require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)
require(RSPlantingDate)

# Double-logistic equation prediction workflow:

#devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies = TRUE)
#overwrite = F
rm(list=ls())
################################################
#### Preparing data
################################################

# The Beck's Hybrids dataset should be named 'cleanedData.csv'.
# We need to remove duplicates of data points for the correct function of our workflow 
temp <- read.csv('inst/cleanedData.csv') %>% 
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
temp <- clean_monotonic(orig_df = temp)

# Plot!
#after = temp %>% filter(ID %in% ids)
#after$check = 'after'
#check = rbind(before,after) %>% mutate(Date = as.Date(Date))
#ggplot(check)+
#  geom_line(aes(x = Date, y = NDVI, color = check), alpha = 0.5) + 
#  geom_point(aes(x = Date, y = NDVI, color = check)) + 
#  facet_wrap(~ID, nrow = 4)

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
#  das2 = das ^ 2,
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date('2016-12-31'))
)

rm(i, overwrite, cornIDs, this, finalIDs, temp, this)

################################################
#### Fitting training dataset to double logistic
################################################

# Fit the double-logistic function to each site.
dlog.data = data.frame(
  ID = NULL,
  mn = NULL,
  # winter NDVI (minimum)
  mx = NULL,
  # maximum NDVI
  sos = NULL,
  # start of season
  rsp = NULL,
  # initial slope
  eos = NULL,
  # end of season
  rau = NULL,
  # ending slope
  dos = NULL,
  # day of sowing
  lat = NULL,
  # latitude
  maxDOY = NULL
  # day of year where NDVI is at maximum following model
) 

# Iterate through all training plots, fit a double logistic, and record the fitted parameters in dlog.data.
pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(trainIDs)) {
  setTxtProgressBar(pb, i/length(trainIDs))
  this = train %>% filter(ID == trainIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  res = fit_double_logistic(x = this$NDVI, t = this$DOY)
  if (any(!is.na(res))) dlog.data = rbind(dlog.data, res)
}

# There are still some outliers here which might be associated with poor fitting data => let's remove them 
ggplot(dlog.data) + geom_point(aes(x = eos, y = dos))
dlog.data = dlog.data %>% filter(eos > 150, eos < 335)

# Consider the relationship between variables and DOS to check for potential use in predicting planting date
ggpairs(dlog.data %>% dplyr::select(mn,mx,sos,rsp,eos,rau,lat,maxDOY,dos))
ggpairs(dlog.data %>% dplyr::select(dos, sos, eos, lat, maxDOY))

int <- lm(dos ~ 1, data=dlog.data)
all <- lm(dos ~ eos*lat*maxDOY, data=dlog.data)

# Perform forward stepwise regression
dlogMod <- step(int, direction='forward', scope=formula(all), trace=0)

################################################
#### Fitting training dataset to Gaussian function
################################################

# fit bell curve to each site 
norm.data = data.frame(ID = NULL,
                       mu = NULL,
                       sigma = NULL,
                       k = NULL,
                       dos = NULL,
                       lat = NULL)

# Iterate through all training plots, fit a normal, and record the fitted parameters in norm.data
pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(trainIDs)) {
  setTxtProgressBar(pb, i/length(trainIDs))
  this = train %>% filter(ID == trainIDs[i]) %>% 
    distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  res = fit_normal(x = this$NDVI, t = this$DOY)
  if (any(!is.na(res))) norm.data = rbind(norm.data, res)
}

# There are still some outliers here which might be associated with poor fitting data => let's remove them 
ggplot(norm.data) + geom_point(aes(x = mu, y = dos))
ggplot(norm.data) + geom_point(aes(x = sigma, y = dos))
norm.data = norm.data %>% filter(mu > 100, mu < 300, sigma < 200)

# Consider the relationship between variables and DOS to check for potential use in predicting planting date
ggpairs(norm.data %>% dplyr::select(mu,sigma,k,dos))

int <- lm(dos ~ 1, data=norm.data)
all <- lm(dos ~ lat*mu*sigma*k, data=norm.data)

# Perform forward stepwise regression
normMod <- step(int, direction='forward', scope=formula(all), trace=0)

################################################
#### Testing the two new functions at each data point 
################################################

# Organize test data
test = test %>% mutate(
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date('2016-12-31'))
)
testIDs = unique(test$ID)

# Predictions
pred.data = data.frame(ID = NULL,
                       obs = NULL,
                       dlog = NULL,
                       norm = NULL)

pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(testIDs)){
  setTxtProgressBar(pb, i/length(testIDs))
  this = test %>% filter(ID == testIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  
  # first dlog
  res = fit_double_logistic(x = this$NDVI, t= this$DOY)
  if (any(!is.na(res))){
    pred1 = predict.lm(dlogMod, res)
  }
  
  # then norm
  res = fit_normal(x = this$NDVI, t = this$DOY)
  if (any(!is.na(res))){
    pred2 = predict.lm(normMod, res)
  }

  pred.data = rbind(pred.data,
                    data.frame(ID = testIDs[i],
                               obs = this$dos[1],
                               dlog = pred1,
                               norm = pred2))
  
}

pred.data = pred.data %>% filter(norm > 90, norm < 182)

pred.data %>% 
  mutate(sqerror1 = (norm-obs)^2,
         sqerror2 = (dlog-obs)^2) %>%
  dplyr::summarize(norm_RMSE = sqrt(mean(sqerror1)),
                   dlog_RMSE = sqrt(mean(sqerror2)))
ggplot(pred.data) + 
  geom_point(aes(x = obs, y = dlog), color = 'blue')+
  geom_point(aes(x = obs, y = norm), color = 'red')+
  geom_abline(slope = 1, intercept = 0)

