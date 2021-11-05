require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)
require(RSPlantingDate)
require(devtools)
require(ncdf4)
require(ncmeta)
require(stars)
require(units)
require(ecmwfr)
require(ggplot2)
require(RNetCDF)
require(raster)

install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies=TRUE)

# Data cleaning and fitting workflow with the double logistic function
rm(list=ls())
set.seed(102396)

##############################
# I. Prepare the data 
##############################

# The Beck's Hybrids dataset should be named 'cleanedData.csv'.
# We need to remove duplicates of data points for the correct function of our workflow. I'm not sure where these doubles came from 
temp <- read.csv('inst/data/organized/cleanedData.csv') %>% 
  mutate(Date = as.Date(Date)) %>%
  distinct(ID, Date,.keep_all = TRUE)

# nrow(distinct(temp, ID)) # 280 sites in Beck's df

# Remove any sites with late DOS values or with DOS values after harvest (dataset errors).
temp <- temp %>% filter(dos < 182, dos < doh)

# Remove any non-corn sites (corn is 'Family 1').
temp = temp %>% filter(Family == 1)
  distinct(ID, Date,.keep_all = TRUE) %>%
  dplyr::select(-X.1, -X)

# Already removed any sites with late DOS values or with DOS values after harvest (dataset errors).
# Already removed any non-corn sites (corn is 'Family 1')

# Store a version of `temp` before data cleaning to use in later comparison. 
before = temp
IDs = unique(temp$ID)
before = before %>% mutate(
  Date = as.Date(Date),
  DOY = as.integer(Date - as.Date(paste0(Year-1,'-12-31'))),
  DAS =  as.integer(Date - as.Date(PLANTED))
)

# 1: First let's remove outliers (MEAN +/- 3*SD)
quants = data.frame(DAS = NULL, mn = NULL, upper = NULL, lower = NULL)
incs = seq(-170, 270, 10)
# determine the bounds 
for (i in 2:length(incs)){
  vals = before %>% filter(DAS < incs[i], DAS >= incs[i-1]) 
  std = sd(vals$NDVI)
  mn = mean(vals$NDVI)
  quants = rbind(quants, 
                 data.frame(mn = mn, 
                            DAS = mean(incs[(i-1):i]),
                            upper = mn + (3*std),
                            lower = mn - (3*std)))
}

# Figure: Demonstrate the bounds of the data where we need to perform removals 
ggplot(before) + 
  geom_point(aes( x= DAS, y = NDVI), alpha = 0.15) +
  geom_ribbon(data = quants, aes( x = DAS, ymin = lower, ymax = upper, y = mn),color = 'red', linetype = 'dashed', fill = 'red',alpha = 0.05) +
  geom_line(data = quants, aes( x= DAS, y = mn), color = 'red', size = 2)

# now remove points outside of the bounds 
mid = before[1,]
for (i in 2:length(incs)){
  vals = before %>% filter(DAS < incs[i], 
                           DAS >= incs[i-1], 
                           NDVI <= quants$upper[i],
                           NDVI >= quants$lower[i])
  mid = rbind(mid, 
              vals)
}
mid = mid[-1,] %>% select(-X.1)

# Figure: Demonstrate the cleaned data with the bounds 
ggplot(mid) + 
  geom_point(aes( x= DAS, y = NDVI), alpha = 0.15) +
  geom_ribbon(data = quants, aes( x = DAS, ymin = lower, ymax = upper, y = mn),color = 'red', linetype = 'dashed', fill = 'red',alpha = 0.05) +
  geom_line(data = quants, aes( x= DAS, y = mn), color = 'red', size = 2)

# 2: Remove any values before April with NDVI values greater than ???
#bad <- mid %>% filter(Date < as.Date('2017-05-01'), NDVI > 0.375)
#mid <- mid %>% filter(!(X %in% bad$X))

# 3: Apply spline cleaning function here 
#temp <- clean_spline(orig_df = mid %>% dplyr::select(-Year,-class,-DOY,-DAS))

temp = mid
IDs = unique(temp$ID)

# 4: Remove any sites with fewer than 7 data points.
finalIDs = c()
for (i in 1:length(IDs)) {
  this <- temp %>% dplyr::filter(ID == IDs[i])
  if (nrow(this) >= 7)
    finalIDs = c(finalIDs, IDs[i])
}
temp <- temp %>% dplyr::filter(ID %in% finalIDs)
save(temp, file = 'inst/data/organized/cleaned_corn_data.Rdata')

# Figure: Generate 16 plots comparing uncleaned and cleaned data.
cornIDs = unique(temp$ID)
#samples = sample(cornIDs, 16)
#compare_cleaned_plots(samples, before %>% dplyr::select(-X.1, -DOY, -class, -Year, -DAS), temp)
rm(quants, this, vals, finalIDs, IDs, incs, mn, std, i, mid)

# 5: Now, we need to extract the needed met information for each location. 
## THIS PART STILL NEEDS TO BE DONE. 

# 6: We split the data into 2 halves: test and train. 
# We can use the train dataset to explore and fit the data. 
# We can use the test dataset to see how the final model works. 
# Randomly select training and test data.
# Setting a seed in the script is important here so we can always get the same datasets 
# if we have to rerun this script. 

trainIDs = sample(cornIDs, length(cornIDs) / 2, replace = FALSE)
temp = temp %>% mutate(
  Date = as.Date(Date)
  #DOY = as.integer(Date - as.Date('2016-12-31'))
)
train = temp %>% filter(ID %in% trainIDs)
test = temp %>% filter(!(ID %in% trainIDs))

save(test, train, file = 'inst/data/organized/test_train_data.Rdata')
rm(cornIDs)

##############################
# II. Estimate emergence for each plot 
##############################

# 1: Fit the double-logistic function to each site
# Store the DLF parameters in the `dlog.data` df: 
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
  lon = NULL,
  # longitude
  maxDOY = NULL
  # day of year where NDVI is at maximum following model
) 

# Iterate through all training plots, fit a double logistic, 
# and record the fitted parameters in dlog.data.
pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(trainIDs)) {
  setTxtProgressBar(pb, i/length(trainIDs))
  print(trainIDs[i])
  this = train %>% filter(ID == trainIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  res = fit_double_logistic(x = this$NDVI, t = this$DOY)
  if (any(!is.na(res))) dlog.data = rbind(dlog.data, res)
}

max(dlog.data$sos) # 215.95
min(dlog.data$sos) # 105.24
mean(dlog.data$sos) # 166.33
sd(dlog.data$sos) # 15.54

#ggpairs(dlog.data %>% dplyr::select(mn,mx,sos,eos,rau,rsp,lat,lon,dos))

# 2: Limit our fits to the ones that are actually reasonable by setting limits on some of the parameters 
bounds = list(
  sos = list(lower = 100,
             upper = 200),
  eos = list(lower = 244,
             upper = 344),
  mn = list(lower = 0,
             upper = 0.375),
  mx = list(lower = 0.5,
             upper = 0.9)
)

dlog.data = dlog.data %>% filter(sos > bounds$sos$lower, sos < bounds$sos$upper,
                                 eos > bounds$eos$lower, eos < bounds$eos$upper, 
                                 mn > bounds$mn$lower, mn < bounds$mn$upper,
                                 mx > bounds$mx$lower, mx < bounds$mx$upper)

# Figure: This figure demonstrates a sample of fitted double logistic functions. 
ids = sample(unique(dlog.data$ID), 16)
fitted = data.frame(ID = NULL, DOY = NULL, pred = NULL)
for (i in 1:length(ids)){
  pars = as.vector(dlog.data %>% filter(ID == ids[i]) %>% dplyr::select(mn, mx, sos, rsp, eos, rau))
  fitted = rbind(fitted, 
                 data.frame(ID = rep(ids[i], 365),
                            DOY = 1:365,
                            preds = predictDLOG(pars, 1:365)))
}
check = temp %>% mutate(Date = as.Date(Date)) %>% filter(ID %in% ids)

# Dashed line shows fitted double logistic function for each plot 
ggplot(check) +
  geom_point(aes(x = DOY, y = NDVI)) + 
  geom_line(data = fitted, aes(x = DOY, y = preds), linetype = 'dashed')+
  scale_linetype_manual(values = c('cleaned'='dashed', 'fitted'='solid')) +
  facet_wrap(~ID, nrow = 4) +
  labs(color = 'Data Cleaning',
       linetype = NULL)

ggplot(site_for_test) +
  geom_point(aes(x=DOY, y=NDVI)) + 
  geom_line(aes(x = DOY, y = NDVI), linetype = 'solid', color = 'red')

# 3: Estimate emergence for each plot based on fitted double logistic curve 
## WE NEED TO FIGURE THIS PART OUT 
# I think we should have one function that estimates emergence based on a fitted double logistic function. 
# The script would just need the fitted parameters and would return a DOY for emergence. 

########################
# Predicting Emergence # ------------------------------------------------------
########################

# We can either use Sentinel data or the fitted double-logistic for this step?
# Test on `site_for_test` (ID 57874)
site_VI <- site_for_test$NDVI # time series VI

# We will use MACD(5, 10, 5). (Gao et al. 2020)
# I.e., S = 5 days ("fast" EMA), L = 10 days ("slow" EMA), 
# K = 5 days ("signal"/"average" EMA)

myema_short <- my
macd_ndvi <- MACD_VI(VI, 5, 10, 5)

##############################
# III. Use climate data to estimate planting date from emergence  
##############################

load('inst/data/finalSiteLocations.Rdata')

# Open the NetCDF so we can extract met data by variable:
met_nc_data <- nc_open('inst/data/met_data.nc')
met_nc_fname <- 'inst/data/met_data.nc'
# "mx2t": max. air temperature
met_nc_mx2t.b <- brick(met_nc_fname, varname="mx2t")
# "mn2t": min. air temperature
met_nc_mn2t.b <- brick(met_nc_fname, varname="mn2t")

# Extract site-specific time-series data from the raster brick 
met_df_71030 <- get_site_met(siteLoc, 71030, met_nc_mx2t.b, met_nc_mn2t.b)
# Remove the NA column from `site_df`
site_df <- site_df[,-c(1)]

# Remove the .nc data variable and close the met. ncfile
rm(met_nc_data)
nc_close('inst/data/met_data.nc')

# Calculate GDD (Growing Degree Days)
# https://ndawn.ndsu.nodak.edu/help-corn-growing-degree-days.html

# Daily Corn GDD (°C) = ((Daily Max Temp °C + Daily Min Temp °C)/2) - 10 °C
# Following constraints:
#   --> If daily Max Temp > 30 °C, it's set equal to 30 °C.
#   --> If daily Max or Min Temp < 10 °C, it's set equal to 10°C. 

gdd_list <- c()
for (i in 1:nrow(site_df)) {
  max_temp <- site_df[i,]$mx2t$mx2t
  min_temp <- site_df[i,]$mn2t$mx2t
  
  if (max_temp > 30) {
    max_temp = 30
  } else if (max_temp < 10) {
    max_temp = 10
  } 
  if (min_temp < 10) {
    min_temp = 10
  }
  
  gdd <- ((max_temp + min_temp)/2) - 10
  gdd_list <- c(gdd_list, gdd)
}
site_df$GDD <- gdd_list

################################################
# IV. Apply the final model to the test dataset to evaluate performance 
################################################

# Organize test data
testIDs = unique(test$ID)

# Predictions
pred.data = data.frame(ID = NULL,
                       obs = NULL,
                       emerg = NULL, 
                       pred = NULL)

pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(testIDs)){
  setTxtProgressBar(pb, i/length(testIDs))
  
  # get the data for this plot 
  this = test %>% filter(ID == testIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  
  # first, fit the double logistic function
  res = fit_double_logistic(x = this$NDVI, t= this$DOY)
  if (any(is.na(res))) next 
  
  # check if the fit is bad 
  res = res %>% filter(sos > bounds$sos$lower, sos < bounds$sos$upper,
                       eos > bounds$eos$lower, eos < bounds$eos$upper, 
                       mn > bounds$mn$lower, mn < bounds$mn$upper,
                       mx > bounds$mx$lower, mx < bounds$mx$upper)
  if (nrow(res) == 0) next
  
  # second, estimate emergence
  
  
  # third, estimate planting date 
  
  
  # save the predictions (need to add those predicted values to the save function)
  pred.data = rbind(pred.data,
                    data.frame(ID = testIDs[i],
                               obs = this$dos[1],
                               emerg = NA, 
                               pred = NA))
}

# Figure: observed vs. predicted planting date 
ggplot(pred.data) + geom_point(aes(x = obs, y = pred)) + 
  geom_abline(slope = 1, linetype = 'dashed') 

# Calculate RMSE of prediction 
pred.data %>% 
  mutate(sqerror1 = (norm-obs)^2,
         sqerror2 = (dlog-obs)^2) %>%
  dplyr::summarize(norm_RMSE = sqrt(mean(sqerror1)),
                   dlog_RMSE = sqrt(mean(sqerror2)))

################################################
# 4. This is just some old code in case we end up needing some of this stuff later on 
################################################

# # This is a code I used to perform forward stepwise regression to determine the best linear model 
# # We might not need to use this for this section, but I'm keeping it here just in case. 
# int <- lm(dos ~ 1, data=dlog.data)
# all <- lm(dos ~ sos*eos*maxDOY*lat*lon, data=dlog.data)
# 
# # Perform forward stepwise regression
# dlogMod <- step(int, direction='forward', scope=formula(all), trace=0)
# summary(dlogMod)
# 
# # Code to group planting dates into 5 categories for easier visualizations
# categorize2 = function(dos){
#   if (is.na(dos)) return(NA)
#   if(dos < 90) cat = 1
#   if(dos >= 90 & dos < 110) cat = 2
#   if(dos >= 110 & dos < 130) cat = 3
#   if(dos >= 130 & dos < 150) cat = 4
#   if(dos >= 150) cat = 5
#   return(cat)
# }
# plotDat$cat2 = sapply(plotDat$dos,categorize2)
# plotDat$cat2 = factor(plotDat$cat2, levels = c(1:5))
# 
# # Code to group residuals into 5 categories for easier visualizations
# categorize = function(res){
#   if (is.na(res)) return(NA)
#   if(res < -20) cat = '< -20'
#   if(res >= -20 & res < 0) cat = '-20-0'
#   if(res >= 0 & res < 20) cat = '0-20'
#   if(res >= 20 & res < 40) cat = '20-40'
#   if(res >= 40) cat = '> 40'
#   return(cat)
# }
# plotDat$cat = sapply(plotDat$value,categorize)
# plotDat$cat = factor(plotDat$cat, levels = c('< -20','-20-0','0-20','20-40','> 40'))




