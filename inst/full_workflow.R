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
rm(quants, this, vals, finalIDs, IDs, incs, mn, std, i, mid)

# 6: We split the data into 2 halves: test and train. 
#  - We can use the train dataset to explore and fit the data. 
#  - We can use the test dataset to see how the final model works.
#
# Note: We need to randomly select training and test data.
#  - Setting a seed in the script is important here so we can always get the 
#    same datasets if we have to rerun this script. 

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
  ID = NULL, # Beck's ID of site
  mn = NULL, # Winter NDVI (minimum)
  mx = NULL, # Maximum NDVI
  sos = NULL, # Start of season
  rsp = NULL, # Initial slope
  eos = NULL, # End of season
  rau = NULL, # Ending slope
  dos = NULL, # Day of sowing
  lat = NULL, # Latitude
  lon = NULL, # Longitude
  maxDOY = NULL # Day of year where NDVI is at maximum following model
) 

# Iterate through all training plots, fit a double logistic, 
# and record the fitted parameters in dlog.data.
pb = txtProgressBar(0, 1, style=3)
for (i in 1:length(trainIDs)) {
  setTxtProgressBar(pb, i/length(trainIDs))
  this = train %>% filter(ID == trainIDs[i]) %>% distinct(DOY, .keep_all = TRUE)
  if (nrow(this) == 0) next
  res = fit_double_logistic(x = this$NDVI, t = this$DOY)
  if (any(!is.na(res))) dlog.data = rbind(dlog.data, res)
}

dlog.data$doe <- NA

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
check <- temp %>% mutate(Date = as.Date(Date)) %>% filter(ID %in% ids)

# Dashed line shows fitted double logistic function for each plot 
ggplot(check) +
  geom_point(aes(x = DOY, y = NDVI)) + 
  geom_line(data = fitted, aes(x = DOY, y = preds), linetype = 'dashed')+
  scale_linetype_manual(values = c('cleaned'='dashed', 'fitted'='solid')) +
  facet_wrap(~ID, nrow = 4) +
  labs(color = 'Data Cleaning',linetype = NULL)

########################
# Predicting Emergence # ------------------------------------------------------
########################

# We can either use Sentinel/Harmonized Sentinel Landsat-8 data or the fitted 
# double-logistic for this step?
pb = txtProgressBar(0, 1, style=3)

for (i in 1:nrow(dlog.data)) {
  setTxtProgressBar(pb, i/nrow(dlog.data))
  
  site_id <- dlog.data[i,]$ID
  site = train %>% filter(ID == site_id)
  pars = as.vector(dlog.data %>% filter(ID == site_id)
                   %>% dplyr::select(mn, mx, sos, rsp, eos, rau))
  
  site_preds = data.frame(ID = rep(site_id, 365),
                            DOY = 1:365,
                            preds = predictDLOG(pars, 1:365))
  
  site_VI <- site_preds$preds # VI time-series for the site. 
  
  # We will use MACD(5, 10, 5). (Gao et al. 2020)
  # I.e., S = 5 days ("fast" EMA), L = 10 days ("slow" EMA), 
  # K = 5 days ("signal"/"average" EMA)
  macd_t <- MACD_VI(site_VI, 5, 10, 10)
  site_preds$macd_t <- macd_t[,1]
  site_preds$ema_macd_c <- macd_t[,2]
  site_preds$macd_div_t <- site_preds$macd_t - site_preds$ema_macd_c
  
  # Settings (0.0, 0.0) for (MACD_threshold, MACD_div_threshold) are standard 
  # settings to track trend changes in stock prices. 
  # We adjust to use settings (0.1, 0.01) here (the fitted DL has uniform tails).
  macd_threshold <- 0.1
  macd_div_threshold <- 0.0001
  
  # The settings of 0.01 and 0.85 for NDVI at the green-up dates represent a 
  # wide range of acceptable green-up conditions since NDVI and EVI2 at green-up
  # time are very low.
  min_VI_greenup <- 0.01
  max_VI_greenup <- 0.90
  VI_increase_doys <- c()
  
  max_ndvi_doy <- site_preds[which.max(site_preds$preds),]$DOY
  
  # We check for momentum points with the various conditions described by 
  # Gao et al. (2020).
  # (ex:
  #   -- Is the NDVI at this doy within a reasonable range of NDVI values?, 
  #   -- Does the MACD begin to increase exponentially?, 
  #   -- Is the MACD within a reasonable range? 
  # )
  for (j in 11:max_ndvi_doy) {
    in_VI_range <-
      site_preds$preds[j] > min_VI_greenup &&
      site_preds$preds[j] < max_VI_greenup
    macd_div_x_intercept <-
      (
        site_preds$macd_div_t[j - 1] < macd_div_threshold &&
          site_preds$macd_div_t[j] > macd_div_threshold
      ) 
    
    macd_within_range <- macd_t[j] < macd_threshold
    if (in_VI_range && macd_div_x_intercept && macd_within_range) {
     VI_increase_doys <- c(VI_increase_doys, j)
    }
  }
  
  if (length(VI_increase_doys) == 0) {
    print(i)
    next
  }
  
  momentums <- data.frame(DOY = VI_increase_doys, NDVI = NA)
  for (k in 1:length(VI_increase_doys)) {
    site_row <- site_preds[VI_increase_doys[k],]
    momentums[k,]$NDVI = site_row$preds
  }

  ## Use a horizontal line to mark `macd_div_t` peak
  emerg_row <- site_preds[VI_increase_doys[1],]
  
  ## Use a horizontal line to mark actual planting date
  site_dos <- (temp %>% filter(ID == trainIDs[i]))[1,]$dos
  
  # Store 'doe' for the site (a.k.a. day of emergence).
  dlog.data[dlog.data$ID == site_id,]$doe <- emerg_row$DOY 
}

# Take a random sample of 10 predicted site emergences to plot.
ids <- sample(dlog.data$ID, 10, replace = FALSE, prob = NULL)

fitted = data.frame(ID = NULL, DOY = NULL, pred = NULL)
for (i in 1:length(ids)){
  pars = as.vector(dlog.data %>% filter(ID == ids[i]) %>% dplyr::select(mn, mx, sos, rsp, eos, rau))
 
  fitted = rbind(fitted, 
                 data.frame(ID = rep(ids[i], 365),
                            DOY = 1:365,
                            preds = predictDLOG(pars, 1:365)))
}

sample_sites <- temp %>% mutate(Date = as.Date(Date)) %>% filter(ID %in% ids)
sample_sites$doe <- NA
for (i in 1:nrow(sample_sites)) {
  dlog_row <- dlog.data %>% filter(ID == sample_sites[i,]$ID)
  sample_sites[i,]$doe <- dlog_row$doe
}

p <- ggplot(sample_sites) +
  geom_point(aes(x = DOY, y = NDVI)) +
  geom_line(aes(x = doe, y = NDVI), linetype = 'solid', color='red') + 
  geom_line(aes(x = dos, y = NDVI), linetype = 'solid', color='green') +
  geom_line(data = fitted, aes(x = DOY, y = preds), linetype = 'dashed') +
  scale_linetype_manual(values = c('cleaned'='dashed', 'fitted'='solid')) +
  facet_wrap(~ID, nrow = 2) +
  labs(color = 'Data Cleaning',linetype = NULL) +
  ylim(0.1, 0.9) + 
  ggtitle("Site DOE Pred Sample")

p + theme(
  plot.title = element_text(color="red", size=14, face="bold.italic",
                            hjust = 0.5))

##############################
# III. Use climate data to estimate planting date from emergence  
##############################

# Open the NetCDF so we can extract met data by variable:
met_nc_data <- nc_open('inst/data/met_data.nc')
met_nc_fname <- 'inst/data/met_data.nc'

met_nc_mx2t.b <- brick(met_nc_fname, varname="mx2t") # "mx2t" (max. air temp.)
met_nc_mn2t.b <- brick(met_nc_fname, varname="mn2t") # "mn2t" (min. air temp.)

# We will calculate the AGDD (accumulated growing degree days) between the 
# observed planting date (dos) and the predicted emergence date (doe).
dlog.data$AGDD <- NA
pb = txtProgressBar(0, 1, style=3)

for (i in 1:nrow(dlog.data)) {
  setTxtProgressBar(pb, i/nrow(dlog.data))
  
  site_id <- dlog.data[i,]$ID
  site_rows <- train %>% filter(ID == site_id)
  site_row <- site_rows[1,]

  met_df <- get_site_met(site_row, met_nc_mx2t.b, met_nc_mn2t.b)
  met_df$doy <- NA
  met_df$year <- NA
  
  for (j in 1:nrow(met_df)) {
    date <- rownames(met_df[j,]$mx2t)
    year <- substr(date, 2, 5)
    month <- substr(date, 7, 8)
    day <- substr(date, 10, 11)
    str_date <- paste(year, month, day, sep="-")
    tmp <- as.POSIXlt(str_date, format = "%Y-%m-%d", tz = "GMT")
    julian_days <- tmp$yday
    met_df[j,]$doy <- julian_days
    met_df[j,]$year <- year
  }
  
  dos_year <- substr(site_row$Date, 1, 4)
  met_df_year <- met_df %>% filter(year == as.integer(dos_year))
  
  # Keep only the weather obsv. within the [DOS, DOE] interval:
  met_df_year <- met_df_year %>% filter(doy <= dlog.data[i,]$doe)
  met_df_year <- met_df_year %>% filter(doy >= site_row$dos)
  
  # If we predicted emergence before the observed planting date, we set the
  # AGDD to 0. 
  if (nrow(met_df_year) == 0) {
    dlog.data[i,]$AGDD <- 0
    next
  }
  
  AGDD_sum <- 0
  for (j in 1:nrow(met_df_year)) {
    max_temp <- met_df_year[j,]$mx2t$mx2t
    min_temp <- met_df_year[j,]$mn2t$mn2t
    
    if (is.na(max_temp) || is.na(min_temp)) {
      next
    }
    
    if (max_temp > 30) {
      max_temp = 30
    } else if (max_temp < 10) {
      max_temp = 10
    } 
    if (min_temp < 10) {
      min_temp = 10
    }
    gdd <- ((max_temp + min_temp)/2) - 10
    AGDD_sum <- AGDD_sum + gdd
  }
  dlog.data[i,]$AGDD <- AGDD_sum
}

# Remove the met. data variable and close the met. ncfile. 
rm(met_nc_data)
nc_close('inst/data/met_data.nc')

write.csv(dlog.data, file = "inst/dlog_data.csv", row.names = FALSE)

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

long_enough <- c()

for (i in 1:length(trainIDs)) {
  id = IDs[i]
  blah = temp %>% filter(ID == id)
  if (nrow(blah) > 20) {
    long_enough <- c(long_enough, id)
  }
}

sentinel_preds <- data.frame(
  ID = c(long_enough)
)
sentinel_preds$dos <- NA
sentinel_preds$doe <- NA

pb = txtProgressBar(0, 1, style=3)

for (i in 1:length(long_enough)) {
  setTxtProgressBar(pb, i/length(long_enough))

  site_id <- long_enough[i]
  sentinel_preds[i,]$ID = site_id
  
  site = train %>% filter(ID == site_id)
 # site_VI <- site$preds # VI time-series for the site. 
  
  # We will use MACD(5, 10, 5). (Gao et al. 2020)
  # I.e., S = 5 days ("fast" EMA), L = 10 days ("slow" EMA), 
  # K = 5 days ("signal"/"average" EMA)
  macd_t <- MACD_VI(site_VI, 5, 10, 10)
  site_preds$macd_t <- macd_t[,1]
  site_preds$ema_macd_c <- macd_t[,2]
  site_preds$macd_div_t <- site_preds$macd_t - site_preds$ema_macd_c
  
  # Settings (0.0, 0.0) for (MACD_threshold, MACD_div_threshold) are standard 
  # settings to track trend changes in stock prices. 
  # We adjust to use settings (0.1, 0.01) here (the fitted DL has uniform tails).
  macd_threshold <- 0.1
  macd_div_threshold <- 0.0001
  
  # The settings of 0.01 and 0.85 for NDVI at the green-up dates represent a 
  # wide range of acceptable green-up conditions since NDVI and EVI2 at green-up
  # time are very low.
  min_VI_greenup <- 0.01
  max_VI_greenup <- 0.90
  VI_increase_doys <- c()
  
  max_ndvi_doy <- site_preds[which.max(site_preds$preds),]$DOY
  
  # We check for momentum points with the various conditions described by 
  # Gao et al. (2020).
  # (ex:
  #   -- Is the NDVI at this doy within a reasonable range of NDVI values?, 
  #   -- Does the MACD begin to increase exponentially?, 
  #   -- Is the MACD within a reasonable range? 
  # )
  for (j in 11:max_ndvi_doy) {
    in_VI_range <-
      site_preds$preds[j] > min_VI_greenup &&
      site_preds$preds[j] < max_VI_greenup
    macd_div_x_intercept <-
      (
        site_preds$macd_div_t[j - 1] < macd_div_threshold &&
          site_preds$macd_div_t[j] > macd_div_threshold
      ) 
    
    macd_within_range <- macd_t[j] < macd_threshold
    if (in_VI_range && macd_div_x_intercept && macd_within_range) {
      VI_increase_doys <- c(VI_increase_doys, j)
    }
  }
  
  if (length(VI_increase_doys) == 0) {
    print(i)
    next
  }
  
  momentums <- data.frame(DOY = VI_increase_doys, NDVI = NA)
  for (k in 1:length(VI_increase_doys)) {
    site_row <- site_preds[VI_increase_doys[k],]
    momentums[k,]$NDVI = site_row$preds
  }
  
  ## Use a horizontal line to mark `macd_div_t` peak
  emerg_row <- site_preds[VI_increase_doys[1],]
  
  ## Use a horizontal line to mark actual planting date
  site_dos <- (temp %>% filter(ID == trainIDs[i]))[1,]$dos
  
  # Store 'doe' for the site (a.k.a. day of emergence).
  dlog.data[dlog.data$ID == site_id,]$doe <- emerg_row$DOY 
}

# --------------------------------------------------------------------------- #
# 12/08/2021 (Alex):

###############################################################################
# ??. Retrieve L8S2 (Landsat 8, Sentinel 2) data for training and test sites 
###############################################################################

# Retrieve Landsat 8 data frame: 
require(L8S2)
require(devtools)
require(remotes)
require(tidyverse)

L8S2.data = data.frame(ID = NULL,
                       Date = NULL,
                       NDVI = NULL, 
                       Satellite = NULL)

blah <- temp %>% filter(temp$Year == "2017")
IDs <- unique(blah$ID)

for (i in 41:94) {
  site_id <- IDs[i]
  site_id <- toString(site_id)
  print(site_id)
  
  this <- temp %>% filter(ID == site_id)
  lat <- round(this[1,]$Latitude, 7)
  lon <- round(this[1,]$Longitude, 7)
  year <- substr(this[1,]$Date, 1, 4)
  
  request_ID = paste("Site", site_id, sep="_")
  mysites <- data.frame(x= c(lon),  # lon.
                        y= c(lat),  # lat.
                        ID= c(request_ID) # Site ID
  )
  
  # This should take a few minutes:
  RS <- DownloadL8S2(mysites, '2017-04-01', '2018-12-31', Indices = c("NDVI"))
  
  # Clean/re-organize the `RS` data frame: 
  RS <- subset(RS, select = -c(...1, Index))
  site_col <- rep(site_id, nrow(RS))
  RS$ID <- site_col
  # Order the data by observation date. 
  RS <- RS[order(RS$date),]
  # Rename the "date" and "Value" columns:
  names(RS)[names(RS) == "date"] <- "Date"
  names(RS)[names(RS) == "Value"] <- "NDVI"
  # Reorder the columns of the `RS` data frame:
  RS <- subset(RS, select=c(3,2,4,1))
  RS <- RS %>% filter(substr(RS$Date,1,4) == "2017")
  #RS <- as.data.frame(RS)
  
  L8S2.data <- rbind(L8S2.data, RS)
}

# Write the `L8S2.data` data frame to a csv file.
write.csv(L8S2.data, "inst/data/L8S2_data.csv", row.names = FALSE)

# Read the `L8S2.data` data frame to a csv file.
L8S2.data <- read.csv(file = "inst/data/L8S2_data.csv")
# Convert the `L8S2.data` data column into a double (as.Date).
L8S2.data$Date <- as.Date(L8S2.data$Date)

## Clean L8S2 data:
ndvi <- c(L8S2.data$NDVI)
ndvi_mean <- mean(ndvi)
ndvi_sd <- sd(ndvi)
idx_remove <- c()
z_scores <- c()

# Remove any NDVI observations +/-3 SDs away from the mean:
for (i in 1:nrow(L8S2.data)) {
  obs <- L8S2.data[i,]$NDVI
  z_score <- abs((obs - ndvi_mean)/(ndvi_sd))
  z_scores <- c(z_score, z_scores)
  
  if (z_score > 3) {
    idx_remove <- c(i, idx_remove)
  }
}

to_remove <- c()
# Remove any sites with fewer than 7 NDVI observations:
for (i in 1:length(IDs)) {
  id <- IDs[i]
  temp <- L8S2.data %>% filter(ID == id)
  if (nrow(temp) < 7) {
    to_remove <- c(id, to_remove)
  }
}
# Remove empty date rows for ID _____:
# nrow(L8S2.data) = 10707
L8S2.data <- L8S2.data[-c(9297:9335),]
# nrow(L8S2.data) = 10668
# 68199
L8S2.data <- L8S2.data %>% filter(ID != 68199)
# nrow(L8S2.data) = 10663
# 68104
L8S2.data <- L8S2.data %>% filter(ID != 68104)
# nrow(L8S2.data) = 10658
# 68220
L8S2.data <- L8S2.data %>% filter(ID != 68220)
# nrow(L8S2.data) = 10652
# 57913
L8S2.data <- L8S2.data %>% filter(ID != 57913)
# nrow(L8S2.data) = 10646


