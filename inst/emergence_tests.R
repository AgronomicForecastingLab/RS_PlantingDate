require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)
require(RSPlantingDate)

## This workflow is doing tests for detecting and significance of emergence for DOS

rm(list=ls())
################################################
#### Load prepared data
################################################

load('inst/data/monotonic_corn_data.Rdata')
load('inst/data/test_train_data.Rdata')

trainIDs = unique(train$ID)

################################################
#### Load prepared data
################################################

fit_double_logistic <- function(x, t) {
  
  # Fit the model based on Beck's et al. equation.
  mod = try(FitDoubleLog(
    x = x,
    t = t,
    weighting = TRUE,
    plot = FALSE
  ),
  silent = T)
  
  # If we aren't able to fit the model, we try 1x more, 
  # then skip if it's still bad.
  if (is.na(mod[2])) {
    mod = try(FitDoubleLog(
      x = x,
      t = t,
      weighting = TRUE,
      plot = FALSE
    ),
    silent = T)
    if (is.na(mod[2])) {
      return(NA)
    }
  }
  
  # Get the fitted model parameters.
  mod = as.vector(mod$params)
  
  # Get the predicted values for the entire year from the fitted dlog model.
  preds = predictDLOG(mod, c(0:365))
  slope = preds - as.numeric(mod[1])
  emerge = which(slope > 0.005)[1]
  
  # Calculate the departure from the winter NDVI
  # (over-time, to help us determine when the crop likely emerged).
  #diffs = preds - as.numeric(mod[1])

  # Saves the data for this plot in the 'dlog.data' data frame.
  res = data.frame(
    ID = this$ID[1],
    mn =  as.numeric(mod[1]),
    mx =  as.numeric(mod[2]),
    sos =  as.numeric(mod[3]),
    rsp =  as.numeric(mod[4]),
    eos =  as.numeric(mod[5]),
    rau =  as.numeric(mod[6]),
    dos = this$dos[1],
    lat = this$Latitude[1],
    lon = this$Longitude[1],
    maxDOY = which.max(preds), 
    emerge = emerge
  )
  return(res)
}

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
  lon = NULL,
  # longitude
  maxDOY = NULL, 
  # day of year where NDVI is at maximum following model
  emerge = NULL
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
dlog.data = dlog.data %>% filter(eos > 150, eos < 335)
ggplot(dlog.data) + geom_point(aes(x = emerge, y = lat))

# Consider the relationship between variables and DOS to check for potential use in predicting planting date
ggpairs(dlog.data %>% dplyr::select(mn,mx,sos,rsp,eos,rau,lat,lon,maxDOY,dos))
ggpairs(dlog.data %>% dplyr::select(dos, sos, eos, lat, maxDOY))

int <- lm(dos ~ 1, data=dlog.data)
all <- lm(dos ~ eos*lat*maxDOY*lon*sos, data=dlog.data)

# Perform forward stepwise regression
dlogMod <- step(int, direction='forward', scope=formula(all), trace=0)