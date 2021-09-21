require(tidyverse)
require(GGally)
require(phenopix)
require(zoo)
require(greenbrown)
require(RSPlantingDate)

# Data cleaning and fitting workflow with the double logistic function

rm(list=ls())
set.seed(102396)

##############################
# Prepare the data 
##############################

# The Beck's Hybrids dataset should be named 'cleanedData.csv'.
# We need to remove duplicates of data points for the correct function of our workflow 
temp <- read.csv('inst/data/cleanedData.csv') %>% 
  mutate(Date = as.Date(Date)) %>%
  distinct(ID, Date,.keep_all = TRUE)

# Remove any sites with late DOS values or with DOS values after harvest (dataset errors).
temp <- temp %>% filter(dos < 182, dos < doh)

# Remove any non-corn sites (corn is 'Family 1').
temp = temp %>% filter(Family == 1)

# Store a version of `temp` before data cleaning to use in later comparison. 
before = temp
IDs = unique(temp$ID)

# Apply spline cleaning function here 
if (!file.exists('inst/data/cleaned_corn_data.Rdata')){
  temp <- clean_WLS(orig_df = temp)
  
  # Remove any sites with fewer than 5 data points.
  finalIDs = c()
  for (i in 1:length(IDs)) {
    this <- temp %>% dplyr::filter(ID == IDs[i])
    if (nrow(this) >= 5)
      finalIDs = c(finalIDs, IDs[i])
  }
  temp <- temp %>% dplyr::filter(ID %in% finalIDs)
  
  save(temp, file = 'inst/data/cleaned_corn_data.Rdata')
} else{
  load('inst/data/cleaned_corn_data.Rdata')
}

# Generate 12 plots comparing uncleaned and cleaned data.
samples = sample(IDs, 12)
comp_plot <- compare_cleaned_plots(samples, before, temp)
comp_plot

cornIDs = unique(temp$ID)

if (!file.exists('inst/data/test_train_data.Rdata')){
  
  # Randomly select training and test data.
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
}else{
  load('inst/data/test_train_data.Rdata')
}


##############################
# Model fitting :: double logistic
##############################

trainIDs = unique(train$ID)

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

# Now we set ranges on our different parameters 
# start of season: XX < sos < XX
# end of season: XX < eos < XX
# maximum NDVI: > 0.6
# winter NDVI: < 0.3
# day of max NDVI: XX < maxDOY < XX

ids = sample(unique(dlog.data$ID), 12)
before = before %>% filter(ID %in% ids)
before$check = 'before'
after = temp %>% filter(ID %in% ids)
after$check = 'after'

fitted = data.frame(ID = NULL, DOY = NULL, pred = NULL)
for (i in 1:length(ids)){
  pars = as.vector(dlog.data %>% filter(ID == ids[i]) %>% dplyr::select(mn, mx, sos, rsp, eos, rau))
  fitted = rbind(fitted, 
                 data.frame(ID = rep(ids[i], 365),
                            DOY = 1:365,
                            preds = predictDLOG(pars, 1:365)))
}


check = rbind(before,after) %>% mutate(Date = as.Date(Date),
                                       DOY = as.numeric(Date - as.Date('2016-12-31')))

ggplot(check) +
  geom_point(aes(x = DOY, y = NDVI, color = check)) + 
  geom_line(data = check %>% filter(check == 'after'), aes(x = DOY, y = NDVI,  linetype = 'cleaned'),color = 'red') +
  geom_line(data = fitted, aes(x = DOY, y = preds, linetype = 'fitted'))+
  scale_color_manual(values = c('before' = 'black', 'after' = 'red'))+
  scale_linetype_manual(values = c('cleaned'='dashed', 'fitted'='solid')) +
  facet_wrap(~ID, nrow = 4) 


##############################
# Employing the fitted model 
##############################
### THIS IS WHERE MY UPDATES ENDED (9/21)

# Consider the relationship between variables and DOS to check for potential use in predicting planting date
ggpairs(dlog.data %>% dplyr::select(mn,mx,sos,rsp,eos,rau,lat,lon,maxDOY,dos))
ggpairs(dlog.data %>% dplyr::select(dos, sos, eos, lat, maxDOY))

int <- lm(dos ~ 1, data=dlog.data)
all <- lm(dos ~ eos*lat*maxDOY*lon*sos, data=dlog.data)

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
                       lat = NULL, 
                       lon = NULL)

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
all <- lm(dos ~ lat*mu*sigma*k*lon, data=norm.data)

# Perform forward stepwise regression
normMod <- step(int, direction='forward', scope=formula(all), trace=0)

################################################
#### Let's check out the model residuals to figure out what's going on 
################################################

plotDat = data.frame(norm = normMod$residuals, 
                     ID = norm.data$ID,
                     dos = norm.data$dos) %>%
  left_join(data.frame(dlog = dlogMod$residuals, 
                       ID = dlog.data$ID),
            by = 'ID') %>%
  reshape2::melt(id = c('ID','dos')) %>%
  left_join(temp %>% dplyr::select(ID, Latitude, Longitude) %>% distinct(.keep_all=TRUE),
            by = 'ID')

# residuals vs. dos
ggplot(plotDat) + 
  geom_point(aes(x = dos, y = value, color = variable), alpha = 0.75) + 
  geom_smooth(method = 'lm', aes(x=dos, y=value, group=variable), color = 'black',se = F) +
  labs(title = 'Model Residuals vs. Observed Planting Date',
       x = 'Observed Planting Date',
       y = 'Residuals (Observed-Predicted)',
       color = 'Model') +
  theme_bw() +
  theme(legend.position = 'top', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) 

# residuals in space
# first classify them into size classes
categorize = function(res){
  if (is.na(res)) return(NA)
  if(res < -20) cat = '< -20'
  if(res >= -20 & res < 0) cat = '-20-0'
  if(res >= 0 & res < 20) cat = '0-20'
  if(res >= 20 & res < 40) cat = '20-40'
  if(res >= 40) cat = '> 40'
  return(cat)
}
plotDat$cat = sapply(plotDat$value,categorize)
plotDat$cat = factor(plotDat$cat, levels = c('< -20','-20-0','0-20','20-40','> 40'))

# looking at how residual direction/magnitude changes spatially 
ggplot(plotDat %>% filter(!is.na(value))) + 
  geom_point(aes(x = Longitude, y = Latitude, color = cat, size = value), alpha = 0.75) + 
  labs(title = 'Spatial Pattern in Residuals',
       x = 'Latitude',
       y = 'Longitude',
       color = 'Residual Size Class') +
  facet_wrap(~variable) + 
  theme_bw() +
  guides(size = "none") + 
  theme(legend.position = 'right', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) +
  scale_color_manual(values = c('< -20'='blue', '-20-0'='deepskyblue', '0-20'='lightblue',
                                '20-40'='lightgoldenrod2','> 40'='goldenrod4'))

categorize2 = function(dos){
  if (is.na(dos)) return(NA)
  if(dos < 90) cat = 1
  if(dos >= 90 & dos < 110) cat = 2
  if(dos >= 110 & dos < 130) cat = 3
  if(dos >= 130 & dos < 150) cat = 4
  if(dos >= 150) cat = 5
  return(cat)
}
plotDat$cat2 = sapply(plotDat$dos,categorize2)
plotDat$cat2 = factor(plotDat$cat2, levels = c(1:5))

# considering how planting date might be related to regions
ggplot(plotDat %>% filter(!is.na(value))) + 
  geom_point(aes(x = Longitude, y = Latitude, color = as.factor(cat2)), alpha = 0.75) + 
  labs(title = 'Regional Pattern in Day of Sowing',
       x = 'Latitude',
       y = 'Longitude',
       color = 'Day of Sowing') +
  facet_wrap(~variable) + 
  theme_bw() +
  theme(legend.position = 'right', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) 

ggplot(plotDat %>% filter(!is.na(value))) + 
  geom_bar(aes(x = cat2, fill = as.factor(cat)), alpha = 0.75) + 
  labs(title = 'Residual Category vs. Day of Sowing Category',
       x = 'Day of Sowing Grouping',
       y = 'Count',
       fill = 'Residual Category') +
  facet_wrap(~variable) + 
  theme_bw() +
  theme(legend.position = 'right', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) 

# consider how other model parameters might relate to residuals => no relationship here really 
ggplot(plotDat %>% filter(!is.na(value)) %>% 
         dplyr::select(-value, -variable) %>%
         distinct(ID, .keep_all = TRUE) %>%
         left_join(dlog.data %>% dplyr::select(ID, mn, mx, sos, rsp, eos, rau, maxDOY), by = c('ID')) %>%
         left_join(norm.data %>% dplyr::select(ID, mu, sigma, k), by = c('ID')) %>%
         reshape2::melt(id = c('ID','dos','Latitude','Longitude','cat','cat2'))) + 
  geom_histogram(aes(x = value, fill = cat), alpha = 0.75) + 
  labs(title = 'Model Parameters and Residual Size',
       x = 'Parameter Value',
       y = 'Count',
       fill = 'Residual Group') +
  facet_wrap(~variable, scales = 'free') + 
  theme_bw() +
  theme(legend.position = 'right', text = element_text(size = 20, family = 'serif'),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 90, family = 'serif', size = 20),
        axis.title.x = element_text(angle = 0, family = 'serif', size = 20),
        strip.text.y = element_text(family = 'serif', size = 20),
        panel.grid.minor = element_line(color = 'gray84'), 
        panel.grid.major = element_line(color = 'gray84')) 


################################################
#### Testing the two new functions at each data point 
################################################

# Organize test data
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

