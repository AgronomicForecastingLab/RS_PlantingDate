# Package workflow demo:

devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies=TRUE)

# Aggregate all adequate NDVI image into a x-y-n array.
dat <- aggregate_dates(9211, 7305, 24)

# Dates list (24)
dates = as.Date(
  c(
    '2017-03-15',
    '2017-03-22',
    '2017-04-24',
    '2017-05-14',
    '2017-06-03',
    '2017-06-10',
    '2017-06-23',
    '2017-07-08',
    '2017-07-13',
    '2017-07-25',
    '2017-07-28',
    '2017-08-02',
    '2017-08-17',
    '2017-08-22',
    '2017-08-29',
    '2017-09-11',
    '2017-09-23',
    '2017-09-26',
    '2017-10-01',
    '2017-10-08',
    '2017-10-13',
    '2017-10-16',
    '2017-10-21',
    '2017-11-10'
  )
)

# List of days since planting at each date in `dates`.
doys = dates - as.Date('2016-12-31')

# Functions from paper for predicting sowing date with two equations.
fA = function(dos, t, N, minNDVI) {
  pred = -0.00013 * (t - dos) ^ 2 + 0.01787 * (t - dos) + minNDVI
  pred[pred < minNDVI] = minNDVI
  sqError = (N - pred) ^ 2
  return(sqrt(mean(sqError, na.rm=TRUE)))
}
A = function(dos, t, minNDVI) {
  pred = -0.00013 * (t - dos) ^ 2 + 0.01787 * (t - dos) + minNDVI
  pred[pred < minNDVI] = minNDVI
  return(pred)
}
fB = function(dos, t, N, minNDVI) {
  pred = -0.00014 * (t - dos) ^ 2 + 0.01931 * (t - dos) + minNDVI
  pred[pred < minNDVI] = minNDVI
  sqError = (N - pred) ^ 2
  return(sqrt(mean(sqError, na.rm=TRUE)))
}
B = function(dos, t, minNDVI) {
  pred = -0.00014 * (t - dos) ^ 2 + 0.01931 * (t - dos) + minNDVI
  pred[pred < minNDVI] = minNDVI
  return(pred)
}

dos.mat <- predict_DOS(dat, doys, fA, fB, 91, 156)




