# Package workflow demo:

devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies=TRUE)

# Aggregate all adequate NDVI image into a x-y-n array.
dat <- aggregate_dates(9211, 7305, 24)