# Minimal runnable example: fit a double-logistic model to one site's NDVI
# series from the bundled cleaned data, then predict NDVI over a DOY range.
#
# Run from R after installing the package:
#   devtools::install_github("AgronomicForecastingLab/RS_PlantingDate",
#                            dependencies = TRUE)
#   source("examples/minimal_double_logistic.R")

library(RSPlantingDate)
library(dplyr)

# Load bundled cleaned NDVI data (object: cleaned_corn_data)
load(system.file("data/organized/cleaned_corn_data.Rdata",
                 package = "RSPlantingDate"))

# Pick the first site and fit a double-logistic model
ids  <- unique(cleaned_corn_data$ID)
site <- cleaned_corn_data %>% filter(ID == ids[1])

fit <- fit_double_logistic(x = site$NDVI, t = site$DOY)
print(fit)

# Predict NDVI across the growing season
doy_grid <- seq(60, 280, by = 1)
preds    <- predictDLOG(params = as.numeric(fit[1, ]), t = doy_grid)

plot(doy_grid, preds, type = "l",
     xlab = "DOY", ylab = "Predicted NDVI",
     main = paste("Double-logistic fit, site", ids[1]))
points(site$DOY, site$NDVI, pch = 19, cex = 0.6)
