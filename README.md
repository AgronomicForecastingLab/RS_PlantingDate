# RS_PlantingDate

R package providing functions to process NDVI time series and estimate corn
planting date at a given spatial resolution.

## Languages and key dependencies

- R (package: `RSPlantingDate`)
- A small Python helper script (`inst/python/hls_retrieval.py`) for HLS retrieval
- R imports: `jsonlite`, `yaml`, `dplyr`, `R6`, `ggplot2`
- R depends: `R6`, `dplyr`, `purrr`
- Workflows additionally use: `tidyverse`, `phenopix`, `greenbrown`, `zoo`,
  `ncdf4`, `ncmeta`, `stars`, `units`, `ecmwfr`, `RNetCDF`, `raster`, `splines`,
  `GGally`, `devtools`

## Install

From R, install directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("AgronomicForecastingLab/RS_PlantingDate", dependencies = TRUE)
```

## Quickstart

A minimal end-to-end example using the included sample data:

```r
library(RSPlantingDate)
library(dplyr)

# Load bundled cleaned NDVI data
load(system.file("data/organized/cleaned_corn_data.Rdata", package = "RSPlantingDate"))

# Pick one site and fit a double-logistic model to its NDVI time series
ids <- unique(cleaned_corn_data$ID)
site <- cleaned_corn_data %>% filter(ID == ids[1])

fit <- fit_double_logistic(x = site$NDVI, t = site$DOY)

# Predict NDVI for a sequence of days using the fitted parameters
preds <- predictDLOG(params = as.numeric(fit[1, ]), t = seq(60, 280, by = 1))
head(preds)
```

For a fuller pipeline (data cleaning -> double-logistic fit -> emergence ->
planting date), see the workflow scripts under `inst/`:

- `inst/full_workflow.R` -- double-logistic workflow
- `inst/quadratic_workflow.R` -- Venancio quadratic-equation workflow
- `inst/data_organization_workflow.R` -- assembling Beck's NDVI data
- `inst/emergence_tests.R`, `inst/draw_sample.R` -- analysis helpers
- `inst/workflow_overview.txt` -- prose overview of the full pipeline

## What's inside

- `R/aggregate_dates.R` -- aggregate per-date NDVI rasters into an `x*y*n` array
- `R/aggregate_met.R` -- aggregate hourly ERA5 met data to daily
- `R/clean_monotonic.R` -- remove non-monotonic NDVI points per site
- `R/clean_spline.R` -- spline-smooth NDVI time series
- `R/clean_WLS.R` -- weighted least-squares moving-window cleaning
- `R/compare_cleaned_plots.R` -- ggplot of before/after cleaning
- `R/dlog.R` -- RMSE objective for double-logistic optimization
- `R/predictDLOG.R` -- predict NDVI from fitted double-logistic parameters
- `R/fit_double_logistic.R` -- fit a double logistic per site (Beck's eq.)
- `R/fit_normal.R` -- fit a normal distribution to NDVI vs. DOY
- `R/EMA_VI.R`, `R/MACD_VI.R` -- EMA and MACD filters on NDVI series
- `R/get_site_met.R` -- extract ERA5 min/max temperature for a site
- `R/predict_DOS.R` -- predict day-of-sowing per pixel via Venancio equations
- `R/smooth_3x3.R` -- 3x3 smoothing window over a DOS matrix
- `inst/data/` -- bundled CSV / Rdata sample inputs (Beck's, L8/S2, met, etc.)
- `inst/python/hls_retrieval.py` -- HLS retrieval helper
- `man/` -- generated roxygen `.Rd` documentation
- `examples/` -- minimal runnable example
- `DESCRIPTION`, `NAMESPACE` -- R package metadata

## Citation

No `CITATION.cff` or paper reference is currently present in the repository.
Authors listed in `DESCRIPTION`: Alexandra Jakucewicz, Marissa Kivi (maintainer:
Hamze Dokoohaki, `hamzed@illinois.edu`).

## License

The `DESCRIPTION` file currently has the placeholder `License: What license is
it under?` and no `LICENSE` file is present in the repository. Please contact
the maintainer regarding terms of use until a license is added.
