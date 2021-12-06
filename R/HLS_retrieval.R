packages <- c('leaflet','rgdal','raster','jsonlite','sp','httr','rasterVis','ggplot2','magrittr', 'RColorBrewer','xml2','dygraphs','xts','lubridate','DT','rmarkdown', 'rprojroot')

# Identify missing (not installed) packages.
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]

# Install new (not installed) packages.
if(length(new.packages)) install.packages(new.packages, repos='http://cran.rstudio.com/') else print('All required packages are installed.')

invisible(lapply(packages, library, character.only = TRUE))

# Create an output directory if it doesn't exist.
wd <- rprojroot::find_rstudio_root_file()
outDir <- file.path(wd, "inst", "data", "HLS_data", fsep="/")
suppressWarnings(dir.create(outDir)) 

search_URL <- 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD/search'

# Now we need to query the HLS data by region and time period of interest.
# First we define our query parameters:

HLS_col <- list("HLSS30.v2.0", "HLSL30.v2.0")











