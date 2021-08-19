# Package workflow demo:

devtools::install_github("AgronomicForecastingLab/noaaGEFSpoint", dependencies=TRUE)
sites <- read.csv("sites.csv")
output_directory <- "Data"
site_file <- system.file("extdata", "noaa_download_site_list.csv", package = "noaaGEFSpoint")
neon_sites <- read_csv(site_file)[1:2,]

# ERA5 Dataset: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

# Set ERA5 dataset download key
wf_set_key(user= "67047", key = "185f8f5d-4c98-46ed-830c-d76a6a0758fa", service = "cds")

# Request ERA5 data
era5_data <- era5_ncfile_request("2019-12-01", "2019-12-03", sites[1:4,])
# Reformat ERA5 data
era5_data_refmt <- reformat_ERA(era5_data)

# Request NOAA GEFS data
gefs_data <- noaa_gefs_request(neon_sites)
# Reformat NOAA GEFS data
gefs_data_refmt <- reformat_NOAA(gefs_data)
# Aggregate NOAA GEFS data
gefs_data_agg <- noaa_gefs_aggregate(gefs_data_refmt)

# Bind ERA5 and NOAA data together
noaa_era_bind <- bind_NOAA_ERA(gefs_data_agg, era5_data_refmt)