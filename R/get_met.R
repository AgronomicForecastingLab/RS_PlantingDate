#' Retrieves hourly total precipitation, maximum temperature, minimum 
#' temperature, soil temperature, and soil moisture data from the CDS reanalysis
#' data product 'ERA5 hourly data on single levels from 1979 to present'. 
#'
#' @param lat Site latitude coordinate.
#' @param lon Site longitude coordinate.
#' @param date blaaah
#' @return A data frame containing cumulative GDD and precipitation for each month of the year
#' @export
get_met <- function(lat, lon, date) {
  request <- list(
    product_type = "reanalysis",
    format = "netcdf",
    variable = c(
      "maximum_2m_temperature_since_previous_post_processing",
      "minimum_2m_temperature_since_previous_post_processing",
      "soil_temperature_level_1",
      "total_precipitation",
      "volumetric_soil_water_layer_1"
    ),
    year = "2017",
    month = "07",
    day = "10",
    time = c(
      "00:00",
      "01:00",
      "02:00",
      "03:00",
      "04:00",
      "05:00",
      "06:00",
      "07:00",
      "08:00",
      "09:00",
      "10:00",
      "11:00",
      "12:00",
      "13:00",
      "14:00",
      "15:00",
      "16:00",
      "17:00",
      "18:00",
      "19:00",
      "20:00",
      "21:00",
      "22:00",
      "23:00"
    ),
    dataset_short_name = "reanalysis-era5-single-levels",
    target = "download.nc"
  )
  return(request)
}