require(L8S2)
require(tidyverse)

remotes::install_github("AgronomicForecastingLab/L8S2/L8S2", dependencies=TRUE)

mysites <- data.frame(x= c(-88.20,-88.20),  #long
                      y= c(40.06,42.06),    #lat
                      ID= c('EnergyFarm',"RandomFarm") #Site ID
)

RS <- DownloadL8S2(mysites,  '2020-01-01', '2021-01-01', Indices = c("NDVI", "EVI", "SAVI", "OSAVI"))
