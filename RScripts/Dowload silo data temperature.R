
# see https://www.longpaddock.qld.gov.au/silo/gridded-data/
# for download instructions

library(ncdf4) # package for netcdf manipulation
library(terra) 

getOption('timeout')
options(timeout = 500)

# web address for monthly gridded rainfall data
silo.min <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/min_temp/"
silo.max <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/max_temp/"

yr.list <- 1889:2021

# download first year
yr <- yr.list[1]
file.min <- paste0(silo.min, yr, ".min_temp.nc")
file.max <- paste0(silo.max, yr, ".max_temp.nc")
download.file(file.min, destfile = file.path("c:/temp/sample_min.nc"), mode = "wb")
download.file(file.max, destfile = file.path("c:/temp/sample_max.nc"), mode = "wb")

temp.min <- terra::rast("c:/temp/sample_min.nc")
temp.max <- terra::rast("c:/temp/sample_max.nc")

mean.temp <- mean((temp.min + temp.max) / 2)
names(mean.temp) <- paste0("yr_", yr) 

# download the other years
for(j in 2:length(yr.list)) {
  yr <- yr.list[j]
  file.min <- paste0(silo.min, yr, ".min_temp.nc")
  file.max <- paste0(silo.max, yr, ".max_temp.nc")
  download.file(file.min, destfile = file.path("c:/temp/sample_min.nc"), mode = "wb")
  download.file(file.max, destfile = file.path("c:/temp/sample_max.nc"), mode = "wb")
  
  temp.min <- terra::rast("c:/temp/sample_min.nc")
  temp.max <- terra::rast("c:/temp/sample_max.nc")
  
  temp <- mean((temp.min + temp.max) / 2)
  names(temp) <- paste0("yr_", yr) 
  
  mean.temp <- c(mean.temp, temp)
}

writeRaster(mean.temp, "./data/gridded_mean_annual_temp_1889_2021.tif", overwrite = T)


