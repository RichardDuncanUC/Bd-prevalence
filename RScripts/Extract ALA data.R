
library(tidyverse)
library(lubridate)
library(galah)
library(terra)
library(sf)

# list of frog species
sp.list <- read.csv("./data/species list.csv")


# download records from ALA and filter
galah_config(email = "richard.duncan@canberra.edu.au", atlas = "Australia")

# total number of records for all species
galah_call() %>%
  galah_identify(sp.list$species) %>% 
  atlas_counts()

# extract the records
dat <- galah_call() %>%
  galah_identify(sp.list$species) %>%
  galah_select(eventID, scientificName, decimalLatitude, decimalLongitude, eventDate, stateProvince, 
               coordinateUncertaintyInMeters, establishmentMeans, basisOfRecord, occurrenceStatus, 
               assertions, outlierLayerCount) %>%
  atlas_occurrences()

head(dat)
glimpse(dat)

# filter, removing uncertain or unwanted records
dat <- dat %>%
  mutate(date = ymd_hms(eventDate),
         year = year(date)) %>%
  filter(!(assertions %in% c("PRESUMED_SWAPPED_COORDINATE", "COORDINATES_CENTRE_OF_STATEPROVINCE", 
                    "COORDINATES_CENTRE_OF_COUNTRY", "PRESUMED_NEGATED_LATITUDE", 
                    "PRESUMED_NEGATED_LONGITUDE"))) %>%
  filter(basisOfRecord != "FOSSIL_SPECIMEN") %>%
  rename(lon = decimalLongitude,
         lat = decimalLatitude,
         species = scientificName) %>%
  dplyr::select(binomial, lon, lat, year) %>%
  drop_na()

glimpse(dat)
  
#-------------------------------------------------------------------------------
# clip to frog shape files
load("./data/Frog shape files.RData")

sf_use_s2(FALSE)

# subset to points within range map polygons
in.range <- list()
for(i in 1:length(sp.list)) {
  temp <- filter(dat, species == sp.list[i])
  temp$loc <- paste(temp$lon, temp$lat)
  xy <- data.frame(lon = temp$lon, lat = temp$lat)
  pnts <- st_as_sf(xy, crs = crs(range.dat), coords = c("lon", "lat"))
  in_pnts <- st_filter(pnts, range.dat$geometry[range.dat$species == sp.list[i]])
  in.xy <- st_coordinates(in_pnts$geometry)
  keep <- paste(in.xy[, 1], in.xy[, 2])
  in.range[[i]] <- filter(temp, loc %in% keep)
}

dat <- in.range[[1]]
for(i in 2:length(sp.list)) {
  dat <- bind_rows(dat, in.range[[i]])
}




# read in shapefile of Australian mainland
mainvec <- vect("c:/users/s429217/onedrive/data/niche contraction/mammals/data/Australia shape vector.gpkg")
mainvec <- st_as_sf(mainvec)

# extract lats and longs
xy <- data.frame(cbind(dat$lon, dat$lat))
xy <- st_as_sf(xy, coords = 1:2, crs = crs(mainvec))

int <- st_filter(xy, mainvec)

keep.xy <- st_coordinates(int$geometry)

keep <- which(paste(dat$lon, dat$lat) %in% paste(keep.xy[, 1], keep.xy[, 2]))

dat <- dat[keep, ]
summary(dat$lon)
summary(dat$lat)

dat <- dat %>%
  arrange(binomial)

glimpse(dat)

write.csv(dat, "C:/Users/s429217/OneDrive/DATA/Niche contraction/amphibian/ALA amphibian data.csv", row.names = F)

