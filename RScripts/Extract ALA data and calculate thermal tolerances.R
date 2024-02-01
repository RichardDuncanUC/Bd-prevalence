
library(tidyverse)
library(lubridate)
library(galah)
library(terra)
library(sf)

# list of frog species
sp.list <- read.csv("./data/species list.csv")$species


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
  dplyr::select(species, lon, lat, year) %>%
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

# read in Anstisa alba locations from tidy prevalence data
all.dat <- read.csv("./data/Tidy Bd prevalence data.csv") |>
  glimpse()

loc <- all.dat |>
  filter(species == "Anstisia alba") |>
  dplyr::select(species, long, lat, year) |>
  rename(lon = long) |>
  glimpse()

dat <- bind_rows(dat, loc) |>
  arrange(species)

table(dat$species)

#-------------------------------------------------------------------------------
# read in gridded mean annual temperature data
# the tif file was created by running "Download silo data temperature.R"
# it is not in the data folder, so needs to be generated

all.temp <- rast("./data/gridded_mean_annual_temp_1889_2021.tif")
all.temp
names(all.temp)

# long-term overall mean annual temperature each cell
mean.mat <- mean(all.temp)

# extract temperature data
a <- terra::extract(mean.mat, cbind(dat$lon, dat$lat), cells = TRUE)
dat$temp <- a[, 2]
dat$cell <- a[, 1]
dat <- dat[!is.na(dat$temp), ]
glimpse(dat)

#-------------------------------------------------------------------------------
# get min, max and median temp for each species
# and temp difference between full and recent records (with recent defined by cutoff)
# thin temps by including only one value per cell
mt <- dat |>
  group_by(species, cell) |>
  summarise(temp = mean(temp, na.rm = T)) |>
  group_by(species) |>
  summarise(med.temp = median(temp, na.rm = T),
            min.temp = min(temp, na.rm = T),
            max.temp = max(temp, na.rm = T),
            n.cell = n()) |>
  glimpse()

write.csv(mt, "./data/Thermal optima.csv", row.names = F)


