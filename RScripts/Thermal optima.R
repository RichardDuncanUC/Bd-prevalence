
library(tidyverse)
library(lubridate)
library(terra)
library(sf)

rm(list = ls())

# unlogit and scale functions
unlogit <- function(x) exp(x) / (1 + exp(x))
sc <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

# list of species
sp.list <- read.csv("./data/species list.csv")

################################################################################
# read in Australian amphibian records
dat <- read.csv("C:/Users/s429217/OneDrive/DATA/Niche contraction/amphibian/ALA amphibian data.csv") |>
  rename(species = binomial)
glimpse(dat)

dat <- filter(dat, species %in% sp.list$species)
table(dat$species)

# read in Anstisa alba locations
loc <- all.dat |>
  filter(species == "Anstisia alba") |>
  dplyr::select(species, long, lat, year) |>
  rename(lon = long)

dat <- bind_rows(dat, loc)
table(dat$species)

################################################################################
# read in frog shape files
load("c:/users/s429217/onedrive/data/frog traits/final data/Frog shape files.RData")

# filter to sp.list
range.dat$species <- paste(range.dat$genus, range.dat$species)
range.dat <- range.dat |>
  filter(species %in% sp.list) |>
  arrange(species)
range.dat

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

#-------------------------------------------------------------------------------
# read in gridded mean annual temperature data
all.temp <- rast("c:/users/s429217/onedrive/data/gridded rainfall data/gridded_mean_annual_temp_1889_2021.tif")
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

write.csv(mt, "c:/users/s429217/onedrive/data/chytrid/tidy data/Thermal optima.csv", row.names = F)
