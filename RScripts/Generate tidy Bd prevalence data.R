
library(tidyverse)
library(lubridate)
library(terra)
library(sf)

rm(list = ls())

#--------------------------------------------------------------------------------------------------------------
# read in Murray et al. 2010 data
# downloaded from https://figshare.com/collections/The_distribution_and_host_range_of_the_pandemic_disease_chytridiomycosis_in_Australia_spanning_surveys_from_1956_2007/3302961
# spans 1956-2007
dat <- read.csv("./data/Chytridiomycosis_data_1956_2007.csv", strip.white = T,
                na.strings = "-9999") |>
  filter(!is.na(Individuals)) |>
  filter(!is.na(Indivs_positive)) |>
  filter(!is.na(Latitude)) |>
  filter(!is.na(Longitude)) |>
  filter(Species != "unknown") |>
  filter(substr(Species, 1, 2) != "sp") |>
  glimpse()

table(dat$Species)
table(dat$Dead_or_sick)

# unique locations
dat$loc <- paste(dat$Longitude, dat$Latitude)

# remove dead or sick individuals
dat <- filter(dat, Dead_or_sick != "noted")

# sum by species, location, year
dat.loc <- dat |>
  mutate(Year = ifelse(is.na(Year), 2020, Year)) |>
  group_by(Species, loc, Year) |>
  summarise(lat = mean(Latitude),
            long = mean(Longitude),
            n.ind = sum(Individuals),
            n.pos = sum(Indivs_positive)) |>
  rename(species = Species) |>
  mutate(source = "Murray et al 2010",
         month = NA) |>
  filter(Year != 2020) |>
  rename(year = Year)

# fix species names
dat.loc$species[dat.loc$species == "Litoria lesueuri sensu lato"] <- "Litoria lesueuri"

glimpse(dat.loc)

table(is.na(dat.loc$month))

#---------------------------------------------------------------------------------------------------------------
# read in Amphibian disease portal data
# downloaded from https://amphibiandisease.org/query/
# specifying disease = Bd, for a box bounding Australia
dat1 <- read.csv("./data/samples_output.csv", strip.white = T) |>
  dplyr::select(eventID, materialSampleID, organismRemarks, individualCount, genus, specificEpithet, 
         infraspecificEpithet, projectId) |>
  mutate(genus = ifelse(genus == "crinia", "Crinia", genus),
         genus = ifelse(genus == "litoria", "Litoria", genus))
  
dat2 <- read.csv("./data/diagnostics_output.csv", strip.white = T) |>
  dplyr::select(diagnosticID, materialSampleID, sampleType, diseaseTested, diseaseTestedPositiveCount, 
         diseaseDetected, testMethod, projectId)

dat3 <- read.csv("./data/events_output.csv", strip.white = T) |>
  dplyr::select(eventID, yearCollected, monthCollected, dayCollected, locationID, decimalLatitude, decimalLongitude, maximumElevationInMeters,
         coordinateUncertaintyInMeters, projectId) |>
  filter(yearCollected != "unknown") |>
  mutate(year = as.numeric(yearCollected))

glimpse(dat1)
glimpse(dat2)
glimpse(dat3)

# join the three files
datp <- full_join(dat1, dat2) |>
  full_join(dat3)

glimpse(datp)

table(datp$diseaseDetected)
table(datp$individualCount)

# filter to exclude uncertain records
datp <- datp |>
  filter(diseaseDetected != "UNKNOWN") |>
  filter(organismRemarks != "Captive") |>
  filter(genus != "Unknown") |>
  filter(!(specificEpithet %in% c("sp.", "spp."))) |>
  filter(!is.na(decimalLatitude)) |>
  filter(!is.na(decimalLongitude)) |>
  mutate(individualCount = as.numeric(ifelse(individualCount == "<unk>", "", individualCount)),
         diseaseDetected = ifelse(diseaseDetected == "FALSE", "false", diseaseDetected),
         diseaseDetected = ifelse(diseaseDetected == "TRUE", "true", diseaseDetected),
         yearCollected = as.numeric(ifelse(yearCollected == "unknown", NA, yearCollected)),
         spp = paste(genus, specificEpithet, infraspecificEpithet),
         spp = trimws(spp)) 

# if individualCount is unknown, then treat as one observation using disease detected as outcome
unk <- ifelse(is.na(datp$individualCount), "no_count", "count")
table(unk)

datp$individualCount[unk == "no_count"] <- 1
datp$diseaseTestedPositiveCount[unk == "no_count"] <- ifelse(datp$diseaseDetected[unk == "no_count"] == "true", 1, 0)

glimpse(datp)
table(datp$year)
table(datp$monthCollected)

# convert to same form as Murray data
datp$loc = paste(datp$decimalLongitude, datp$decimalLatitude)

datp.loc <- datp |>
  group_by(spp, loc, year, monthCollected) |>
  summarise(lat = mean(decimalLatitude), 
            long = mean(decimalLongitude),
            n.ind = sum(individualCount),
            n.pos = sum(diseaseTestedPositiveCount)) |>
  rename(species = spp) |>
  mutate(source = "Amphibian disease portal") |>
  filter(species != "Litoria wilcoxii/jungguy") |>
  rename(month = monthCollected)

table(is.na(datp.loc$month))

#-------------------------------------------------------------------------------
glimpse(datp.loc)
glimpse(dat.loc)

# fix some names
datp.loc$species[datp.loc$species == "Litoria verreauxii verreauxii"] <- "Litoria verreauxii"

# merge the two data files
all.dat <- bind_rows(dat.loc, datp.loc)
glimpse(all.dat)

# check n.ind >= n.pos
all.dat$n.ind <- ifelse(all.dat$n.ind < all.dat$n.pos, all.dat$n.pos, all.dat$n.ind)

glimpse(all.dat)
table(all.dat$year, exclude = NULL)

#-------------------------------------------------------------------------------
# update nomenclature
all.dat$species[all.dat$species == "Geocrinia alba"] <- "Anstisia alba"
all.dat$species[all.dat$species == "Geocrinia rosea"] <- "Anstisia rosea"
all.dat$species[all.dat$species == "Geocrinia vitellina"] <- "Anstisia vitellina"
all.dat$species[all.dat$species == "Hylarana daemeli"] <- "Papurana daemeli"
all.dat$species[all.dat$species == "Lechriodus fletcheri"] <- "Platyplectrum fletcheri"
all.dat$species[all.dat$species == "Litoria alboguttata"] <- "Cyclorana alboguttata"
all.dat$species[all.dat$species == "Litoria andirrmalin"] <- "Litoria andiirrmalin"
all.dat$species[all.dat$species == "Litoria australis"] <- "Cyclorana australis"
all.dat$species[all.dat$species == "Litoria brevipes"] <- "Cyclorana brevipes"
all.dat$species[all.dat$species == "Litoria Brevipes"] <- "Cyclorana brevipes"
all.dat$species[all.dat$species == "Litoria genimaculata"] <- "Litoria serrata"
all.dat$species[all.dat$species == "Litoria novaehollandiae"] <- "Cyclorana novaehollandiae"
all.dat$species[all.dat$species == "Neobatrachus sudelli"] <- "Neobatrachus sudellae"
all.dat$species[all.dat$species == "Nyctimystes dayi"] <- "Litoria dayi"
all.dat$species[all.dat$species == "Litoria burrowsi"] <- "Litoria burrowsae"
all.dat$species[all.dat$species == "Litoria lesueurii"] <- "Litoria lesueuri"
all.dat$species[all.dat$species == "Litoria manya"] <- "Cyclorana manya"
all.dat$species[all.dat$species == "Litoria paraewingii"] <- "Litoria paraewingi"
all.dat$species[all.dat$species == "Litoria platycephala"] <- "Cyclorana platycephalus"
all.dat$species[all.dat$species == "Litoria verreauxii alpina"] <- "Litoria verreauxii"
all.dat$species[all.dat$species == "Litoria verrucosa"] <- "Cyclorana verrucosa"
all.dat$species[all.dat$species == "Litoria watjulemensis"] <- "Litoria watjulumensis"

# remove unknown species
all.dat <- all.dat |>
  filter(species != "Crinia destructor")

#-------------------------------------------------------------------------------
# deal with duplicate data
# compare location/years
# first, create new location/years
all.dat$long <- round(all.dat$long, 3)
all.dat$lat <- round(all.dat$lat, 3)
all.dat$loc.year <- paste(all.dat$long, all.dat$lat, all.dat$year, sep = ",")
all.dat$loc <- paste(all.dat$long, all.dat$lat, sep = ",")

#--------------------------------------------------------------------------
# which datasets contain which location/year
a <- tapply(all.dat$n.ind, list(all.dat$loc.year, all.dat$source), sum)
a[is.na(a)] <- 0
b <- ifelse(a > 0, 1, 0)
head(b)

table(rowSums(b))

# for the same locations, remove the data set with the least individuals,
# use Murray et al. 2010 if equal
# subset to locations with both
sub.a <- a[rowSums(b) == 2, ]
head(sub.a)

# locations with two entries
rem <- data.frame(loc.year = rownames(sub.a)) 

# use Murray et al 2010
# if m >= a, then drop a, else drop m
 rem$drop = ifelse(sub.a[, 2] >= sub.a[, 1], "drop_a", "drop_m")

# else use the Amphibian portal 
# rem$drop = ifelse(sub.a[, 1] >= sub.a[, 2], "drop_m", "drop_a")

table(rem$drop)

all.dat <- full_join(all.dat, rem) |>
  mutate(drop = ifelse(is.na(drop), "keep", drop)) |>
  filter(!(source == "Murray et al 2010" & drop == "drop_m")) |> 
  filter(!(source == "Amphibian disease portal" & drop == "drop_a")) 

glimpse(all.dat)
all.dat
table(all.dat$drop)

table(is.na(all.dat$month))

# number of individuals
sum(all.dat$n.ind)

# number of species and locations
length(levels(factor(all.dat$species)))
length(levels(factor(all.dat$loc.year)))
length(levels(factor(all.dat$loc)))
sum(all.dat$n.ind)
sum(all.dat$n.pos)

#-------------------------------------------------------------------------------
# read in Australian Society of Herpetology (ASH) names
name <- read.csv("./data/ASH frog names.csv",
                 strip.white = T, na.strings = c(" ", "", "NA")) 

# subset to species only (remove subspecies)
name$species <- paste(unlist(lapply(strsplit(name$species_ASH, split = " "), function(x) x[1])),
                      unlist(lapply(strsplit(name$species_ASH, split = " "), function(x) x[2])))

name <- name |>
  dplyr::select(species, family) |>
  unique()

glimpse(name)

all.dat <- left_join(all.dat, name) |>
  mutate(genus = unlist(lapply(strsplit(species, split = " "), function(x) x[1])))

# check everything has matched and has a family
table(all.dat$family, exclude = NULL)

# check for duplicate species x location x year
table(duplicated(cbind(all.dat$species, all.dat$loc, all.dat$year)))

# combine these
all.dat <- all.dat |>
  group_by(species, loc, year, lat, long, family, genus) |>
  summarise(n.ind = sum(n.ind),
            n.pos = sum(n.pos))

#-------------------------------------------------------------------------------
# check for duplicates
table(duplicated(cbind(all.dat$species, all.dat$loc, all.dat$year)))

# total number of individuals per species and number of locations
all.dat <- all.dat |>
  group_by(species) |>
  mutate(total.ind = sum(n.ind),
         total.loc = n()) |>
  glimpse()
  
################################################################################
# filter to include species with >= xx individuals sampled at yy locations
all.dat <- filter(all.dat, total.ind >= 50 & total.loc >= 10)

# exclude cane toads
all.dat <- filter(all.dat, species != "Rhinella marina")
all.dat

################################################################################
# first year infection recorded
first.year <- min(all.dat$year[all.dat$n.pos > 0])
first.year

# exclude observations prior to first.year
all.dat <- all.dat |>
  filter(year >= first.year)

# year range
range(all.dat$year)

################################################################################
# read in list of chytrid impacted frogs from Scheele 2017
imp <- read.csv("./data/Scheele et al 2017 chytrid impacted frogs.csv") |>
  mutate(species = paste(genus, spp),
         decline = "yes") |>
  dplyr::select(species, decline) |>
  unique() |>
  arrange(species)

# merge with infection data
all.dat <- left_join(all.dat, imp)
all.dat$decline[is.na(all.dat$decline)] <- "no"
glimpse(all.dat)

################################################################################
# read in gridded temperature data created using "Download silo temperature data.R"
# read in gridded mean annual temperature data
# this file is not included in the data (too large) but needs to generated to run this
all.temp <- rast("./data/gridded_mean_annual_temp_1889_2021.tif")
all.temp
names(all.temp)

# subset to study years
study.years <- paste0("yr_", min(all.dat$year):max(all.dat$year))
mat <- all.temp[[which(names(all.temp) %in% study.years)]]

# overall mean annual temperature for the study period
mean.mat <- mean(mat)

# variation in mean annual temperature at a location
diff <- max(mat) - min(mat)
diff.val <- terra::extract(diff, cbind(all.dat$long, all.dat$lat))[, 1]
range(diff.val, na.rm = T)

#-------------------------------------------------------------------------------
# overall mean temperature data for chytrid locations
all.dat$temp.overall <- terra::extract(mean.mat, cbind(all.dat$long, all.dat$lat))[, 1]
glimpse(all.dat)

# temperature data specific to each year
# number years
n.year <- all.dat$year - 1977

all.dat$temp.year <- NA
for(i in 1:nrow(all.dat)) {
  all.dat$temp.year[i] <- terra::extract(mat[[n.year[i]]], cbind(all.dat$long[i], all.dat$lat[i]))[, 1]
}
glimpse(all.dat)

# plot spatial versus inter-annual temperature variation
par(mfrow = c(1, 1))
plot(all.dat$temp.year ~ all.dat$temp.overall)

#-------------------------------------------------------------------------------
# output file
write.csv(all.dat, "./data/Tidy Bd prevalence data.csv", row.names = F)

# output species list
list.spp <- data.frame(species = levels(factor(all.dat$species)))
write.csv(list.spp, "./data/Species list.csv", row.names = F)

