# Used this script preliminarily to explore ecodata data

# Explore data in ecodata package
# already loaded and edited bottom temp, and calanus stage
str(early.SST)
# New possible data to edit

# Chlorophyll has annual from 1998-2019
head(chl_pp)
str(chl_pp)

range(chl_pp$Time)

Primary <- chl_pp

# cold pool index (maybe, but which value to use)
head(cold_pool)
str(cold_pool)
levels(as.factor(cold_pool$Var))


head(energy_density)
str(energy_density)       

# ok there are only a few years
herring.ED <- energy_density%>%
  filter(Species==c("Atl. Herring"))

# gulf stream index
head(gsi)
range(gsi$Time)
# ranges from 1993-2019 and it's only anomaly for all regions

# marine heatwave
head(heatwave)
# perfectly formatted from 1982-2020
range(heatwave$Time)

#ichthyo nah its just diversity
head(long_term_sst) # this is for all regions
head(nefsc_survey)


# Get raw survey data from Sean Lucey code and edit to get what I want

install.packages("Rtools") # doesn't work
devtools::install_github("slucey/RSurvey/Survdat")

library(data.table); library(rgdal); library(Survdat)
library(dplyr);library(sf);library(tidyr)
#-------------------------------------------------------------------------------
nefsc_survey_rdata <- 'Survdat.RData'
#get spatial polygons for Ecological Production Units (EPUs) that are used to clip SST data.
EPU <- readOGR('gis', layer='EPU_extended')


map.crs <- CRS("+proj=longlat +lat_1=35 +lat_2=45 +lat_0=40 +lon_0=-77 +x_0=0
               +y_0=0 +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")


strata <- ecodata::epu_sf %>% as("Spatial")

write.csv(x = survdat, file = "Lucey_NEFSC_Survdat.csv", col.names = T)

strata <- ecodata::epu_sf %>% as("Spatial")

survdat.EPU <- Survdat::poststrat(survdat, strata, 'EPU')
setnames(survdat.EPU, 'newstrata', 'EPU')

write.csv(x = survdat.EPU, file = "Survdat_EPU.csv", col.names = T)

