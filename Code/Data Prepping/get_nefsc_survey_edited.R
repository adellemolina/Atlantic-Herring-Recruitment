 #Set up for finding survey proportions

#This script was adapted from code written by Sean Lucey at the NEFSC. Note that you
#will need the latest survdat file to get this to run successfully.

#-------------------------------------------------------------------------
#Required packages
devtools::install_github("slucey/RSurvey/Survdat")

library(data.table); library(rgdal); library(Survdat)
library(dplyr);library(sf);library(tidyr);library(sp)
#-------------------------------------------------------------------------------
# Set directory and load raw data, from ecodata
data.dir <- here::here("data-raw")
load(file.path(data.dir, nefsc_survey_rdata ))
nefsc_survey_rdata <- survdat

# transform and reproject
#get spatial polygons for Ecological Production Units (EPUs) that are used to clip SST data.
EPU <- readOGR('gis', layer='EPU_extended')
# Set projection
crs <- CRS("+proj=longlat +lat_1=35 +lat_2=45 +lat_0=40 +lon_0=-77 +x_0=0
               +y_0=0 +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Reproject
EPU_sf <- spTransform(EPU, crs)

#Post stratify data and rename columns
survdat.EPU <- Survdat::poststrat(survdat, EPU_sf, 'EPU')
setnames(survdat.EPU, 'newstrata', 'EPU')

# Save this version
write.csv(survdat.EPU, file = "Survdat EPU.csv")

# ok diff versions of data as csv, the one called lucey survdat only goes up to 2018 and is not clipped
demo  <- read.csv(file.choose(), header = T, blank.lines.skip = T) # survey temperature anomalies
range(survdat.EPU$LAT)

# old methods to reproject
strata <- ecodata::epu_sf %>% as("Spatial") # this spits an error that it wasn't retransformed
strata <- map.crs %>% as("Spatial")
epu_sf <- as(EPU, "sf")

# Raw data has 2,892,498 rows
# Data I just clipped has 2,625,840 (correct projection and years)
# Data I loaded in (from a csv I previously created or saved, but not sure where or how) has 2,523,356 --> correct projection but wrong years, I must have done the projection and saved it using old version of the data (from seans git), now I can only find the 2019 one 
range(demo$YEAR)
range(nefsc_survey_rdata$LAT)
range(NEFSC.RAW$YEAR)

# lol just found the code I used to save this