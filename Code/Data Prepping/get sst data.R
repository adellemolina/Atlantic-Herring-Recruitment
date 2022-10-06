
##      Title:            Get, prepare, and plot various environmental variables
##      Script Purpose:   Trim some temp datasets, saved as csv's, which can be loaded in separetely 
##      Author:           Adelle Molina
##      Created:          6/21/22
##      Updated:          9/1/22
##      To Do:            Plot mean annual temp in each of the four regions, calculate mean seasonal temp in GoM
##      Notes:            Go through data folder and keep what i need, explore data I haven't looked at, delete what isn't necessary
##      Data Questions    1. What is what is bcodmo?
##                        2. ichthyo is diversity
##                        3. ilong_term_sst is for all regions
##                        4. Energy Density
##      Comments:         Used this script to prep time series, don't really need anymore
# For herring, only 2017 & 2021
herring.ED <- energy_density%>%
  filter(Species==c("Atl. Herring"))

#libraries
library(ncdf4)
library(dplyr)
library(readr)
library(tidyr)
library(sp)
library(rgdal)
library(raster)
library(stringr)
library(lubridate)
library(plyr)
library(usethis)


# OISST -------------------------------------------------------------------

# set working directory to load spatial polygons
setwd("c:/users/adell/documents/github/tech-doc")

#get spatial polygons for Ecological Production Units (EPUs) that are used to clip SST data.
EPU <- readOGR('gis', layer='EPU_extended')


map.crs <- CRS("+proj=longlat +lat_1=35 +lat_2=45 +lat_0=40 +lon_0=-77 +x_0=0
               +y_0=0 +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#find long term daily mean SSTs, create empty vectors to fill for each of the 4 regions
MAB_sst_daily_mean <- NULL
GB_sst_daily_mean <- NULL
GOM_sst_daily_mean <- NULL
SS_sst_daily_mean <- NULL

# now set new wd that has the sst data in it
setwd("c:/users/adell/documents/sst_data/")

# Ok so put all the data in folders called 1,2, or 3 and you can run this all at once or separately
for (dir. in 3){
  
  #Loop through directories
  setwd(paste0('c:/users/adell/documents/sst_data/',dir.))
  print(getwd())
  
  for (f in 1:length(list.files())){
    
    if (!str_detect(list.files()[f],".nc")){
      print(paste(list.files()[f],"is not a raster")) #Based on file type
      next
    }
    
    for (j in c("MAB","GB","GOM","SS")){
      
      sub_region <- EPU[EPU@data$EPU == j,]
      y <- as.numeric(str_extract(list.files()[f],"[0-9]+")) #get year
      
      for (i in 1:365){
        print(paste(j,y,i))
        daily_mean <- raster(paste0(list.files()[f]), band = i) #get band
        
        #set crs
        daily_mean@crs <- sub_region@proj4string 
        
        
        #rotate to lon scale from 0-360 to -180-180
        daily_mean <- rotate(daily_mean)
        
        #mask raster with spatialpolygon
        daily_mean_clipped <- mask(daily_mean, sub_region)
        
        
        #add mean value to data.frame
        assign(paste0(j,"_sst_daily_mean"),
               rbind(get(paste0(j,"_sst_daily_mean")),
                     c(mean(daily_mean_clipped@data@values, na.rm = T),y,i)))
        
      }
    }
    
    
  }
}

# Create separate object for each region
mab3 <- data.frame(EPU = "MAB",
                  year = MAB_sst_daily_mean[,2],
                  day = MAB_sst_daily_mean[,3],
                  Value = MAB_sst_daily_mean[,1])
gb3 <- data.frame(EPU = "GB",
                 year = GB_sst_daily_mean[,2],
                 day = GB_sst_daily_mean[,3],
                 Value = GB_sst_daily_mean[,1])

gom3 <- data.frame(EPU = "GOM",
                  year = GOM_sst_daily_mean[,2],
                  day = GOM_sst_daily_mean[,3],
                  Value = GOM_sst_daily_mean[,1])

ss3 <- data.frame(EPU = "SS",
                  year = SS_sst_daily_mean[,2],
                  day = SS_sst_daily_mean[,3],
                  Value = SS_sst_daily_mean[,1])

# combine and save each of the three objects after each run
early.SST <- rbind(mab, gb, gom, ss)
save(early.SST, file = "early_sst.Rdata")

mid.SST <- rbind(mab2, gb2, gom2, ss2)
save(mid.SST, file = "mid_sst.Rdata")

late.SST <- rbind(mab3, gb3, gom3, ss3)
save(late.SST, file = "late_sst.Rdata")
str(late.SST)

# save as a csv
write.csv(late.SST,
          file="OISST.csv")

# BOTTOM TEMP from SURVEY ------------------------------------------------------

# Load in bottom temp datasets from ecodata folder 
GB <- read.csv(file.choose(), header = T, blank.lines.skip = T)
GOM <- read.csv(file.choose(), header = T, blank.lines.skip = T)
MAB <- read.csv(file.choose(), header = T, blank.lines.skip = T)
SS <- read.csv(file.choose(), header = T, blank.lines.skip = T)


get_bottom_temp <- function(save_clean = F){
  
  SS <- SS%>% mutate(EPU = "SS")
  GOM <- GOM%>% mutate(EPU = "GOM")
  GB <- GB%>% mutate(EPU = "GB")
  MAB <- MAB%>% mutate(EPU = "MAB")
  
  bottom_temp <- rbind(SS, GOM, GB, MAB) %>% #bind all
    dplyr::mutate(Units = "degreesC", Time = as.Date(format(date_decimal(Time), "%Y-%b-%d"), "%Y-%b-%d"),
                  Var, Var = plyr::mapvalues(Var, from = c("Tsfc_anom",#Rename variables
                                                           "Tsfc_ref",
                                                           "Tbot_anom",
                                                           "Tbot_ref"),
                                             to = c("sst anomaly in situ",
                                                    "reference sst in situ (1981-2010)",
                                                    "bottom temp anomaly in situ",
                                                    "reference bt in situ (1981-2010)"))) %>%
    dplyr::group_by(Time = year(Time), EPU, Var, Units) %>%
    dplyr::summarise(Value = mean(Value)) %>%
    as.data.frame()
  
  
  if (save_clean){
    usethis::use_data(bottom_temp, overwrite = T)
  } else {
    return(bottom_temp)
  }
}
get_bottom_temp(save_clean = T)

# calculate annual mean temp by adding the reference to the anomaly

bottom_temp_edited<-bottom_temp %>%
  select(Time, EPU, Var, Value) %>%
  dplyr::filter(Var%in%c("bottom temp anomaly in situ", "reference bt in situ (1981-2010)"))%>% 
  dplyr::group_by(Time, EPU) %>%
  mutate(ActualBT = Value + Value[Var==c("bottom temp anomaly in situ")])%>%
  dplyr::filter(!Var==c("bottom temp anomaly in situ"))


# NEFSC data ---------------------------------------------------------

# use this to to extract salinity and temp againn
NEFSC <- read.csv(file.choose(), header = T, blan.lines.skip = T)
NEFSC <- NEFSC %>%
  mutate_if(is.integer, as.numeric)%>%
  mutate_if(is.character, as.factor)
head(NEFSC) 

NEFSC.anmean <- NEFSC
# going to have to trim it to the EPU, which I'm not sure how to do

# Other Temperature Data --------------------------------------------------

# sabas bottom temp, 
 


# Salinity Data -----------------------------------------------------------


# These data include in situ regional time series of both surface and bottom salinity anomalies
# on the Northeast Continental Shelf. Raw data is split into four files by EPU (SS, GOM, GB, and MAB).

library(dplyr)
library(tidyr)
library(lubridate)

raw.dir <- here::here("data-raw") #input raw
# Shoot dont have this data

get_oceansal_insitu <- function(save_clean = F){
  ss <- read.csv(file.path(raw.dir,"EcoSS_core_Stopbot.csv")) %>% mutate(EPU = "SS")
  gom <- read.csv(file.path(raw.dir,"EcoGoM_core_Stopbot.csv")) %>% mutate(EPU = "GOM")
  gb <- read.csv(file.path(raw.dir,"EcoGB_core_Stopbot.csv")) %>% mutate(EPU = "GB")
  mab <- read.csv(file.path(raw.dir,"EcoMAB_core_Stopbot.csv")) %>% mutate(EPU = "MAB")
  
  oceansal_insitu <- rbind(ss, gom, gb, mab) %>% #bind all
    dplyr::rename(Time = decimal.year, Var = variable.name, Value = salinity) %>% #rename
    dplyr::mutate(Units = "PSU", Time = lubridate::as.Date(format(date_decimal(Time), "%Y-%b-%d"), "%Y-%b-%d"),
                  Var, Var = plyr::mapvalues(Var, from = c("Ssfc_anom",
                                                           "Ssfc_ref",
                                                           "Sbot_anom",
                                                           "Sbot_ref"),
                                             to = c("surface salinity anomaly in situ",
                                                    "reference surface salinity in situ (1981-2010)",
                                                    "bottom salinity anomaly in situ",
                                                    "reference bottom salinity in situ (1981-2010)"))) %>%
    dplyr::group_by(Time = year(Time), EPU, Var, Units) %>%
    dplyr::summarise(Value = mean(Value)) %>%
    as.data.frame()
  
  if (save_clean){
    usethis::use_data(oceansal_insitu)
  } else {
    return(oceansal_insitu)
  }
  
}


# Copepod  ------------------------------------------------------------

# First manually load in the r data file, and save it as an object
copepod <- calanus_stage
str(copepod)
copepod <- copepod%>%
  mutate_if(is.character, as.factor)%>%
  dplyr::filter(EPU==c("GOM"))%>% 
  dplyr::group_by(Year, season, Var) %>%
  dplyr::summarise(Anmean = mean(Value)) %>%
  as.data.frame()

#save as a csv
