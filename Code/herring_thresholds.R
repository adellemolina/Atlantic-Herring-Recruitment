##      Title:            Thermal threshold Indicator
##      Script Purpose:   Develop larval/juvenile herring thermal threshold indicator
##      Author:           Adelle Molina
##      Created:          8/12/24
##      Updated:          11/10/24
##      Notes:            Play around with other ideas section of BT one more time and save BT var plots, 
##                        compare diff dates of mean fall SST's...I ran it but didn't "look at it"
##      To Do:            Populate new script w final indicator creation 
#load("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/Workspaces/Busythresh.RData")
#save(OISST.EPU.Daily, strata.map, stress.ts,lethal.ts, stress.ts2, lethal.ts2, opt.start.ts, opt.dur.ts, optBT.prop.ts, optBT.days.ts, MeanBT.ts, mean.SST.ts, file = "Thresholds.Figs.RData")
# Libraries ---------------------------------------------------------------
#libraries
library(raster)
library(sf)
library(terra)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(ecodata)
#remotes::install_github("NEFSC/NEFSC-Spatial")
#library(NEFSC-Spatial)

# Play with daily mean OISST timeseries ----------------------------------------------------
# Previously manually downloaded daily means and calculated average daily temp as average of all pixels in each EPU 
# plot daily temperature 
OISST.EPU.Daily <- OISST%>%
  #dplyr::filter(EPU==c("GOM"))%>%
  ggplot(aes(x=day, y = Value))+
  geom_point(aes(col = year),na.rm = T) + 
  geom_smooth(color="black")+
  facet_wrap(~EPU)+
  scale_x_continuous(breaks = seq(0,365,30),minor_breaks = seq(0,365,5)) +
  ylab("OISST")+
  theme_bw()+
  geom_hline(yintercept=c(0, 16), col="red")+
  geom_hline(yintercept=c(0, 21), col="purple")+
  geom_vline(xintercept=c(90, 181, 273))

# Calculate daily means and deviation from that mean
OISST_anoms <- OISST %>%
  dplyr::group_by(EPU, day) %>%
  dplyr::mutate(daily_mean = mean(Value),
                anom = Value-daily_mean)
# plot
OISST_anoms%>%
  ggplot(aes(x=day, y = daily_mean))+
  geom_line()+
  #geom_point(aes(col = year),na.rm = T) + # change color ramp
  facet_wrap(~EPU)+
  scale_x_continuous(breaks = seq(0,365,30),minor_breaks = seq(0,365,5)) +
  ylab("OISST")+
  theme_bw()+
  geom_hline(yintercept=c(0, 16,20))+
  geom_vline(xintercept=c(90, 181, 273))

# how does this compare to the anoms from ecodata
str(ecodata::seasonal_oisst_anom)
ecodata::seasonal_oisst_anom%>%
  #dplyr::filter(EPU==c("GOM"))%>%
  ggplot(aes(x=Time, y = Value))+
  geom_point(aes(col=EPU), na.rm=T)+
  geom_line(aes(col = EPU),na.rm = T) +
  #geom_smooth(color="black")+
  facet_wrap(~Var)+
  ylab("Seasonal OISST Anomaly")+
  theme_bw()
# Maybe include one or a few of these in the tree
# can be a specific one like fall or winter in GoM

# they also have a gridded one
# is better to look at anomalies or thresholds in values for indicators
ecodata::plot_seasonal_sst_anomaly_gridded(seasonal_sst_anomaly_gridded)

# quick calculations (for non-gridded data)
# number of days where the mean daily temperature for all cells within each EPU is above the threshold
T_thresh <- OISST%>%
  #dplyr::filter(EPU==c("GB"))%>%
  dplyr::group_by(year, EPU)%>%
  dplyr::summarize(number.opt = sum(Value>=16),
                   number.limit = sum(Value>=20),
                   prop.opt = number.opt/365)

ggplot(T_thresh, aes(x=year, col=EPU))+
  geom_line(aes(y=prop.opt))+
  geom_point(aes(y=prop.opt))+
  theme_bw()

# Set up shapefiles (EPU)  ----------------------------------------------------

# Load shapefile (all EPUs)
EPU.areas <- st_read(here::here("Data/shapefiles/EPU_extended.shp"), quiet = TRUE)

# reproject the shapefile to match the raster
EPU.areas.repro <- sf::st_transform(EPU.areas, crs = "+proj=longlat +datum=WGS84 +no_defs  ")

# Plot it
ggplot(EPU.areas.repro) +
  geom_sf(fill = "#69b3a2", color = "white") +
  theme_void()

# make terra format
EPU.vec<- vect(EPU.areas.repro)
plot(EPU.vec) # Check it

# create shapefile for just GoM
GOM_strat <- ecodata::epu_sf%>%filter(EPU == "GOM") 
GOM_strat <- sf::st_transform(GOM_strat, crs = "+proj=longlat +datum=WGS84 +no_defs  ")
GOM_vec<- terra::vect(GOM_strat)
plot(GOM_vec)
ggplot(GOM_strat) +
  geom_sf(fill = "#69b3a2", color = "white") +
  theme_void()

# Set up shapefiles (bottom trawl strata) -----------------------------------------------------
bt_strata <- st_read(here::here('Data/shapefiles/NES_BOTTOM_TRAWL_STRATA.shp'),quiet = TRUE)
herring_spring <- c(01010, 01020, 01030, 01040, 01050, 01060, 01070, 01080, 01090, 
                    01100, 01110, 01120, 01130, 01140, 01150, 01160, 01170, 01180, 
                    01190, 01200, 01210, 01220, 01230, 01240, 01250, 01260, 01270, 
                    01280, 01290, 01300, 01360, 01370, 01380, 01390, 01400, 01610, 
                    01620, 01630, 01640, 01650, 01660, 01670, 01680, 01690, 01700, 
                    01710, 01720, 01730, 01740, 01750, 01760)
herring_fall <- c(01050, 01060, 01070, 01080, 01090, 01100, 01110, 01120, 01130, 
                  01140, 01150, 01160, 01170, 01180, 01190, 01200, 01210, 01220, 
                  01230, 01240, 01250, 01260, 01270, 01280, 01290, 01300, 01360, 
                  01370, 01380, 01390, 01400)

her_strata_fall <- bt_strata %>% 
  sf::st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs  ") %>% 
  dplyr::select(STRATUMA, geometry) %>%
  filter(STRATUMA %in% c('01050', '01060', '01070', '01080', '01090', '01100', '01110', '01120', '01130', 
                         '01140', '01150', '01160', '01170', '01180', '01190', '01200', '01210', '01220', 
                         '01230', '01240', '01250', '01260', '01270', '01280', '01290', '01300', '01360', 
                         '01370', '01380', '01390', '01400'))%>%
  dplyr::summarise(geometry = sf::st_union(geometry))

# make vector for the loops
strata.vec <- terra::vect(her_strata_fall)

# object with total area of the shapefile (denominator in the proportions)
fall_strat_area <- expanse(terra::vect(her_strata_fall),unit = "km")

# Plot strata
plot(her_strata_fall)
plot(strata.vec)
projection(her_strata_fall)
str(her_strata_fall)

ggplot(her_strata_fall)+
  geom_sf()

# load in states map
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))%>%
  st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs  ")
plot(usa)
projection(usa)

xlim<- c(75,65)
ylim<- c(36,46)

# from sarah for bathymetry
bathy <- marmap::getNOAA.bathy(-81,-58, 27, 46)
st_crs(bathy)

# plot it all together
strata.map <- ggplot(her_strata_fall)+
  geom_sf(color = "#2b2b2b", fill = '#cbdbcd', linewidth=1)+
  labs(x = '', y = '') +
  theme_bw() +
  geom_sf(data = usa,color = 'gray20', fill = "lightyellow", size=0.125)+
  coord_sf(xlim = c(-75,-65), ylim = c(38,46), crs =  st_crs(4326), datum= st_crs(4326), expand = FALSE)+
  geom_contour(data = bathy, aes(x=x,y=y,z=-1*z), breaks=c(50,100,150,200, Inf),size=c(0.3),col = 'darkgrey')

# Play with data and functions (test one year) -----------------------------------------------

################# load netcdf to play with one year
test.brick = brick(here::here("test.nc")) # SST
test.brick_BT = brick(here::here("test3.nc")) # BT

# Check it out
test.brick
crs(test.brick)
plot(test.brick)
extent(test.brick)
extent(test.brick_BT)
crs(test.brick_BT)
plot(test.brick_BT)

# Rotate the lat long
test.dat2 <- raster::rotate(test.brick) # raster way
plot(test.dat2[[1]])
test.dat <- terra::rotate(terra::rast(test.brick)) # terra way
plot(test.dat[[1]])

# now crop and plot
test <-  terra::crop(test.dat, EPU.vec, mask= T) # terra way
plot(test[[1]])

test2<- raster::crop(test.dat2, extent(EPU.areas.repro)) # raster way
test2 <- raster::mask(x = test2, mask = EPU.areas.repro) 
plot(test2[[1]]) 

# Using joe's functions (I tweaked them a bit)
head(test.brick)
# Create masked object to input in the next functions
test.mask <- make_temp_mask(file.in = test.brick, min.val = 15, max.val = 20, file.shp = strata.vec) 
plot(test.mask)
GoM_mask <- make_temp_mask(file.in = here::here("sst.day.mean.2004.nc"), min.val = 15, max.val = 20, file.shp = strata.vec )
plot(GoM_mask)

# practice chopping off dates for SST syntax
test.shape <- make_temp_mask(file.in = "sst.day.mean.2004.nc", min.val = 10, max.val = 15, file.shp = strata.vec)
plot(test.shape$sst_1)
time(test.shape)

# make new object, or modify existing, to just keep dates of choice
test.shape.fall <- test.shape[[time(test.shape) >= "2004-08-01"]]
time(test.shape.fall) # check if it worked

# Compute number of days
make_temp_ndays(test.mask, out.df = F) # doesn't work with out.df=T
mask.day <- (make_temp_ndays(GoM_mask, out.df = F)) # doesn't work with out.df=T

# compute total area
gom.area <- make_shp_area(terra::rotate(rast(in.file)), GOM_vec, "km") # total area
make_shp_area(terra::rotate(rast(in.file)), EPU.vec, "km")

#make_shp_area(GoM_mask, GOM_vec, "km") # masked area, I think, should prob be todal

make_shp_temp_area(GoM_mask, GOM_vec, "km")
make_shp_temp_area(terra::rotate(rast(in.file)), GOM_vec, "km")

# calculate percent of area outside the range
#obj <-data.frame(year=2, prop= make_shp_temp_area(test.mask, EPU.vec, "km")/make_shp_area(test.mask, EPU.vec, "km"))
# 37 percent of the entire area was outside the "optimal range in the one year of data I am using as a test
# 
# SST Stress Larval Threshold Days  -------------------------------------------------
# make sure shapefile was set up properly above (It's in two forms, the spatvec and the sf dataframe for plotting)
# NOTE that this has all dates so this also includes summer

#set wd to the folder with the data
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/OISST")

# set years and create empty df
years <- 1983:2023
stress_ndays <- data.frame()

# run the loop for all years for the stress threshold
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  GoM_mask <- make_temp_mask(file.in = name, min.val = 16, max.val = 25, file.shp = her_strata_fall ) # would it change if I used 30 (prob should set to whatever the max is)
  output <- terra::expanse(GoM_mask, "km")/fall_strat_area 
  result <- output %>%
    dplyr::summarize(days = sum(area>=0.5))%>%
    mutate(year=j) 
  stress_ndays <- rbind(stress_ndays, result)
  message(paste("done with", j))
}
stress_ndays # view it
#write.csv(stress_ndays, file = "Data/stress_days_fallstrat.csv", row.names = F)

# plot
stress.ts <- ggplot2::ggplot(stress_ndays, aes(x=year, y = days))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("# of days w/ 50% of fall strata >16")+
  xlab("Year")+
  theme_bw()

# Repeat for just fall dates in fall strata
years <- 1983:2023
fall_stress_ndays <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  strata_mask <- make_temp_mask(file.in = name, min.val = 16, max.val = 21, file.shp = strata.vec) 
  strata_mask <- strata_mask[[time(strata_mask) >= paste0("",j, "-09-01")]] #  Start on 1st-sept through end of december
  output <- terra::expanse(strata_mask, "km")/fall_strat_area 
  result <- output %>%
    dplyr::summarize(days = sum(area>=0.5))%>%
    mutate(year=j) 
  fall_stress_ndays <- rbind(fall_stress_ndays, result)
  message(paste("done with", j))
}
fall_stress_ndays
# note that this doesn't include "lethal" days 
write.csv(fall_stress_ndays, file = "Results/Fall_stress_duration.csv", row.names = F)

stress.ts2 <- ggplot2::ggplot(fall_stress_ndays, aes(x=year, y = days))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("# of days")+
  xlab("Year")+
  ggtitle("# of days in spawning season w/ 50% of fall strata >16")+
  theme_bw()
ggsave("S-D_Stress_Duration.png", width = 8, height = 8, units = 'in')
  
# SST Lethal Larval Threshold Days --------------------------------------------------------
# 2 diff ideas
# 1. fall (for sm larvae, 1 year lag) --> explore this here 
# 2. spring/summer (for large larvae and juveniles) but might need to use the spring strata for this --> maybe explore later

# Repeat as above but for lethal limit
lethal_ndays <- data.frame() 
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  lethal_mask <- make_temp_mask(file.in = name, min.val = 21, max.val = 30, file.shp = strata.vec)
  lethal_mask <- lethal_mask[[time(lethal_mask) >= paste0("",j, "-09-01")]] #  Start on 1st-sept through end of december
  output <- terra::expanse(lethal_mask, "km")/fall_strat_area 
  result <- output %>%
    dplyr::summarize(days = sum(area>=0.5))%>%
    mutate(year=j) 
  lethal_ndays <- rbind(lethal_ndays, result)
  message(paste("done with", j))
}
lethal_ndays # This loop formerly used all days of year (which are mostly summer with a little fall and spring)
# the data exported in the code below is the one for all dates
#write.csv(lethal_ndays, file = "Data/lethal_days_fallstrat.csv", row.names = F)

# Repeated for cropped dates "fall"
# remember this represents number of days that are possibly lethal in your first year of life but only in fall strata on days you'd be in those strata
write.csv(lethal_ndays, file = "Results/Fall_lethal_duration.csv", row.names = F)

lethal.ts2 <- ggplot2::ggplot(lethal_ndays, aes(x=year, y = days))+
  geom_point(na.rm=T)+geom_line(na.rm = T) +
  ylab("# of days")+
  xlab("Year")+
  ggtitle("# of days in spawning season w/ 50% of fall strata >21")+
  theme_bw()
ggsave("Results/S-D_Lethal_Duration.png", width = 8, height = 8, units = 'in')

png("Results/Figures/SST.days.indicator.ts.png", width = 8, height = 8, units = 'in', res = 300)
ggpubr::ggarrange(stress.ts2, lethal.ts2, nrow=2, widths=c(.7,1))
dev.off()
# ah ok so I had not saved the plot for less dates, for all it looked very different, here it's mostly a lot of zeros, not super useful
#gridExtra::grid.arrange(stress.ts2, lethal.ts2)
#ggsave("Durations.Fig.11.18.png", width = 8, height = 8, units = 'in')

# SST Optimal Window start and duration -----------------------------------------------------------------

# Try something slightly different
# when does surface temperature become habitable in fall aka when does it start, small larvae, one year lag)

years <- 1983:2023
lethal_fall_days <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")

  # crop to the fall strata and fall dates 
  data.rast = terra::rotate(terra::rast(name))
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec)
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-09-01")]] #  Start on 1st-september through end of december
  # this has mean daily temp (also min and max) for all cells in fall strata from sept-dec
  # the layers are called sst_dayofyear and each one has temp values for each cell in the shapefile
  
  #turn into a data frame
  #data.df <- terra::as.data.frame(data.crop, xy=TRUE, row.names = T) 

  # do calculations 
  # 1st number of days above limit
  # then also figure out start date of suitable temps for larvae: what is the first day (which doy) when lethal area=0
  
  # select only pixels with temp > 21 # keep lethal blob
  data.crop [data.crop  < 21] <- NA # keep pixels with lethal temps, remove optimal habitat
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  hot.area <- terra::expanse(data.crop, "km", usenames=T) # area of lethal blob in each day, 
  
  # Which day is the first with 0 area or first day of suitable fall temperatures
  result <- hot.area %>%
    dplyr::summarize(start.day = first(layer[area==0]), # this is the start of suitable temps in the "fall"  
                     duration=sum(area==0))%>% # this marks the duration of suitable habitat
    mutate(year=j) 
  lethal_fall_days <- rbind(lethal_fall_days, result)
  message(paste("done with", j))
}
lethal_fall_days
lethal_fall_days <- lethal_fall_days%>%
  mutate(blob.end = as.integer(start.day))

opt.start.ts <- ggplot2::ggplot(lethal_fall_days, aes(x=year, y = blob.end))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  scale_x_continuous(breaks = seq(1980,2025,5))+
  ylab(" DOY")+
  xlab("Year")+
  ggtitle("Start day of optimal fall temp window (SST < 21 C)")+
  theme_bw()+
  geom_hline(yintercept=c(265, mean(lethal_fall_days$blob.end)), col="red") + # red line for the ts. mean 
  geom_hline(yintercept=265, col="blue") # blue line for the  fall equinox

ggsave("Fall_Optimal_Day.png", width = 8, height = 8, units = 'in')

opt.dur.ts <- ggplot2::ggplot(lethal_fall_days, aes(x=year, y = duration))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("Number of days with SST < 21 C")+ 
  xlab("Year")+
  theme_bw()+
  ggtitle("Duration of optimal temp from Sept-Dec")
# lethal summer temperatures exist later into the year aka fall is starting later
# worse for those that spawn early, better for those that spawn later 
#ggsave("S-D_Optimal_Duration.png", width = 8, height = 8, units = 'in')
#write.csv(lethal_fall_days, file = "Results/Duration.Optimal.SST.Sept-Dec.csv", row.names = F)
#write.csv(lethal_fall_days, file = "Results/sublethal_SST.csv", row.names = F)

# Explore Ecodata seasonal gridded bottom temps ---------------------------

# take a look at the seasonal gridded bottom temps
str(ecodata::bottom_temp_seasonal_gridded)
plot(ecodata::bottom_temp_seasonal_gridded)
ecodata::plot_bottom_temp_seasonal_gridded(
  shadedRegion = NULL,
  report = "MidAtlantic",
  scale = "celsius") # this is the average bottom temp at each point for the season 
# is this plot the mean across years or just one year...?

# pull out fall
fall.bt <- ecodata::bottom_temp_seasonal_gridded %>%
  dplyr::filter(Var == "fall")

# now what? 
head(fall.bt)
#Can I just get average temperature on the bottom in fall for fall strata

# Play with Hubert BT data (test one year) ------------------------------------------------

# All the netcdfs for hubert's bt product are in the folder
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/hubert_bt_revised")

# use one file to "check it out"& test w the Joe functions (this works but not sure its "right)
testbt <- make_temp_mask(file.in = "bottom_temp_1994.nc", min.val = 7, max.val = 13, file.shp = strata.vec )
str(testbt) # this has a layer for each day in the year  
# layer has NA if it doesn't meet the "requirements"  
plot(testbt[[1]]) #plot map with the temps in areas that do meet the requirement
time(testbt)
testbt <- testbt[[time(testbt) >= "1994-08-01"]] #  crop to just the fall dates
time(testbt) #something wrong here b/c some years include the first of the next year plus need to deal w leap year)
names(testbt) <- as.integer(lubridate::yday(time(testbt))) # rename the layers to be just day of year 
head(testbt) # could handle end of year stuff here using new names

suit.prop <- terra::expanse(testbt, "km")/fall_strat_area # proportion optimal habitat in each layer, but the layer names got all messed up
mean(suit.prop$area)
range(suit.prop$area)
max(suit.prop$area)

# summarize for the year
suit.results <- suit.prop %>%
  dplyr::summarize(mean.suitable.hab = mean(area, na.rm=T), # average proportion of suitable habitat in fall for the year
                   var = var(area, na.rm = T),
                   max = max(area, na.rm = T)) # variance
 
# TESTING VALUES/METHODS:
# if you use 7-13 get a tightish range of small values
# but if you use 5-20 you get basically same proportion in all days (and in all years, so that's too wide)
# the question is, does that change from year to year (mean was 0.075 in 1994)

# to account for leap year need to basically split up the data and then put back together after summing

# so number of days could just be sum of the layers that meet some condition (I've used proportion of area greater than 50)
# or you can get an area for each day and take the mean 
#(but what about the issue of zeros, current runs don't have this issue but early attempts did b/c it included all dates)
# count them, number of days with no suitable habitat

# BT Average Optimal Area Prop & duration  -------------------------------------------
# Repeat same as below but "correctly" --> the only difference is whether I use the rotate version of Joes function
# does it make sense to average across all the days
# note that unlike SST I used 75% of the area instead of 50 b/c at fifty it was always the same, indicating that the bottom temperatures are well within this range
years <- 1984:2020 # only goes to 2020
BT_opthab <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("bottom_temp_",j, ".nc")
  BT_mask <- make_temp_mask2(file.in = name, min.val = 7, max.val = 13, file.shp = strata.vec) # select pixels that meet criteria
  BT_mask <- BT_mask[[time(BT_mask) >= paste0("",j, "-08-01")]] #  crop to just the fall dates
  BT_prop <- terra::expanse(BT_mask, "km")/fall_strat_area # total area of optimal bt/full area aka proportion optimal in each day 
  BT_results <- BT_prop %>%
    dplyr::summarize(mean.suit = mean(area, na.rm=T), # average proportion of suitable habitat in fall
                     var = var(area, na.rm = T), # variance
                     range = max(area, na.rm = T)-min(area, na.rm = T),#range of suitable habitat
                     days = sum(area>=0.75),
                     days.50 = sum(area>=0.50))%>%  # duration of high cover of suitable habitat
    mutate(year=j)  
  BT_opthab <- rbind(BT_opthab, BT_results)
  message(paste("done with", j))
}
head(BT_opthab) 

optBT.prop.ts <- ggplot2::ggplot(BT_opthab, aes(x=year, y = mean.suit))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("suitable area/total strata area")+
  xlab("Year")+
  ggtitle("Proportion of optimal bottom habitat between 7-13 C")+
  theme_bw()
#ggsave("Fall_Optimal.Prop_BT.png", width = 8, height = 8, units = 'in')

optBT.days.ts <- ggplot2::ggplot(BT_opthab, aes(x=year, y = days))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("# of days ")+
  xlab("Year")+
  ggtitle("Duration of optimal BT (>75% of the area between 7-13 C)")+
  theme_bw()

png("Results/Figures/BT.indicators.ts.png", width = 8, height = 8, units = 'in', res = 300)
ggpubr::ggarrange(optBT.prop.ts, optBT.days.ts, nrow=2, widths=c(.7,1))
dev.off()

# BT Average Optimal Area Prop ("old, wrong" version)-----------------------------
# set up the loop
years <- 1985:2020 # only goes to 2020
hatch_area <- data.frame()

# run the loop
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("bottom_temp_",j, ".nc")
  BT_mask <- make_temp_mask(file.in = name, min.val = 9, max.val = 13, file.shp = strata.vec) 
  BT_mask[[time(BT_mask) >= paste0("",j, "-08-01")]] #  crop to just the fall dates (this isn't right I don't think)
  output <- terra::expanse(BT_mask, "km")/fall_strat_area
  result <- output %>%
    dplyr::summarize(suitable.hab = mean(area, na.rm=T),
                     var = range(area, na.rm = T))%>% # add something to account for zeros, or remove them...? or c
    mutate(year=j) 
  hatch_area <- rbind(hatch_area, result)
  message(paste("done with", j))
}
# So this is the average suitable area proportion for "fall" aka august-december in each year for temps 9-13
hatch_area # but it's very small (is that b/c it's zero inflated? no b/c the rotation is unnecessary) 

ggplot2::ggplot(hatch_area, aes(x=year, y = suitable.hab))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("proportion of suitable bottom temperatures in fall")+
  xlab("Year")+
  theme_bw()

# OTHER IDEAS (play with one more jawn here, ok I did I tried something else and nope) -------------------------------------------------------------
# 1. Count number of days above a threshold area (as I did for stress and lethal SST) --> I did not do this, it might be the same for all days
# 2. Percentage of time in the optimal window (number of days) 
#       dodson used perentage of time (weeks) during the spawning period that were within it (what was the "threshold")
#       did this above by calculating the number of days during the spawning period that have 75% of the area within the optimal range
# 3. Suitable habitat area (mean of the area that is suitable for the days in that season that are optimal (as above)
#     
# 4. use same approach as above for lethal SST when is the area of unsuitable habitat=0? --> aka when does optimal temp start and how long does it persist
#   but essentially temps are "optimal" the whole time, the area of lethal or unsuitable temperatures seems to never be = 0
#   using 5-20 and same with 7-13 I got NA --> this means there is never a time when there are no optimal temps
# So will need to reconceptualize this b/c same approach isn't working
# but I did compute the average area of suboptimal temps (sort of same same but opposite as the proportion approach above)
# should I instead use just a threshold hot value....how many days are too hot instead?
getwd()
years <- 1984:2020
hatching_fall_days <- data.frame()

# run the loop
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("bottom_temp_",j, ".nc")
  
  # crop to the fall strata and fall dates 
  data.rast = terra::rast(name)
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec)
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-08-01")]] #  Start on 1st-august 
 
  # select only pixels with temp outside of 7-13 or 5-20, set all optimal habitat areas to NA 
  data.crop [data.crop  > 7 & data.crop < 13] <- NA # this might be too restrictive maybe need to use some sort of percent of area thing
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  lethal.area <- terra::expanse(data.crop, "km", usenames=T) # area of suboptimal blob in each day, if 0 that means its "all suitable/optimal"

  # Which day is the first with 0 area
  result <- lethal.area %>%
    dplyr::summarize(start.day = first(layer[area==0]),   
                     duration=sum(area==0),
                     mean.area = mean(area))%>% # add a calculation for the average lethal area
    mutate(year=j) 
  hatching_fall_days <- rbind(hatching_fall_days, result)
  message(paste("done with", j))
}

ggplot2::ggplot(hatching_fall_days, aes(x=year, y = duration))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("Days")+
  xlab("Year")+
  ggtitle("Duration of suitable bottom temp")+
  theme_bw()
str(hatching_fall_days)

# this shows how much habitat is inaccessible in each year, the size of the region outside of their optimal range
ggplot2::ggplot(hatching_fall_days, aes(x=year, y = mean.area))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("Area (km")+
  xlab("Year")+
  ggtitle("Area  of sub-optimal bottom temp")+
  theme_bw()

# Repeat similar process using same approach as for SST stressful that I am using above... a max cutoff and look at warm blob isntead of like a window
# Repeat as above but for lethal limit
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/hubert_bt_revised")
years <- 1984:2020
hatching_fall_days2 <- data.frame() 
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("bottom_temp_",j, ".nc")
  
  # crop to the fall strata and fall dates 
  data.rast = terra::rast(name)
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec)
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-08-01")]] #  Start on 1st-august 
  
  # select only pixels with temp less than 20, set all areas that are too hot to NA 
  data.crop [data.crop  > 20] <- NA 
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  cool.area <- terra::expanse(data.crop, "km") /fall_strat_area # calculate the proportional area of nonlethal temperature

  result <- cool.area %>%
    dplyr::summarize(days = sum(area>=0.75))%>%
    mutate(year=j) 
  hatching_fall_days2 <- rbind(hatching_fall_days2, result)
  message(paste("done with", j))
}
ggplot2::ggplot(hatching_fall_days2, aes(x=year, y = days))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("Days")+
  xlab("Year")+
  ggtitle("Number of days with 50% of BT are not too hot")+
  theme_bw()
# Calculate average BT in fall dates/areas --------------------------------
# for hubert bt
getwd()
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/hubert_bt_revised")
years <- 1984:2020 
mean.her.fall.BT <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("bottom_temp_",j, ".nc")
  data.rast = terra::rast(name)
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec) # crop to fall strata
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-08-01")]] #  Crop to Aug-Dec 
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  result <- terra::global(data.crop, "mean", na.rm=T) # average temperature across all pixels for each day
  result <- result %>%
    dplyr::mutate(DOY = as.numeric(row.names(result)))
  an_results <- result %>%
    dplyr::summarize(avg.fall.BT = mean(mean, na.rm=T), # average bottom temp in fall herr strata
                     var = var(mean, na.rm = T), # variance
                     max = max(mean, na.rm = T),
                     min = min(mean, na.rm = T))%>% 
    mutate(year=j)  
  mean.her.fall.BT <-  rbind(mean.her.fall.BT, an_results)
  message(paste("done with", j))
}

MeanBT.ts <-ggplot2::ggplot(mean.her.fall.BT, aes(x=year, y = avg.fall.BT))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("T (C)")+
  xlab("Year")+
  ggtitle("Average BT (Aug-Dec, fall strata)")+
  theme_bw()+ 
  geom_errorbar(aes(ymin=avg.fall.BT-var, ymax=avg.fall.BT+var), width=.2, 
                position=position_dodge(0.05))
ggsave("MeanBT_spawnseason.png", width = 8, height = 8, units = 'in')

# look at the range, min, max, add a plot to do all at once faceted, scales fixl jawn
# might want to keep or use one of these
ggplot2::ggplot(mean.her.fall.BT, aes(x=year, y = var))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("dT (C)")+
  xlab("Year")
# very noticeably a major shift at around 2010. The variance suddenly almost doubles and there seems to be a big shift to higher minima and maxima at that point
# 2018 was the lwest of the last decade, one year before the only small recruitment peak in the new model

# play with doing spatial weighting...not necessary
# first get a spatially weighted average temp for all those days 
#temps*area/total area
# temps=data.crop
area <- expanse(data.crop,unit = "km")
#total area = fall_strat_area
data.sum <- data.crop*area/fall_strat_area
#area$area[1]/fall_strat_area # this is always the same so can prob just take the mean

# Calculate average SST in fall dates/areas  (aug - dec) --------------
# Did this for aug-dec but perhaps it should be september to december
# or maybe even september-november since that seems to more closely match the survey timing
# check WD
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/OISST")
years <- 1983:2023
mean.her.fall.SST <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  data.rast = terra::rotate(terra::rast(name))
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec) # crop to fall strata
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-08-01")]] #  Crop to Aug-Dec 
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  result <- terra::global(data.crop, "mean", na.rm=T) # average temperature across all pixels for each day
  result <- result %>%
    dplyr::mutate(DOY = as.numeric(row.names(result)))
  an_results <- result %>%
    dplyr::summarize(avg.fall.SST = mean(mean, na.rm=T), # average sea temp in fall herr strata
                     range = sd(mean, na.rm = T), # variance
                     max = max(mean, na.rm = T),
                     min = min(mean, na.rm = T))%>% 
    mutate(year=j)  
  mean.her.fall.SST <-  rbind(mean.her.fall.SST, an_results)
  message(paste("done with", j))
}

# Repeat for different fall dates --> lol I ran this but didn't change the dates...idiot
years <- 1983:2023
mean.her.09.12.SST <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  data.rast = terra::rotate(terra::rast(name))
  data.crop = terra::crop(terra::mask(data.rast,strata.vec, touches = T),strata.vec) # crop to fall strata
  data.crop <- data.crop[[time(data.crop) >= paste0("",j, "-09-01")]] #  Crop to Sept-Dec 
  names(data.crop) <- lubridate::yday(time(data.crop)) # rename the layers to be just day of year
  result <- terra::global(data.crop, "mean", na.rm=T) # average temperature across all pixels for each day
  result <- result %>%
    dplyr::mutate(DOY = as.numeric(row.names(result)))
  an_results <- result %>%
    dplyr::summarize(avg.fall.SST = mean(mean, na.rm=T), # average sea temp in fall herr strata
                     range = sd(mean, na.rm = T), # variance
                     max = max(mean, na.rm = T),
                     min = min(mean, na.rm = T))%>% 
    mutate(year=j)  
  mean.her.09.12.SST <-  rbind(mean.her.09.12.SST, an_results)
  message(paste("done with", j))
}

# compare results
mean.SST.ts <- ggplot2::ggplot(mean.her.fall.SST, aes(x=year, y = avg.fall.SST))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("T (C)")+
  xlab("Year")+
  ggtitle("Average SST (Aug-Dec, fall strata)")+
  theme_bw()

mean.SST.ts2 <- ggplot2::ggplot(mean.her.09.12.SST, aes(x=year, y = avg.fall.SST))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("T (C)")+
  xlab("Year")+
  ggtitle("Average SST (Sept-Dec, fall strata)")+
  theme_bw()
# same shape ust scaled a bit diff
  # +geom_errorbar(aes(ymin=avg.fall.SST-var, ymax=avg.fall.SST+var), width=.2, position=position_dodge(0.05))
# Combine indicator time series -------------------------------------------
# reformat, and combine mean BT with SST as well as stress duration jawn maybe is the best and add lags for inclusion in BRT
# not stress duration...we want the optimal window duration
# uh suddenly I'm thinking the mean SSt should not go in b/c like it is weird...
#OISST  <- read.csv(here::here('Data/Inputs/OISST.csv'), header = T, blank.lines.skip = T)

thermal.indicators<- merge(mean.her.fall.BT, mean.her.fall.SST, by = "year",  all.x = TRUE, all.y = T)
thermal.indicators <- thermal.indicators%>%
  merge(Duration_Optimal_SST_Sept_Dec, by = "year",  all.x = TRUE, all.y = T)%>%
  select(year, duration, avg.fall.BT, avg.fall.SST)%>%
  rename(days.optimal.SST=duration)

# Play with spatial plotting ---------------------------------------------------

ggplot() + 
  geom_sf(data = test2,color = 'gray20', fill = '#cbdbcd') +
  geom_contour(data = EPU.areas.repro,
               aes(x=x,y=y,z=-1*z),
               breaks=c(50,100,150,200, Inf),
               size=c(0.3),
               col = 'darkgrey') +
  stat_summary_2d(aes(),
                  binwidth=c(0.16666,0.16666)) + 
  scale_fill_gradientn(colors = jet.colors(20)) +
  coord_sf(xlim = c(-75,-65.5), ylim = c(36,44), datum = sf::st_crs(4326))  +
  labs(x = '', y = '', fill = 'Bottom temperature (°C)') +
  theme_bw() 

ggplot() +  
  geom_tile(data=test2, aes(x=x, y=y, fill=value), alpha=0.8) + 
  geom_polygon(data=EPU.areas.repro, aes(x=long, y=lat, group=group), 
               fill=NA, color="grey50", size=0.25) +
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"))

test_df <-
  as.data.frame(test2, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() %>%
  #--- change the variable names ---#
  data.table::setnames(
    paste0("tmax_Jan_09.", 1:5),
    seq(ymd("2009-01-01"), ymd("2009-01-05"), by = "days") %>%
      as.character()
  )

ggplot(test_df)+
  geom_raster(aes(x = x, y = y, fill = `X1981.09.01`)) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "bottom")

test_df_long <- 
  test_df %>%
  tidyr::pivot_longer(
    c(-x, -y),
    names_to = "date",
    values_to = "dailyt")

ggplot(test_df_long)+
  geom_raster(aes(x = x, y = y, fill = dailyt)) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "bottom")+
  facet_wrap(date ~ .) +
  coord_equal()

