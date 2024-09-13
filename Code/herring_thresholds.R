##      Title:            Thermal threshold Indicator
##      Script Purpose:   Develop larval/juvenile herring thermal threshold indicator
##      Author:           Adelle Molina
##      Created:          8/12/24
##      Updated:          8/19/24
##      Notes:            
##      To Do:            Run the loop for all years, woop now some fxs not working

# Libraries ---------------------------------------------------------------
#libraries
library(raster)
library(sf)
library(terra)
library(ggplot2)
library(ncdf4)
library(dplyr)
library(ecodata)
remotes::install_github("NEFSC/NEFSC-Spatial")
library(NEFSC)

# Play with daily mean OISST timeseries ----------------------------------------------------
# Previously manually downloaded daily means and calculated average daily temp as average of all pixels in each EPU 
# plot daily temperature 
OISST%>%
  #dplyr::filter(EPU==c("GOM"))%>%
  ggplot(aes(x=day, y = Value))+
  geom_point(aes(col = year),na.rm = T) + 
  geom_smooth(color="black")+
  facet_wrap(~EPU)+
  scale_x_continuous(breaks = seq(0,365,30),minor_breaks = seq(0,365,5)) +
  ylab("OISST")+
  theme_bw()+
  geom_hline(yintercept=c(0, 16), col="red")+
  geom_vline(xintercept=c(90, 181, 273))

# Calculate daily means and deviation from that mean
OISST_anoms <- OISST %>%
  dplyr::group_by(EPU, day) %>%
  dplyr::mutate(daily_mean = mean(Value),
                anom = Value-daily_mean)

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

# plot up 
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
# hmmm maybe instead of my crudely calculated regional oisst's I use these in the tree
# can be a specific one like fall or winter in GoM
# they also have a gridded one so maybe its better to look at anomalies than thresholds in values
ecodata::plot_seasonal_sst_anomaly_gridded(seasonal_sst_anomaly_gridded)

# will need to redo with gridded data but for quickly calculate number of days above thresholds
# so this is the number of days where the mean daily temperature for all cells within each EPU
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
plot(EPU.vec)
write('CURL_SSL_BACKEND=openssl', file = "~/.Renviron", append = TRUE)


# should prob repeat with the herring spring and fall strata instead of using all EPUs, here use just GoM
GOM_strat <- ecodata::epu_sf%>%filter(EPU == "GOM") # just GOM, try to see if this works
GOM_strat <- sf::st_transform(GOM_strat, crs = "+proj=longlat +datum=WGS84 +no_defs  ")
GOM_vec<- terra::vect(GOM_strat)
plot(GOM_vec)
ggplot(GOM_strat) +
  geom_sf(fill = "#69b3a2", color = "white") +
  theme_void()

# Set up shapefiles (from bottom trawl) -----------------------------------------------------

bt_strata <- st_read(here::here('Data/shapefiles/NES_BOTTOM_TRAWL_STRATA.shp'),quiet = TRUE, )

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
  sf::st_transform(proj4string = "+proj=longlat +datum=WGS84 +no_defs  ") %>% 
  dplyr::select(STRATUMA, geometry) %>%
  filter(STRATUMA %in% c('01050', '01060', '01070', '01080', '01090', '01100', '01110', '01120', '01130', 
                         '01140', '01150', '01160', '01170', '01180', '01190', '01200', '01210', '01220', 
                         '01230', '01240', '01250', '01260', '01270', '01280', '01290', '01300', '01360', 
                         '01370', '01380', '01390', '01400'))%>%
  dplyr::summarise(geometry = sf::st_union(geometry))

################# load netcdfs
test.brick = brick(here::here("test.nc")) # start with just one year

# Check it out
test.brick
crs(test.brick)
plot(test.brick)
extent(test.brick)

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

# Calculate Threshold Areas --> test one year first-----------------------------------------------
# Used joe's functions (I tweaked them a bit)

# Create masked object to input in the next functions
test.mask <- make_temp_mask(file.in = test.brick, min.val = 15, max.val = 20, file.shp = EPU.vec ) 
plot(test.mask)
GoM_mask <- make_temp_mask(file.in = here::here("test.nc"), min.val = 15, max.val = 20, file.shp = GOM_vec )
plot(GoM_mask)

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
# Calculate proportion outside optimal range time series -------------------------------------------------

# write a loop to do this for all the files
setwd("C:/Users/adelle.molina/Desktop/Atlantic_Herring/Thresholds/OISST")

# make sure shapefile was set up properly above
years <- 1983:2023
stress_ndays <- data.frame()
for(j in years) {
  message(paste("starting", j)) 
  name <- paste0("sst.day.mean.",j, ".nc")
  GoM_mask <- make_temp_mask(file.in = name, min.val = 16, max.val = 25, file.shp = GOM_vec )
  output <- terra::expanse(GoM_mask, "km")/73801.24 
  result <- output %>%
    dplyr::summarize(days = sum(area>=0.5))%>%
    mutate(year=j) 
  stress_ndays <- rbind(stress_ndays, result)# maybe this should create a list, with one data frame for each year
  # then here you count up the number of rows and print that as the number of days
  # mutate sum_if GoM_prop$area >=0.5
  message(paste("done with", j))
}
stress_ndays

write.csv(stress_ndays, file = "Data/stress_days.csv", row.names = F)

stress.ts <- ggplot2::ggplot(stress_ndays, aes(x=year, y = days))+
  geom_point(na.rm=T)+
  geom_line(na.rm = T) +
  ylab("# of days w/ 50% of GoM >16")+
  xlab("Year")+
  theme_bw()

dev.off()
# Abby code Example--------------------------------------------------

for(j in years) {
  message(paste("starting", j)) # will need to use code like the one I have in the folder for manually downloaded data
  # download data ----
  #dir.create(here::here("data-raw","gridded", "sst_data"), recursive = TRUE)
  #url <- paste0("https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.", j, ".nc")
  #download.file(url, destfile = "j.nc")
  
  name <- paste0("sst.day.mean.",j, ".nc")
  #name <- "test.nc"
  
  data <- ecopull::nc_to_raster(nc = name, varname = 'sst') # converts to NAD83
  data <- raster::rotate(data)
  message("converted to raster...")
  
  # make sure all days are there ----
  
  if(raster::nlayers(data) < 365) {
    message(j, " does not have a full year of data! skipping!")
  } else {
    
    # crop to gom ----
    ndays <- raster::nlayers(data) # account for leap years
    
    gom_temp <- raster::mask(x = data[[1:180]], 
                             mask = EPU.vec # gom
    )
    gom_temp2 <- raster::mask(x = data[[181:ndays]], 
                              mask = EPU.vec # gom
    )
    message("cropped to EPU...")

    # calculate total area ----
    raster_areas <- raster::area(gom_temp, na.rm = TRUE) %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y))
    raster_areas2 <- raster::area(gom_temp2, na.rm = TRUE) %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y))
    
    total_area <- raster_areas %>%
      colSums(na.rm = TRUE) %>%
      unique() # always the same (should do earlier/simpler)
    
    temps_df <- gom_temp %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y))
    temps_df2 <- gom_temp2 %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y))
    
    # calculate weighted mean temp
    weighted_mean_temp <- (temps_df * raster_areas / total_area) %>%
      colSums(na.rm = TRUE)
    weighted_mean_temp2 <- (temps_df2 * raster_areas2 / total_area) %>%
      colSums(na.rm = TRUE)
    
    
    # calculate area 18-25.6C ----
    gom_temp@data@values[which(gom_temp@data@values < 18)] <- NA
    # also remove areas that are too warm
    gom_temp@data@values[which(gom_temp@data@values > 25.6)] <- NA
    
    gom_temp2@data@values[which(gom_temp2@data@values < 18)] <- NA
    # also remove areas that are too warm
    gom_temp2@data@values[which(gom_temp2@data@values > 25.6)] <- NA
    
    warm_area <- raster::area(gom_temp, na.rm = TRUE) %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y)) %>%
      colSums(na.rm = TRUE)
    warm_area2 <- raster::area(gom_temp2, na.rm = TRUE) %>%
      raster::as.data.frame(xy = TRUE) %>%
      dplyr::select(-c(x, y)) %>%
      colSums(na.rm = TRUE)
    
    area_data <- tibble::tibble(names = c(names(weighted_mean_temp), names(weighted_mean_temp2)),
                                warm_prop = c(warm_area, warm_area2)/total_area,
                                mean_temp = c(weighted_mean_temp, weighted_mean_temp2)
    )
    
    result <- area_data %>%
      dplyr::mutate(Year = stringr::str_extract(names, pattern = "\\d{4}"), 
                    names = stringr::str_remove(names, pattern = "X"), 
                    DOY = lubridate::as_date(names))
    
    # first and last days with mean temp > 18
    #this_first <- gom_prop %>%
      #dplyr::filter(mean_temp >= 18) %>%
      #dplyr::arrange(DOY)
    
    #first <- c(first, lubridate::yday(this_first$DOY[1]))
    #last <- c(last, lubridate::yday(this_first$DOY[nrow(this_first)]))
    #message("calculated mean temp...")
    
    # number of days with >=75% of area >18
    #this_n_days <- gom_prop %>%
      #dplyr::filter(warm_prop >= 0.75) %>%
      #dplyr::arrange(DOY)
    
    #n_days <- c(n_days, length(this_n_days$DOY))
    #message("calculated proportion...")
  }
  message(paste("done with", j))
}

# PLotting (from Sarah) ---------------------------------------------------

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


tail(test2)
# From examples...not sure what's going on here -----------------------------------------------------------

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

#--- take a look ---#
head(test_df)

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

# Plotting from sarah salois code
jet.colors <-colorRampPalette(c("blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# select just years with study fleet bottom temps
test.dat <- sfob.env %>% filter(year>2006 & depth > 50) 
yrs = sort(unique(sf.bt$year))    

#for(i in 1:length(yrs)){
tempmap <- function(dat, yrs){
  #dat %>% filter(year == yrs) %>% 
    ggplot() + 
    geom_sf(data = EPU.areas %>% st_as_sf(),color = 'gray20', fill = '#cbdbcd') +
    geom_contour(data = EPU.areas,
                 aes(x=x,y=y,z=-1*z),
                 breaks=c(50,100,150,200, Inf),
                 size=c(0.3),
                 col = 'darkgrey') +
    stat_summary_2d(aes(x=start_lon, y=start_lat, z = bottomt),
                    binwidth=c(0.16666,0.16666)) + 
    scale_fill_gradientn(colors = jet.colors(20)) +
    coord_sf(xlim = c(-75,-65.5), ylim = c(36,44), datum = sf::st_crs(4326))  +
    labs(x = '', y = '', fill = 'Bottom temperature (°C)') +
    theme_bw() 
}


for(i in 1:length(yrs)){
  cat("\n######",  as.character(yrs[i]),"\n")
  print(yrmap(yrs[i])) 
  cat("\n")   
}

names(test2)


