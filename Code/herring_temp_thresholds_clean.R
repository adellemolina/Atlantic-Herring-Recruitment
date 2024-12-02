##      Title:            Thermal threshold Indicator Clean
##      Script Purpose:   Run loops to create larval/juvenile herring thermal threshold indicators
##      Author:           Adelle Molina
##      Created:          11/10/24
##      Updated:          11/11/24

ggplot(WHAM, aes(x=year, y=logR))+# plot obs vs pred
  theme_bw()+  
  geom_point(na.rm=T, col="blue", size=2)+
  geom_point(data=fitted, mapping=aes(x=year, y=brt.pred)) +
  geom_line(data=fitted, mapping=aes(x=year, y=brt.pred)) +
  #geom_line(data=fitted.test, mapping=aes(x=year, y=brt.pred), col="blue") +
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +  ylab("log Recruitment")
# Set up shapefiles -------------------------------------------------------
bt_strata <- st_read(here::here('Data/shapefiles/NES_BOTTOM_TRAWL_STRATA.shp'),quiet = TRUE)
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

# make an object with total area of the shapefile (denominator in the proportions)
fall_strat_area <- expanse(terra::vect(her_strata_fall),unit = "km")

# OISST Indicators --------------------------------------------------------

##      STRESS
#set wd to the folder with the data
setwd("C:/Users/adelle.molina/Desktop/Atlantic_Herring/Thresholds/OISST")

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

# OPTIMAL/Not lethal
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
lethal_fall_days <- lethal_fall_days%>%
  mutate(blob.end = as.integer(start.day))


# Bottom Temp Indicators --------------------------------------------------

# set wd
setwd("C:/Users/adelle.molina/Documents/GitHub/Atlantic-Herring-Recruitment/data-raw/gridded/hubert_bt_revised")

