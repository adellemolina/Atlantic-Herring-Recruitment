# Simple clean threshold code example reproducible

# set up shapefile
GOM_strat <- ecodata::epu_sf%>%filter(EPU == "GOM") # just GOM, try to see if this works
GOM_strat <- sf::st_transform(GOM_strat, crs = "+proj=longlat +datum=WGS84 +no_defs  ")
GOM_vec<- terra::vect(GOM_strat)

# set up netcdf 
in.file = here::here("test.nc")

GoM_test <- make_temp_mask(file.in = in.file, min.val = 15, max.val = 20, file.shp = GOM_vec )
plot(GoM_test)

make_shp_temp_area(GoM_test, GOM_vec, "km")

new.cut <- terra::expanse(GoM_test, "km") # this gives the area for each day that is in the range
new.cut$prop <- new.cut$area/53083.24
# now add an iff then

ggplot() +
  tidyterra::geom_spatraster(data = GoM_test) +
  facet_wrap(~lyr, ncol = 20) # plots all days, edit to plot only when not NA

ggplot() +
  tidyterra::geom_spatraster(data = GoM_test, aes(fill = sst_29))

                  