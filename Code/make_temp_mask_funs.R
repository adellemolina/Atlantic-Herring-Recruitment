library(terra)

# Load data 
in.file = here::here("test.nc")
# rasterize
temp.rast = terra::rast(in.file)

#crs(temp.rast, proj=T) # "+proj=longlat +datum=WGS84 +no_defs"
#cat(crs(temp.rast), "\n") #"WGS 84 (CRS84)"
#ext(temp.rast)
#ext(temp.rast) <- c(-76.9999978360966, -65.6666664403039, 35.8326971184431, 44.6666720853874)

#set threshold values
min.val = 15
max.val = 25

# load in shapefile
shape.file = here::here( "Data/shapefiles/EPU_extended.shp")
shape.vec <- terra::vect(shape.file)

#check it
plot(shape.vec)
crs(shape.vec) #"+proj=longlat +datum=NAD83 +no_defs"
cat(crs(shape.vec), "\n") #"EPSG",6269

# reproject to match the oisst data
file.shp = terra::project(shape.vec,"+proj=longlat +datum=WGS84 +no_defs+ellps=WGS84")
#file.shp = terra::project(shape.vec, in.file)
crs(file.shp, proj=T) #"+proj=longlat +datum=WGS84 +no_defs"
cat(crs(file.shp), "\n") 
plot(file.shp)
file.shp

#Function to extract value range from GLORYS data (or any gridded data)

## Function to create a mask of a gridded file based on a min and max value
# file.in: full file name for desired netcdf file
# file.shp: shape file read in by terra as a spatvector. Try terra::vect("file.name.shp")
# min.val/max.val: min and max value from gridded variable you want to crop
# Returns a Terra spatRaster of just the cells that are in that data range. 
make_temp_mask = function(file.in, file.shp, min.val, max.val){
  
  data.rast = terra::rotate(terra::rast(file.in))
  
  data.crop = terra::crop(terra::mask(data.rast,file.shp, touches = T),file.shp)
  
  data.window =  terra::clamp(data.crop, lower = min.val, upper = max.val,values = F)
  
  # plot(subset(data.window,300))
  
  return(data.window)
  
}

# make a second version that doesn't rotate
make_temp_mask2 = function(file.in, file.shp, min.val, max.val){
  
  data.rast = terra::rast(file.in)
  
  data.crop = terra::crop(terra::mask(data.rast,file.shp, touches = T),file.shp)
  
  data.window =  terra::clamp(data.crop, lower = min.val, upper = max.val,values = F)
  
  # plot(subset(data.window,300))
  
  return(data.window)
  
}

# apply it to the test nc file for the whole epu shapefile
masked<- make_temp_mask(in.file, file.shp, min.val, max.val)

## Function to calculate sum of days within temperature range
# data: a Terra SpatRaster object (returned from make_temp_mask())
# out.df: T/F flag. Turn on true if you want an X,Y,sum data frame
# Returns a SpatRaster of the sum of days in the input file that are within the min/max range above
make_temp_ndays = function(data, out.df = T){
  
  data.binary = (data * 0)+1
  
  data.binary.n = sum(data.binary,na.rm=T)
  
  if(out.df == T){
    data.n.df = as.data.frame(data.binary.n,xy =T)
    
    data.n.df$sum[which(data.n.df$sum == 0)]= NA  
  }else{
    return(data.binary.n)
  }
  
}
# run it for the test year
masked.days <- make_temp_ndays(masked, out.df = F)
plot(masked.days) # map showing the number of days in the range in the whole area
# min max means how many days in that pixel are in the range

#Function to calculate the total area for input SpatVector shapefile
# data: Input spatRaster from make_temp_mask
# file.shp: Input SpatVect object for area of interest
# units: "m" or "km" 
make_shp_area = function(data,file.shp,units){
  
  data.rast = terra::rotate(terra::rast(data))
  data.rast = terra::subset(data,1)
  
  data.crop = terra::crop(terra::mask(data.rast,file.shp, touches = T),file.shp)
  
  shp.area = expanse(data.crop,unit = units)$area
  
  return(shp.area)
}

#apply to test year to calculate total area in all 4 EPUs
#fall_strat_area <- make_shp_area(terra::rotate(rast(in.file)), file.shp, "km")

#Function to the area that is satisfying the min/max conditions from make_temp_mask()
# data: make_temp_mask(file.in,file.shp,min.val,max.val)
# area.shp: SpatVector of region of interest
# units: "m" or "km" 
make_shp_temp_area = function(data, area.shp,units){
  
  data.binary = (data *0)+1
  # plot(subset(data.binary,300))
  
  ntime = length(terra::time(data.binary))
  
  data.sum =app(data.binary,sum,na.rm=T)
  data.less = (data.sum== ntime)
  data.all = data.sum * data.less
  # plot(data.all)
  
  values(data.all)[values(data.all) == 0] = NA
  # plot(data.all)
  # plot(area.shp,add =T)
  # values(data.mask)
  
  area.mask = expanse(data.all, unit = units)$area
  
  return(area.mask)
  
}

# apply functions
make_shp_temp_area(terra::rotate(rast(in.file)), file.shp, "km") # not sure if supposed to use on raw data file or on outputs from the other thing
# this gives a huge number (larger than that from make shape area)
make_shp_temp_area(masked, file.shp, "km") # not sure if supposed to use on raw data file or on outputs from the other thing
# this gives 0
plot(mas)
