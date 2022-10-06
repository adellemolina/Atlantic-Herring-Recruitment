

st_join(x = NEFSC.ctd, EPU.final[,EPU.final$EPU == "GOM"])
gom.fart = EPU.final[EPU.final$EPU == "GOM",]
epu.SIMP <- EPU.final%>% 
  select(EPU)

# old map stuff that didn't work trying to export spatial polygons n shit wasn't necessry just add geom and back to df
str(NEFSC.ctd)
class(NEFSC.ctd)
is.projected(NEFSC.ctd) # see if it is projected, it isn't
proj4string(NEFSC.ctd) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
writeOGR(NEFSC.ctd, "gis/ctdproj", "ctdproj", driver="ESRI Shapefile", overwrite_layer = T) # save shapefile
surveyctd <- readOGR("gis", "ctdproj") # load it back in

class(surveyctd)
geometry(EPU)
IDKAGAIN <- over(EPU, geometry(surveyctd))
class(IDKAGAIN)
COMBINED <- rbind(surveyctd, EPU)
names(surveyctd)
# cHECK THAT THEY MATCH
proj4string(surveyctd)
proj4string(EPU)
is.projected(surveyctd)
range(coordinates(surveyctd))
range(coordinates(EPU))
range(coordinates(NEFSC.ctd))
class(surveyctd)
geometry(EPU)
IDKAGAIN <- over(EPU, geometry(surveyctd))
class(IDKAGAIN)
COMBINED <- rbind(surveyctd, EPU)
names(surveyctd)