##      Title:            Get, prepare, combine, and plot environmental & biological variables
##      Script Purpose:   Prepare multivariate time series for analysis
##      Author:           Adelle Molina
##      Created:          7/22/22
##      Updated:          8/19/24
##      Notes:            Change temp source data (to what)
##      To Do:            
# Libraries ---------------------------------------------------------------
#remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
library(ecodata)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(timetk)
library(corrplot)
library(reshape2)
library(ggstatsplot)
library(data.table)
library(gbm)
library(dismo)
library(gridExtra)

# 1. Environmental data ---------------------------------------------------

# Bring in the Regional Areas (compute area weight)
EPUwt <- ecodata::epu_sf%>%
  as.data.frame()%>%
  dplyr::mutate(Weight=Shape_Area/sum(Shape_Area))%>%
  dplyr::select(EPU, Weight)

#     A. OISST -------------------------------------------------------------------
# Manually downloaded daily mean SST (1982-present) from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# In "get sst data" script" computed daily mean by EPU, combined, and exported, load that in
# perhaps something like daily max or daily var is better...
OISST  <- read.csv(here::here('Data/OISST.csv'), header = T, blank.lines.skip = T)

# Calculate basin wide annual mean (no area weighting)
OISST.Annual  <- OISST %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(OISST = mean(Value))%>%
  dplyr::select(year,OISST)%>%
  unique()

# Calculate annual means by region 
OISST.Reg <- OISST %>%
  dplyr::select(year, EPU, day, Value) %>%
  dplyr::group_by(year, EPU) %>%
  dplyr::summarise(Annual.Mean = mean(Value))

#combine with basin wide mean into one object
OISST.All <- OISST.Reg %>%
  pivot_wider(names_from = EPU, values_from = Annual.Mean) %>%
  dplyr::rename(OISST.GB=GB,
                OISST.GOM=GOM,
                OISST.MAB=MAB,
                OISST.SS=SS)%>%
  merge(OISST.Annual, by = "year", all.x = TRUE, all.y = T)

#     B. BT  -------------------------------------------------------
#str(bottom_temp) #in situ anomalies from survey 
#str(bottom_temp_comp) # high resolution BT model 
#levels(as.factor(bottom_temp_comp$Var))
# Select the hi res annual bottom temp 
BT <- bottom_temp_comp %>%
  dplyr::filter(Var==c("Annual_Bottom Temp"))

# join to weights
BT <- BT%>%
  left_join(EPUwt, by="EPU")

BT.An <- BT%>%
  dplyr::group_by(Time)%>%
  dplyr::mutate(BT = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::rename(year=Time) %>%
  dplyr::select(year,BT) %>%
  unique()

BT.Reg <- BT%>%
  dplyr::select(Time, EPU, Value) %>%
  dplyr::rename(year=Time) %>%
  dplyr::group_by(year, EPU)%>%
  dplyr::summarise(Annual.Mean = mean(Value))
  
# Calculate area weighted, annual basin mean and add to an object with columns for each EPU
BT.All <- BT.Reg %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Annual.Mean)  %>%
  dplyr::rename(BT.GB=GB,
                BT.GOM=GOM,
                BT.MAB=MAB,
                BT.SS=SS)%>%
  merge(BT.An, by = "year", all.x = TRUE, all.y = T)

#     C. Plot   -------------------------------------------------------
# Oisst (regional & annual --> patterns similar across region)
ts.OISST.Reg <- ggplot(OISST.Reg, aes(x = year, y = Annual.Mean, col = EPU))  + 
  theme_bw() +  
  geom_point(na.rm=T, size=1.5)+
  geom_line(na.rm = T, linewidth = 1)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("OISST")

ts.OISST.An <- ggplot(OISST.All, aes(x = year, y = OISST))  + 
  theme_bw() +  
  geom_point(na.rm=T, size=1.5)+
  geom_line(na.rm = T, linewidth = 1)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("OISST")

ts.BT.Reg <- ggplot(BT.Reg, aes(x = year, y = Annual.Mean, col = EPU))  +
  theme_bw()+
  geom_point(na.rm=T, size=1.5)+
  geom_line(na.rm = T, linewidth = 1)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("BT")

ts.BT.An <- ggplot(BT.An, aes(x=year, y=BT))+
  geom_point(na.rm=T, size=1.5)+
  geom_line(na.rm = T, linewidth = 1)+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("BT") 

# Now combine and save
tiff("Figures/ts.Temp.tiff", width = 6, height = 6, units = 'in', res = 300)
ggarrange(ts.OISST.Reg, ts.OISST.An, ts.BT.Reg, ts.BT.An, ncol=2, nrow = 2, widths=c(1, .7))
dev.off()

# 2. Environmental Indices -----------------------------------------------------------------
# Slope water proportions
#levels(as.factor(slopewater$Var))
WSW <- ecodata::slopewater %>%
  dplyr::filter(Var%in%c("WSW proportion ne channel"))%>% # select the warm slope water variable
  dplyr::select(Time, Value)%>%
  dplyr::rename(year=Time,
                WSW = Value)

# cold pool index (from ecodata) only MAB, various indices from 1958-2021 (recently changed calculation methods)
#str(cold_pool)
#levels(as.factor(cold_pool$Var))
CP <- ecodata::cold_pool %>%
  dplyr::filter(Var%in%c("cold_pool_index"))%>% # select the cold pool index variable
  dplyr::select(Time, Value)%>%
  dplyr::rename(year=Time,
                CP = Value)

# gulf stream index (from ecodata) from 1954-2020 for all regions (updated)
GSI <- ecodata::gsi
GSI$Year.Mo <- as.character(GSI$Time)

#split up time column into year and month
GSI <- GSI %>%
  as.data.frame()%>%
  tidyr::separate(Year.Mo, c("Year","Month"), sep = 4, remove = T, fill = "left", convert = T, extra = "merge") 

# compute annual average from monthly means --> might want to do this for just certain months?
GSI.Annual <- GSI %>%
  dplyr::group_by(Year) %>%
  summarise(GSI = mean(Value))%>%
  rename(year=Year)

# marine heatwave intensity (from ecodata) by epu, multiple variables 
### NOTE THAT THERES A NEW ONE CALLED ESP that may be better
str(ecodata::heatwave)
levels(as.factor(ecodata::heatwave$Var)) # 6 indices, cumulative intensity and max intensity

# add area weights and calculate basin wide annual average
# select cumulative intensity index
HW <- ecodata::heatwave %>% 
  dplyr::filter(Var%in%c("cumulative intensity-SurfaceDetrended"))%>% 
  left_join(EPUwt, by="EPU")

HW.An <- HW %>%
  dplyr::select(Time, EPU, Value, Weight) %>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(HW = weighted.mean(Value, Weight, na.rm = TRUE))%>% 
  dplyr::rename(year=Time)

# Pivot and rename for merge
HW.Reg <- HW %>%
  dplyr::select(Time, EPU, Value) %>%
  dplyr::group_by(Time) %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Value)  %>%
  dplyr::rename(HW.GB=GB,
                HW.GOM=GOM,
                HW.MAB=MAB,
                year = Time) %>%
  dplyr::select(!HW.MAB)

# Winter NAO from Hurrell website (use this one)
#https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
winterNAO_PC <- read.delim(here::here('Data/WinterNAO.txt'), header = T, sep = "\t", dec = ".")
winterNAO_PC <- winterNAO_PC %>%
  tidyr::separate(Hurrell.PC.Based.North.Atlantic.Oscillation.Index..DJFM.,
           c("year","NAOIndex"), remove = T,  sep="   ", convert = T)

#     A. Plot -------------------------------------------------------------------
# Cold Pool (annual), Gulf Stream (annual), Heatwave Index (annual and regional), Western Slope Water (annual), WinterNao
c <- ggplot(CP, aes(x = year, y = CP))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  ylab("Cold Pool Index")

g <- ggplot(GSI.Annual, aes(x = year, y = GSI))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  ylab("Gulf Stream Index")

h <- ggplot(HW, aes(x = Time, y = Value))  + 
  theme_bw() +  
  xlab("Year") +
  facet_wrap(~EPU)+
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  ylab("Heatwave Index")

h2 <- ggplot(HW.An, aes(x = year, y = HW))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  ylab("Heatwave Index")

w <- ggplot(WSW, aes(x=year, y=WSW))+
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  ylab("Warm Slope Water %")

nao <- ggplot(winterNAO_PC, aes(x = year, y = NAOIndex))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2025,5),minor_breaks = seq(0,2025,1), limits = c(1965,2023))+
  ylab("Winter NAO Index")

# Save these as one combined plot
tiff("Figures/Environmental Indices.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(c, g, h2, h, w, nao, ncol=2, nrow=3, widths=c(.7,1))
dev.off()

# 3. Upper Trophic ---------------------------------------------------------

# A. Predators 
# (downloaded from stock assessment portal and manually combined)
# Haddock SSB (Retrospective Adjusted from 2022 Assessment (for both GoM and GB)
# Should I add mackerel back in (don't think so) but def try with Micah's index instead of both Had separate
predators <- read.csv(here::here('Data/Pred_Assessment.csv'))
#predators <- predators%>% 
  #dplyr::select(year, HadSSB_GB, HadSSB_GoM)%>% 
  #dplyr::mutate(HadSSB_GoM=log(HadSSB_GoM),HadSSB_GB=log(HadSSB_GB))

#Plot
ts.had.GoM  <- ggplot(predators, aes(x = year, y = HadSSB_GoM))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("log GoM Haddock SSB (metric tons)")

ts.had.GB  <- ggplot(predators, aes(x = year, y = HadSSB_GB))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("GB Haddock SSB (metric tons)")

#ts.mackSSB  <- ggplot(predators, aes(x = year, y = MackSSB))  + 
  #theme_bw() +  
  #geom_line(linewidth = 1, na.rm = T) +
  #xlab("Year") +
  #ylab("Atlantic Mackerel SSB (metric tons)")

tiff("Figures/ts.Predators.tiff", width = 5, height = 5, units = 'in', res = 300)
ggarrange(ts.had.GoM, ts.had.GB)
dev.off()

#B. Competitors
# Make object from zooplankton anomalies indicator and add area weights
zoodat <- ecodata::zoo_abundance_anom %>% 
  full_join(EPUwt, by="EPU")%>%
  dplyr::rename(year=Time) %>%
  mutate(Value=as.numeric(Value))

# Cnidaria by EPU
jelly.Reg <- zoodat%>% 
  dplyr::filter(Var==c("Cnidaria"))%>%
  dplyr::select(year, EPU, Value)%>%
  dplyr::group_by(year) %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Value)  %>%
  dplyr::rename(jelly.GB=GB,
                jelly.GOM=GOM,
                jelly.MAB=MAB)

# Cnidaria basin wide --> this doesn't make sense b/c you shouldn't average an anomalie
jelly.An <- zoodat%>% 
  dplyr::filter(Var==c("Cnidaria"))%>%
  dplyr::group_by(year) %>%
  dplyr::summarize(Jelly.Abun = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::select(year, Jelly.Abun)

# Plot
ts.jelly.An <- ggplot(jelly.An, aes(x = year, y = Jelly.Abun))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Cnidaria Abundance Anoms")

ts.jelly.Reg <- zoodat%>% 
  dplyr::filter(Var==c("Cnidaria"))%>%
  dplyr::select(year, EPU, Value)%>%
  ggplot(aes(x = year, y = Value, group=EPU, col=EPU))  + 
  theme_bw() + 
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Cnidaria Abundance Anoms")

tiff("Figures/ts.Competitors.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(ts.jelly.An, ts.jelly.Reg)
dev.off()

# 4. Lower Trophic --------------------------------------------------------

#     A. Chlorophyll
# New indicator includes production by species 
# this must be cropped to seasonal shapefiles and includes each season
# for recruitment....chl from fall area in fall and winter is prob the most important
str(ecodata::ESP_seasonal_chl)
levels(ecodata::ESP_seasonal_chl$ESP)
#atlantic_herring_spring"
#atlantic_herring_fall" 

# Pull out chlorophyll data for the fall shapefile
Chl.fall <- ecodata::ESP_seasonal_chl %>%
  dplyr::filter(ESP==c("atlantic_herring_fall"))

# Split out into the seasons
Chl.fall.season <- Chl.fall %>% 
  dplyr::select(Time, Var, Value)%>%
  tidyr::pivot_wider(names_from = Var, values_from = Value) %>%
  dplyr::rename(Chl.fall.w = winter,
                Chl.fall.sp = spring,
                Chl.fall.su = summer,
                Chl.fall.f = fall) %>%
  dplyr::rename(year=Time)

ts.Chl.falldist  <- ggplot(Chl.fall, aes(x = Time, y = Value, group = Var, col = Var))  + 
  theme_bw() +  
  geom_line(na.rm = T, linewidth = 1) +
  xlab("Year") +
  ylab("Primary Production (gC/m2/day)")

#     B. Zooplankton 
#         i. Calanus stages -------------------------------------------------------
# Load calanus stage data, save as object, modify and summarize
head(ecodata::calanus_stage)
levels(as.factor(ecodata::calanus_stage$Var))

# split up by season and stage
c5fall <- ecodata::calanus_stage %>%
  dplyr::filter(Var==c("c5-Fall"))

# Pull out GoM timeseries of fall stage c5 
c5fall_GoM <- c5fall%>%
  dplyr::filter(EPU==c("GOM"))%>%
  #mutate(Value=log(Value)) %>%
  dplyr::rename(year=Time,
                c5fall.gom=Value)%>%
  dplyr::select(year, c5fall.gom)

#### Plots
 # plot of all data in stage stuff
ts.Cal.Stages <- ggplot(ecodata::calanus_stage, aes(x=Time, y = Value, group=EPU, color=EPU))+
  geom_line(size=1, na.rm=T)+
  facet_grid(.~Var)+
  theme_bw()

#  Stage C5 calanus in fall with colors for EPU
ts.Cal.Fall <- ggplot(c5fall, aes(x=Time, y = Value, group=EPU, color=EPU))+
  geom_line(na.rm=T)+
  geom_point(na.rm=T)+
  #facet_grid(.~Var)+
  ggtitle("Fall")+
  theme_bw()

# Plot just ts of calanus c5 in GoM fall
ts.C5.FalGom <- ggplot(c5fall_GoM, aes(x=year, y = c5fall.gom))+
  geom_line(na.rm=T)+
  geom_point(na.rm=T)+
  theme_bw()+
  xlab("Year") +
  ylab("Calanus finmarchicus C5 GoM Fall")

# plot together and save
#tiff("TS.Seas.Reg.CalStages.tiff", width = 8, height = 8, units = 'in', res = 300)
tiff("Figures/ts.prod&calanus.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ts.Chl.falldist, ts.Cal.Stages, ts.Cal.Fall, ts.C5.FalGom)
# Summer: lots of variation in all stages, but most noticeable in c5 summer, lots of interannual var that is diff in different regions
# Spring: most variation is in GoM, random spring blooms (mostly in copepodites)
# Fall: very low var in other stages, c5 shows lots of var in GoM and steady decline since 2000 --> use just that
dev.off()

#         ii. All Zooplankton  -------------------------------------------------------------
str(ecodata::zoo_strat_abun) # region and taxa (euphausids & cnidaria)
levels(as.factor(ecodata::zoo_strat_abun$EPU)) # this is only mab now 
str(ecodata::zoo_abundance_anom) # abundance anomalies for various zoo taxa
levels(as.factor(ecodata::zoo_abundance_anom$Var)) #calfin, smCopepods, LgCopepods, cnidaria (used above already)

# Pull out taxa of interest (from object made previously for cnidaria) & Pivot to a column for each region
calfin.Reg <- zoodat%>% 
  dplyr::filter(Var==c("Calfin"))%>%
  dplyr::select(year, EPU, Value)%>%
  dplyr::group_by(year) %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Value)  %>%
  dplyr::rename(cfin.GB=GB,
                cfin.GOM=GOM,
                cfin.MAB=MAB)

# Calculate area weighted basin wide annual mean
calfin.An <- zoodat%>% 
  dplyr::filter(Var==c("Calfin"))%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(cfin = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::select(year, cfin)

# Repeat for small and large copepod variables
smcope.Reg <- zoodat%>% 
  dplyr::filter(Var==c("SmCopepods"))%>%
  dplyr::select(year, EPU, Value)%>%
  dplyr::group_by(year) %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Value)  %>%
  dplyr::rename(smcope.GB=GB,
                smcope.GOM=GOM,
                smcope.MAB=MAB)

smcope.An <- zoodat%>% 
  dplyr::filter(Var==c("SmCopepods"))%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(Sm.cope = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::select(year, Sm.cope)

lgcope.Reg <- zoodat%>% 
  dplyr::filter(Var==c("LgCopepods"))%>%
  dplyr::select(year, EPU, Value)%>%
  dplyr::group_by(year) %>%
  tidyr::pivot_wider(names_from = EPU, values_from = Value)  %>%
  dplyr::rename(lgcope.GB=GB,
                lgcope.GOM=GOM,
                lgcope.MAB=MAB)

lgcope.An <- zoodat%>% 
  dplyr::filter(Var==c("LgCopepods"))%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(Lg.cope = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::select(year, Lg.cope)

#             Plot Zooplankton Time Series --------------------------------------------
# Plot all taxa w colors for EPU
ts.ZooReg.Taxa <- ggplot(na.omit(zoodat), aes(year, Value, col=EPU))+
  geom_line(na.rm = T)+
  theme_bw() +  
  facet_wrap(~Var)+
  xlab("Year") +
  ylab("Abundance Anomalies")

# Plot regionally faceted and annual average for calfin
ts.calfin.Reg <- zoodat%>% 
  dplyr::filter(Var==c("Calfin"))%>%
  dplyr::select(year, EPU, Value)%>%
  ggplot(aes(x = year, y = Value, group=EPU, col=EPU))  + 
  theme_bw() + 
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("C. fin Abundance Anoms")

ts.calfin.An <- ggplot(calfin.An, aes(x = year, y = cfin))  + 
  theme_bw() + 
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("C. fin Abundance Anoms")

# Repeat for Large Copepods 
ts.LgCope.Reg <- zoodat%>% 
  dplyr::filter(Var==c("LgCopepods"))%>%
  dplyr::select(year, EPU, Value)%>%
  ggplot(aes(x = year, y = Value, group=EPU, col=EPU))  + 
  theme_bw() + 
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Lg Cope Abundance Anoms")

ts.LgCope.An <- ggplot(lgcope.An, aes(x = year, y = Lg.cope))  + 
  theme_bw() + 
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Sm Cope Abundance Anoms")

# Repeat for Small
ts.SmCope.Reg <- zoodat%>% 
  dplyr::filter(Var==c("SmCopepods"))%>%
  dplyr::select(year, EPU, Value)%>%
  ggplot(aes(x = year, y = Value, group=EPU, col=EPU))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Sm Cope Abundance Anoms")

ts.SmCope.An <- ggplot(smcope.An, aes(x = year, y = Sm.cope))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Sm Cope Abundance Anoms")

# Combine and save
#tiff("TS.ZooAbun.tiff", width = 8, height = 8, units = 'in', res = 300)
#ggarrange(ts.ZooReg, ts.ZooAn, ts.SmCope, ts.LgCope, nrow=3, ncol=2)
tiff("ts.Zoo.Anoms.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ts.calfin.An, ts.calfin.Reg, ts.SmCope.An, ts.SmCope.Reg,ts.LgCope.An, ts.LgCope.Reg, nrow=3, ncol=2)
dev.off()

# Load Recruitment Data (from ASAP model) ---------------------------------------------------
herring = structure(list(
  year = c(1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 
           1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 
           2015, 2016, 2017, 2018, 2019, 2020, 2021),
  recruits = c(4.864754e+006, 3.913303e+006, 8.418418e+006, 4.027622e+006, 5.201490e+006, 2.475948e+006, 1.467439e+007, 2.562144e+006, 2.866499e+006, 3.507926e+006, 2.132604e+006, 2.565756e+006, 5.527145e+006, 5.771979e+006, 9.022276e+005, 3.973672e+006, 2.069876e+006, 1.921318e+006, 1.359330e+006, 4.504452e+006, 3.675050e+006, 3.223908e+006, 4.289285e+006, 6.276145e+006, 7.311144e+006, 
               7.794697e+006, 6.605043e+006, 3.360277e+006, 3.367537e+006, 3.707747e+006, 1.060684e+007, 4.650978e+006, 4.515168e+006, 2.864188e+006, 7.535893e+006, 2.209808e+006, 2.103963e+006, 4.495514e+006, 5.469233e+006, 2.949779e+006, 2.163047e+006, 4.575145e+006, 1.458040e+006, 2.773981e+006, 1.036147e+007, 1.995011e+006, 1.941889e+006, 6.115507e+006, 1.370991e+006, 1.316059e+006, 
               7.049114e+005, 3.435258e+005, 8.597470e+005, 6.927995e+005, 1.570989e+006, 8.637931e+005, 2.144506e+006),
  R.no.devs = c(3.051852e+006, 3.051852e+006, 3.052789e+006, 3.052839e+006, 3.052007e+006, 3.049923e+006, 3.047958e+006, 3.044345e+006, 3.042181e+006, 3.048009e+006, 3.044520e+006, 3.038618e+006, 3.033022e+006, 3.019695e+006, 2.998879e+006, 3.014066e+006, 3.009700e+006, 3.010080e+006, 3.001299e+006, 3.003036e+006, 3.002436e+006, 3.016398e+006, 3.022155e+006, 3.030740e+006, 3.030537e+006, 
                3.033298e+006, 3.033234e+006, 3.037982e+006, 3.044090e+006, 3.045528e+006, 3.043841e+006, 3.041957e+006, 3.040509e+006, 3.048207e+006, 3.045852e+006, 3.043410e+006, 3.044128e+006, 3.045739e+006, 3.042645e+006, 3.039240e+006, 3.040606e+006, 3.041675e+006, 3.037564e+006, 3.036944e+006, 3.037121e+006, 3.028771e+006, 3.026175e+006, 3.036045e+006, 3.040205e+006, 3.036605e+006, 
                3.042903e+006, 3.039353e+006, 3.030069e+006, 3.020596e+006, 2.999424e+006, 2.986980e+006, 2.984785e+006),
  logR.dev = c(4.662674e-001, 2.486331e-001, 1.014366e+000, 2.771040e-001, 5.331457e-001, -2.084930e-001, 1.571632e+000, -1.724412e-001, -5.948348e-002, 1.405364e-001, -3.559993e-001, -1.691496e-001, 6.001118e-001, 6.478591e-001, -1.201127e+000, 2.764007e-001, -3.743519e-001, -4.489552e-001, -7.920529e-001, 4.054424e-001, 2.021429e-001, 6.653083e-002, 3.501500e-001, 7.279491e-001, 8.806600e-001, 
               9.437933e-001, 7.782042e-001, 1.008297e-001, 1.009797e-001, 1.967501e-001, 1.248378e+000, 4.245765e-001, 3.954175e-001, -6.226842e-002, 9.058968e-001, -3.200731e-001, -3.693919e-001, 3.893366e-001, 5.864113e-001, -2.987703e-002, -3.405391e-001, 4.082298e-001, -7.339628e-001, -9.056839e-002, 1.227184e+000, -4.175076e-001, -4.436381e-001, 7.002721e-001, -7.963913e-001, -8.360980e-001, 
               -1.462495e+000, -2.180138e+000, -1.259703e+000, -1.472469e+000, -6.467147e-001, -1.240685e+000, -3.306184e-001),
  SR.std.resid = c(5.600442e-001, 2.986388e-001, 1.218378e+000, 3.328358e-001, 6.403733e-001, -2.504257e-001, 1.887722e+000, -2.071230e-001, -7.144694e-002, 1.688014e-001, -4.275987e-001, -2.031694e-001, 7.208078e-001, 7.781581e-001, -1.442700e+000, 3.319911e-001, -4.496425e-001, -5.392502e-001, -9.513524e-001, 4.869860e-001, 2.427984e-001, 7.991167e-002, 4.205730e-001, 8.743560e-001, 1.057780e+000, 
                   1.133611e+000, 9.347185e-001, 1.211089e-001, 1.212890e-001, 2.363210e-001, 1.499455e+000, 5.099683e-001, 4.749448e-001, -7.479199e-002, 1.088093e+000, -3.844469e-001, -4.436849e-001, 4.676409e-001, 7.043518e-001, -3.588598e-002, -4.090292e-001, 4.903340e-001, -8.815791e-001, -1.087837e-001, 1.473998e+000, -5.014778e-001, -5.328637e-001, 8.411125e-001, -9.565635e-001, -1.004256e+000, 
                   -1.756636e+000, -2.618612e+000, -1.513057e+000, -1.768615e+000, -7.767836e-001, -1.490215e+000, -3.971132e-001),
  SSB = c(9.753591e+005, 1.289719e+006, 1.312431e+006, 1.016283e+006, 6.492272e+005, 4.840740e+005, 3.296430e+005, 2.766896e+005, 4.872840e+005, 3.348256e+005, 2.187105e+005, 1.644389e+005, 1.031168e+005, 6.488890e+004, 8.900813e+004, 8.044368e+004, 8.112356e+004, 6.783045e+004, 7.010820e+004, 6.930487e+004, 9.436356e+004, 1.107674e+005, 1.492925e+005, 1.480771e+005, 1.664762e+005, 
          1.659991e+005, 2.108224e+005, 3.223843e+005, 3.681304e+005, 3.155974e+005, 2.721486e+005, 2.460963e+005, 5.001770e+005, 3.802812e+005, 3.044773e+005, 3.234533e+005, 3.759482e+005, 2.865523e+005, 2.270181e+005, 2.476787e+005, 2.666653e+005, 2.059334e+005, 1.990849e+005, 2.009925e+005, 1.382856e+005, 1.260161e+005, 1.899175e+005, 2.412469e+005, 1.955221e+005, 2.923659e+005, 
          2.286043e+005, 1.453515e+005, 1.057933e+005, 6.552942e+004, 5.344101e+004, 5.174865e+004, 5.656580e+004),
  SPR = c(2.073434e-001, 2.296337e-001, 2.192815e-001, 2.258184e-001, 2.103850e-001, 2.320578e-001, 2.317367e-001, 2.624353e-001, 2.829608e-001, 2.424901e-001, 2.518418e-001, 2.497410e-001, 2.239522e-001, 2.179583e-001, 2.558753e-001, 2.464404e-001, 2.508879e-001, 2.729209e-001, 2.929284e-001, 2.730506e-001, 2.597250e-001, 2.340829e-001, 2.177183e-001, 2.092220e-001, 2.227996e-001, 
                  1.998174e-001, 1.885841e-001, 1.911144e-001, 1.853427e-001, 1.637801e-001, 1.648042e-001, 1.629544e-001, 1.960076e-001, 1.621234e-001, 1.478445e-001, 1.759988e-001, 1.748467e-001, 1.579626e-001, 1.652241e-001, 1.751374e-001, 1.801409e-001, 1.617166e-001, 1.839529e-001, 1.785504e-001, 1.764589e-001, 1.642347e-001, 1.487846e-001, 1.508338e-001, 1.525817e-001, 1.717259e-001, 
                  1.667766e-001, 1.435470e-001, 1.587506e-001, 1.902731e-001, 1.950997e-001, 1.856265e-001, 1.593317e-001)),
  .Names = c("year", "recruits", "R.no.devs", "logR.dev", "SR.std.resid", "SSB", "SPR"),
  row.names = c("1965", "1966", "1967", "1968", "1969", "1970", "1971", "1972", "1973", "1974", "1975", "1976", "1977", "1978", "1979", "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989",
                "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014",
                "2015", "2016", "2017", "2018", "2019", "2020", "2021"),
  class = "data.frame")

#Calculate Spawning Success Rs (lnRt/ssbt-1)
herring$Rs <- NA
for(i in 2:nrow(herring)){
  herring$Rs[i] <- log(herring$recruits[i]/herring$SSB[i-1])
}

# Add lagged SSB column
herring$SSB_1 <- NA
for(i in 2:nrow(herring)){
  herring$SSB_1[i] <- herring$SSB[i-1]
}

# save as csv
#write.csv(herring, file = "ASAPdat.csv", row.names = F)

# Plot herring recruitment indices  -------------------------------------------------------------------
rec <- ggplot(herring, aes(x = year, y = recruits))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Recruitment")

nodev <- ggplot(herring, aes(x = year, y = R.no.devs))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Recruitment (No devs)")
#AR1 means?

logdevs <- ggplot(herring, aes(x = year, y = logR.dev))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Log Recruitment Deviations")

srresid <- ggplot(herring, aes(x = year, y = SR.std.resid))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Stock-Recruit Residuals (Standardized)")

ssb <- ggplot(herring, aes(x = year, y = SSB))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("SSB")

SPR <-  ggplot(herring, aes(x = year, y = SPR))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Spawning potential ratio")

Rs <- ggplot(herring, aes(x = year, y = Rs))  + 
  theme_bw() +  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Spawning Success (Rt/SSBt-1")

# Combine and save
#tiff("Recruitment Indices.tiff", width = 8, height = 6, units = 'in', res = 300)
ggarrange(ssb, rec, SPR, srresid,  nodev, Rs, logdevs)
#dev.off()

# NEW DATA ----------------------------------------------------------------
######## Load in Sarah's plankton indices 
springcalfinindex <- readRDS(here::here('Data/springcalfinindex.rds'))
springlgcopepodALLindex <- readRDS(here::here('Data/springlgcopepodALLindex.rds'))
springlg_wide <- springlgcopepodALLindex |> tidyr::pivot_wider(names_from = Var, values_from = Value)
fallsmcopepodALLindex <- readRDS(here::here('Data/fallsmcopepodALLindex.rds'))

#fallsmcopepodSOEindex <- readRDS(here::here('Data/fallsmcopepodSOEindex.rds'))
#fallsm_wide_SOE <- fallsmcopepodSOEindex |> tidyr::pivot_wider(names_from = Var, values_from = Value)

# Look at plots to decide which of these to select 

# 1. fall small (1 year lag to represent food in early larval stage) 
# Sarah has both the soe version and her version with diff spp and strata?
# Should see if this matches my calanus c5 in gom fall w one year lag b/c represents similar process
# plot
ggplot(fallsmcopepodALLindex, aes(x=Time, y =Value, col=EPU, group=EPU))+
  geom_line(na.rm=T)+
  geom_point(na.rm=T)+
  facet_wrap(~Var)

# Select areas and vars
fall.sm <- fallsmcopepodALLindex %>%
  dplyr::filter(Var =="Fall Small copepods ALL Abundance Index Estimate") %>%
  dplyr::filter(EPU=="her_fa")  %>%
  dplyr::rename(year=Time,
                fall.sm = Value) %>%
  #dplyr::mutate(fall.sm=log(fall.sm)) %>%
  dplyr::select(year, fall.sm)

# replace the 0 with na
fall.sm[fall.sm==0] <- NA
#as.na(fall.sm) <- sapply(fall.sm, is.infinite)

# 2. small spring (really winter) - don't have, maybe should although by winter they are late larvae
# 3. large spring - try calanus spring for this

# plot
springcalfinindex  %>%
  dplyr::filter(Var=="Spring Calanus finmarchicus Abundance Index Estimate")  %>%
  ggplot(aes(x=Time, y =Value, col=EPU, group=EPU))+
  geom_line(na.rm=T)+
  geom_point(na.rm=T)+
  facet_wrap(~EPU)

#also try sarahs large b/c it includes more spp
springlgcopepodALLindex  %>%
  dplyr::filter(Var=="Spring Large copepods ALL Abundance Index Estimate" )  %>%
  ggplot(aes(x=Time, y =Value, col=EPU, group=EPU))+
  geom_line(na.rm=T)+
  geom_point(na.rm=T)+
  facet_wrap(~EPU)
# these look pretty much the same

# select area and index
spring.lg <-   springlgcopepodALLindex  %>%
  dplyr::filter(Var=="Spring Large copepods ALL Abundance Index Estimate" )  %>%
  dplyr::filter(EPU=="her_sp")  %>%
  dplyr::rename(year=Time,
                spring.lg = Value)  %>%
  #dplyr::mutate(spring.lg=log(spring.lg)) %>%
  dplyr::select(year, spring.lg)

# and new model outputs from mm 192
# import updated recruitment timeseries
#mm192 <- readRDS(here::here('Data/mm192/mm192.rds'))
# somehow need to extract the first column of the NAA object 
mm192 <- readRDS(here::here('Data/mm192/res_tables/stock_1_region_1_NAA_table.rds'))
str(mm192)
names(mm192_1)
mm192_1 <- as.data.frame(mm192)
mm192_1$year <- as.numeric(row.names(mm192_1))
WHAM <- mm192_1  %>%
  dplyr::select(year, '1') %>%
  dplyr::rename(recrt = '1')%>%
  dplyr::mutate(logR = log(recrt),
                logRdev = logR - mean(logR))

ggplot(WHAM, aes(x=year, y=logRdev))+
  theme_bw()+  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("log Recruitment deviations")

ggplot(WHAM, aes(x=year, y=recrt))+
  theme_bw()+  
  geom_line(na.rm=T, linewidth=1) +
  geom_point(na.rm=T, size=1.5)+
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +
  ylab("Recruitment")

# "Smart" Lags ------------------------------------------------------------

##################### Temperature data 
# SST: 0 lag for larvae and juveniles (just GoM & GB)
# SST: up to 4 years for prespawn adults (all areas)
# BT: 1 year lag for egg & YSL temperatures (GoM/GB only)
Temp_lags <- merge(BT.All, OISST.All, by = "year", all.x = TRUE, all.y = T)%>%
  arrange(year) %>%
  mutate(OISST_1 = lag(OISST,1),
         #OISST_2 = lag(OISST,2),
         #OISST_3 = lag(OISST,3),
         #OISST_4 = lag(OISST,4),
         #BT_1 = lag(BT,1),
         #BT_2 = lag(BT,2),
         BT.GOM_1 = lag(BT.GOM, 1),
         BT.GB_1 = lag(BT.GB, 1))%>%
  dplyr::select(c(year, BT.GB_1 ,BT.GOM_1,OISST.GB, OISST.GOM)) # only select what we want

###################### Environmental indices
#(not including slope water, gsi or cold pool, or hw (missing years) maybe should, b/c one year lag on water bodies may impact early larval drift in fall
# winter NAO: Up to 5 years (see Groger)
# HW: no lag (juvenile conditions), 1 year lag (early larvae) (just GB/GoM)
Env_lags <- merge(HW.Reg, winterNAO_PC, by = "year", all.x = TRUE, all.y = T) %>%
  mutate(NAOIndex_1= lag(NAOIndex,1),
         NAOIndex_2= lag(NAOIndex,2),
         NAOIndex_3= lag(NAOIndex,3),
         #NAOIndex_4= lag(NAOIndex,4),
         #NAOIndex_5= lag(NAOIndex,5),
         HW.GB_1 = lag(HW.GB, 1),
         HW.GOM_1 = lag(HW.GOM,1))

###################### Predation
# benthic egg predation, 1 year lag for both haddock stocks
# predation and competition with larvae, 1 year for jellies in GB & GOM
Pred_lags <- predators %>%
  mutate(HadSSB_GB_1 = lag(HadSSB_GB,1),
         HadSSB_GoM_1 = lag(HadSSB_GoM,1))

###################### Food 
# Large copepods: 1 - 3 year lag for adults (only spring, add large fall to represent pre-spawning feeding)
# Large copepods: 0 year for juveniles and late larvae in spring
# Small copepods: 1 year lag for early larvae in the fall
Food_lags <-  merge(spring.lg, fall.sm, by = "year", all.x = TRUE, all.y = T) %>%
  arrange(year) %>%
  mutate(#spring.lg_1 = lag(spring.lg, 1),
         #spring.lg_2 = lag(spring.lg, 2),
         #spring.lg_3 = lag(spring.lg, 3),
         fall.sm_1 = lag(fall.sm, 1))

Food_lags%>% 
  reshape2::melt(id.vars="year", variable.name="series")%>%
  ggplot(aes(year,value)) +
  geom_line(na.rm = T)+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free")

# Calanus finmarchicus stage 5 in fall in the GoM (species and stage specific early larval food in one area)
c5fall_GoM_lag <- c5fall_GoM %>%
  mutate(c5fall.gom_1 = lag(c5fall.gom, 1))

# Combine  Physical Variables ---------------------------------------------

# Full list of ALL possible w lags
#env.full.lag <- list(Env_lags, Temp_lags, Prod_lags, WSW, CP) # work on this
#might want to also include the regionals but I don't have lags for those ugh

# Full list of ALL possible without lags
#env.full <- list(BT.All, OISST.All, HW.An, HW.Reg,GSI.Annual, CP, WSW, winterNAO_PC) 
#env.full <- env.full %>% reduce(full_join, by='year')  %>% arrange(year)

# Pared down list
#env.vars <- list(BT.An, BT.An_1, OISST.Annual, OISST.An_1, HW.An, GSI.Annual, WSW) # no cold pool?

# Chop off extra years
#env.vars <- env.vars %>% reduce(full_join, by='year')

# Combine Biological Variables ------------------------------------------------------

# Full list of ALL possible w lags --> work on this
eco.full.lag <- list(predators, predators_1, 
                        Chla.An, PPD.An, 
                        jelly, smallcal, smallcal_1, ZooAbun.An, ZooAbun.An_1, cfin.GoM.C5, cfin.GoM.C5_1)
eco.full.lag <- eco.full.lag %>% reduce(full_join, by='year')

# Full list of ALL possible without lags --> not updated with Sarah's indices
eco.full <- list(Chl.fall.season, predators,
                 jelly.An, jelly.Reg,
                 c5fall_GoM, calfin.An, calfin.Reg,
                 lgcope.An, lgcope.Reg,
                 smcope.An, smcope.Reg)

eco.full <- eco.full %>% reduce(full_join, by='year') %>%
  arrange(year)

# Join variables & herring data  -----------------------------------------
# adding new data - combine only smart lags and possibly important vars
smart_vars <- list(Temp_lags, Food_lags, Pred_lags, c5fall_GoM_lag , Env_lags) 
smart_vars <- smart_vars %>% reduce(full_join, by='year')

##### 
#vars.full <- merge(eco.full, env.full, by = "year", all.x = TRUE, all.y = T)
#variables <- merge(eco.vars, env.vars, by = "year", all.x = TRUE, all.y = T)

# Export and save
#write.csv(variables, file = "Data/Variables_smartlag.csv", row.names = F)

# Combine the data called herring and variables into one dataframe
#dat.full <- merge(herring, vars.full, by = "year",  all.x = TRUE, all.y = T)
#combined.dat <-merge(herring, variables, by = "year",  all.x = TRUE, all.y = T)
newdat <- merge(WHAM, smart_vars, by = "year",  all.x = TRUE, all.y = T)

# chop it down
#combined.dat <- combined.dat %>% filter(between(year,1982,2019))
#dat.full <- dat.full %>% filter(between(year,1982,2021))
newdat <- newdat %>% filter(between(year,1987,2023))

# plot all
newdat%>% 
  reshape2::melt(id.vars="year", variable.name="series")%>%
  ggplot(aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free")
# Trim off columns we don'nt need/want (mostly SS and MAb stuff)
#newdat <- newdat%>%dplyr::select(-HW.MAB, -BT.MAB, -BT.SS,-OISST.MAB,  -OISST.SS, -fall.sm)
# Export to reload easily into quarto/rmd
#write.csv(combined.dat4, file = "Data/Combined_lags.csv", row.names = F)

# look at correlations
do <- get_upper_tri(cor(as.data.frame(newdat), use="na.or.complete"))
cormat <- reshape2::melt(do, na.rm=T)
cor1 <- ggplot(data = cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
cor1     


ggplot2::ggplot(ecodata::storminess, aes(x=Year, y=Value, col=Var))+
  geom_point(na.rm=T)+
  geom_line(na.rm=T) +
  facet_wrap(~Var)


