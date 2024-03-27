##      Title:            Get, prepare, combine, and plot environmental & biological variables
##      Script Purpose:   Prepare multivariate time series for analysis
##      Author:           Adelle Molina
##      Created:          7/22/22
##      Updated:          2/19/24

#Change these to the following and add all data to the folder
#read.csv(here::here('data/catch_data/gt_data_model_cpue.csv'))

# Libraries ---------------------------------------------------------------
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

# Bring in the Regional Areas (to use for area weighted means)
load("Data/StratAreas.Rdata") # from techdoc

# Calculate the weights
strat.area <- strat.area%>%
  dplyr::mutate(Weight=Area/sum(Area))

# A. OISST -------------------------------------------------------------------
# Manually downloaded daily mean SST (1982-present) from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# In "get sst data" script" computed daily mean by EPU, combined, and exported, load that in --> now cannot find that
OISST  <- read.csv(file.choose(), header = T, blank.lines.skip = T)

# join to weights
OISST <- OISST%>%
  left_join(strat.area, by="EPU")

# Calculate basin wide annual mean
OISST.Annual  <- OISST %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(OISST = weighted.mean(Value, Weight, na.rm = TRUE))

# Calculate annual means by region
OISST.Regional <- OISST %>%
  dplyr::select(year, EPU, day, Value) %>%
  dplyr::group_by(year, EPU) %>%
  dplyr::summarise(Annual.Mean = mean(Value))

# Combine these into one with columns for each region
OISST.reg <- pivot_wider(OISST.Regional, names_from = EPU, values_from = Annual.Mean) 
OISST.reg <- OISST.reg %>%
  dplyr::rename(OISST.GB=GB,
                OISST.GOM=GOM,
                OISST.MAB=MAB,
                OISST.SS=SS)
# B. BT  -------------------------------------------------------
str(bottom_temp) #in situ anomalies from survey 
str(bottom_temp_comp) # high resolution BT model 
levels(as.factor(bottom_temp_comp$Var))

# Select the hi res annual bottom temp (by epu)
BT <- bottom_temp_comp %>%
  dplyr::filter(Var==c("Annual_Bottom Temp"))

# join to weights
BT <- BT%>%
  left_join(strat.area, by="EPU")

# Calculate area weighted, annual basin mean
BT.An <- BT%>%
  dplyr::select(Time, EPU, Value, Area, Weight)%>%
  dplyr::group_by(Time)%>%
  dplyr::summarise(BT = weighted.mean(Value, Weight, na.rm = TRUE))%>%
  dplyr::rename(year=Time)

# Plot Temperature  -------------------------------------------------------
# Oisst (regional & annual --> patterns similar across region)
ts.OISST.Reg <- ggplot(OISST.Regional, aes(x = year, y = Annual.Mean, col = EPU))  + 
  theme_bw() +  
  geom_line(na.rm = T, size = 1.5)+
  xlab("Year") +
  ylab("OISST")

ts.OISST.An <- ggplot(OISST.Annual, aes(x = year, y = OISST))  + 
  theme_bw() +  
  geom_line(na.rm = T, size = 1.5)+
  xlab("Year") +
  ylab("OISST")

ts.bt.An <- ggplot(BT.An, aes(x = year, y = BT))  + 
  theme_bw() +   
  geom_line(na.rm = T, size = 1.5) +
  xlab("Year") +
  ylab("BT")

ts.BT.Reg <- ggplot(BT, aes(x=Time, y=Value, col=EPU, group=EPU))+
  geom_line(na.rm = T, size = 1.5)+
  theme_bw() +
  xlab("Year") +
  ylab("BT") 

# Now combine and save
tiff("ts.Temp.tiff", width = 6, height = 6, units = 'in', res = 300)
ggarrange(ts.OISST.Reg, ts.OISST.An, ts.BT.Reg, ts.bt.An, ncol=2, nrow = 2, widths=c(1, .7))
dev.off()

# Environmental Indices -----------------------------------------------------------------

# Slope water proportions
str(slopewater)
levels(as.factor(slopewater$Var))
WSW <- slopewater %>%
  dplyr::filter(Var%in%c("WSW proportion ne channel"))%>% # select the warm slope water variable
  dplyr::select(Time, Value)%>%
  dplyr::rename(year=Time,
                WSW = Value)

# cold pool index (from ecodata) only MAB, various indices from 1958-2021 (recently changed calculation methods)
str(cold_pool)
levels(as.factor(cold_pool$Var))
CP <- cold_pool %>%
  dplyr::filter(Var%in%c("cold_pool_index"))%>% # select the cold pool index variable
  dplyr::select(Time, Value)%>%
  dplyr::rename(year=Time,
                CP = Value)

# gulf stream index (from ecodata) from 1954-2020 for all regions (updated)
GSI <- gsi
GSI$Year.Mo <- as.character(GSI$Time)

#split up time column into year and month
GSI <- GSI %>%
  as.data.frame()%>%
  separate(Year.Mo, c("Year","Month"), sep = 4, remove = T, fill = "left", convert = T, extra = "merge") 

# compute annual average
GSI.Annual <- GSI %>%
  dplyr::group_by(Year) %>%
  summarise(GSI = mean(Value))%>%
  rename(year=Year)

# marine heatwave intensity (from ecodata) by epu, multiple variables 
str(heatwave)
HW <- heatwave
levels(as.factor(heatwave$Var)) # two indices, cumulative intensity and max intensity

# add area weights and calculate basin wide annual average
HW <- HW %>% left_join(strat.area, by="EPU")
HW.An <- HW %>%
  dplyr::filter(Var%in%c("cumulative intensity"))%>% # select cumulative intensity index
  dplyr::group_by(Time) %>%
  dplyr::summarise(HW = weighted.mean(Value, Weight, na.rm = TRUE))%>% 
  dplyr::rename(year=Time)

# annual regional means
HW.Reg <- heatwave %>%
  dplyr::filter(Var%in%c("cumulative intensity"))%>%
  dplyr::group_by(Time, EPU) %>%
  dplyr::summarise(HW = mean(Value, na.rm=T))%>%
  dplyr::rename(year=Time)%>%
  as.data.frame()

# Pivot and rename for merge
HW.reg <- HW.Reg%>%
  tidyr::pivot_wider(names_from = EPU, values_from = HW)  %>%
  dplyr::rename(HW.GB=GB,
                HW.GOM=GOM,
                HW.MAB=MAB)
# lots of N/A in this set b/c not available for all regions in all years

# Winter NAO from Hurrell website (use this one)
#https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
winterNAO_PC <- read.delim(file.choose(), header = T, sep = "\t", dec = ".")
winterNAO_PC <- winterNAO_PC %>%
  separate(Hurrell.PC.Based.North.Atlantic.Oscillation.Index..DJFM.,
           c("year","NAOIndex"), remove = T,  sep="   ", convert = T)

# Plot environmental indices -------------------------------------------------------------------
# Cold Pool (annual), Gulf Stream (annual), Heatwave Index (annual and regional), Western Slope Water (annual), WinterNao
c <- ggplot(CP, aes(x = year, y = CP))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T) +
  ylab("Cold Pool Index")

g <- ggplot(GSI.Annual, aes(x = year, y = GSI))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T) +
  ylab("Gulf Stream Index")

h <- ggplot(HW.Reg, aes(x = year, y = HW))  + 
  theme_bw() +  
  xlab("Year") +
  facet_wrap(~EPU)+
  geom_line(na.rm=T) +
  ylab("Heatwave Index")

h2 <- ggplot(HW.An, aes(x = year, y = HW))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T) +
  ylab("Heatwave Index")

w <- ggplot(WSW, aes(x=year, y=WSW))+
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T) +
  ylab("Warm Slope Water %")

nao <- ggplot(winterNAO_PC, aes(x = year, y = NAOIndex))  + 
  theme_bw() +  
  xlab("Year") +
  geom_line(na.rm=T) +
  xlim(1965,2023)+
  ylab("Winter NAO Index")

# Save these as one combined plot
tiff("Environmental Indices.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(c, g, h2, h, w, nao, ncol=2, nrow=3, widths=c(.7,1))
dev.off()

# Load and modify Biological Data ---------------------------------------------------------

# Predator data 
# Atlantic mackerel SSB from 2021 assessment (downloaded from stock assessment portal)
# Haddock SSB (from John with recent update --> now using SSb from 2022 assessment)
#predators <- read.csv(file.choose(), header = T, blank.lines.skip = T)
#predators <- predators%>%dplyr::rename(year=Year)
predators <- read.csv(here::here('Data/Pred_Assessment.csv'))

# Primary Production (chlorophyll a annual area weighted, primary production regional and area weighted basin wide)
Primary <- chl_pp
levels(as.factor(Primary$Var))

# First split up the time variable and add area weights
Primary$Year <- as.character(Primary$Time)
PP <- Primary%>%
  separate(Year, c("Extra","year"), sep = 2, remove = T, fill = "left", convert = T, extra = "merge")
PP <- PP %>% left_join(strat.area, by="EPU")
PP$year <- as.integer(PP$year, na.rm=T)

# Chlorophyll a  from 1998-2020 regional
Chla <- PP%>%
  dplyr::filter(Var==c("ANNUAL_CHLOR_A_MEDIAN"))%>% 
  dplyr::select(year, EPU, Value, Var, Weight)%>%
  as.data.frame()

# Calculate annual basin mean
Chla.An <- Chla%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(ChlA=weighted.mean(Value, Weight, na.rm = TRUE))

# Primary Production Rate regional
PPD.Reg <- PP%>%
  dplyr::filter(Var==c("ANNUAL_PPD_MEDIAN"))%>% 
  dplyr::select(year, EPU, Value, Var, Weight)%>%
  as.data.frame()

# Pivot and rename for merging 
PPD.reg <- PPD.Reg %>% 
  dplyr::select(year, EPU, Value)%>%
  pivot_wider(names_from = EPU, values_from = Value) %>%
  dplyr::rename(PP.GB=GB,
                PP.GOM=GOM,
                PP.MAB=MAB)

# annual basin mean
PPD.An <- PPD.Reg%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(PPD=weighted.mean(Value, Weight, na.rm = TRUE))

# Plot production variables and predation -----------------------------------------------
ts.ChlAn  <- ggplot(Chla.An, aes(x = year, y = ChlA))  + 
  theme_bw() +  
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("Chlorophyll a (mg/m3)")

ts.PPDAn  <- ggplot(PPD.An, aes(x = year , y = PPD))  + 
  theme_bw() +  
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("Primary Production (gC/m2/day)")

ts.PPDReg  <- ggplot(PPD.Reg, aes(x = year, y = Value, group = EPU, col = EPU))  + 
  theme_bw() +  
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("Primary Production (gC/m2/day)")

# Combine and save
tiff("ts.PrimaryProd.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(ts.ChlAn, ts.PPDAn, ts.PPDReg, widths=c(1, .7))
dev.off()

# Predators
ts.hadSSB.gom  <- ggplot(predators, aes(x = year, y = HadSSB_GoM))  + 
  theme_bw() +  
  geom_line(size = 1, na.rm = T) +
  xlab("Year") +
  ylab("GoM Haddock SSB (metric tons)")

ts.hadSSB.gb  <- ggplot(predators, aes(x = year, y = HadSSB_GB))  + 
  theme_bw() +  
  geom_line(size = 1, na.rm = T) +
  xlab("Year") +
  ylab("GB Haddock SSB (metric tons)")

ts.mackSSB  <- ggplot(predators, aes(x = year, y = MackSSB))  + 
  theme_bw() +  
  geom_line(size = 1, na.rm = T) +
  xlab("Year") +
  ylab("Atlantic Mackerel SSB (metric tons)")

tiff("ts.Predators.tiff", width = 5, height = 5, units = 'in', res = 300)
ggarrange(ts.mackSSB, ts.hadSSB.gom, ts.hadSSB.gb)
dev.off()

# PLANKTON  ------------------------------------------------------------

# 1. Calanus stages -------------------------------------------------------
# Load calanus stage data, save as object, modify and summarize
load(file='Data/RM_20210120_CalanusStage.Rda')

calanus <- CalanusStage%>%
  dplyr::rename(EPU=epu)%>%
  left_join(strat.area, by="EPU")

# split up by season
Cal.Fall <- calanus%>%
  dplyr::filter(Units==c("No. per 100m^-3"))%>%
  dplyr::filter(season==c("Fall"))
Cal.Summer <- calanus%>%
  dplyr::filter(Units==c("No. per 100m^-3"))%>%
  dplyr::filter(season==c("Summer"))
Cal.Spring <- calanus%>%
  dplyr::filter(Units==c("No. per 100m^-3"))%>%
  dplyr::filter(season==c("Spring"))

# Pull out GoM timeseries of fall stage c5 
cfin.GoM.C5 <- Cal.Fall%>%
  dplyr::filter(Var==c("c5"))%>%
  dplyr::filter(EPU==c("GOM"))%>%
  mutate(Value=log(Value))
  
#           Plot Calanus finmarchicus time series -----------------------------------
#  seasonal trends for all stages with colors for EPU
ts.Cal.Fall <- ggplot(Cal.Fall, aes(x=Year, y = Value, group=EPU, color=EPU))+
  geom_line(size=1, na.rm=T)+
  facet_grid(.~Var)+
  ggtitle("Fall")+
  theme_bw()

ts.Cal.Spring <- ggplot(Cal.Spring, aes(x=Year, y = Value, group=EPU, color=EPU))+
  geom_line(size=1, na.rm=T)+
  theme_bw()+
  facet_grid(.~Var)+
  ggtitle("Spring")

ts.Cal.Summer <- ggplot(Cal.Summer, aes(x=Year, y = Value, group=EPU, color=EPU))+
  geom_line(size=1, na.rm=T)+
  theme_bw()+
  facet_grid(.~Var)+
  ggtitle("Summer")

ts.C5.FalGom <- ggplot(cfin.GoM.C5, aes(x=Year, y = Value))+
  geom_line(size=1, na.rm=T)+
  theme_bw()+
  xlab("Year") +
  ylab("Calanus finmarchicus C5 GoM Fall")

# plot together and save
tiff("TS.Seas.Reg.CalStages.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ts.Cal.Summer, ts.Cal.Spring, ts.Cal.Fall, ts.C5.FalGom)
# Summer: lots of variation in all stages, but most noticeable in c5 summer, lots of interannual var that is diff in different regions
# Spring: most variation is in GoM, random spring blooms (mostly in copepodites)
# Fall: very low var in other stages, c5 shows lots of var in GoM and steady decline since 2000 --> use just that
dev.off()

# 2.All Zooplankton Stratified Abundance -------------------------------------------------------------
str(zoo_strat_abun) # region and taxa (euphausids, cnidaria, large and small calanoida)
levels(as.factor(zoo_strat_abun$Var))

# Add weights
zoodat <- zoo_strat_abun %>% full_join(strat.area, by="EPU")
head(zoodat)
# Calculate area weighted basin wide annual mean for each taxa
smallcal <- zoodat%>% 
  dplyr::filter(Var==c("SmallCalanoida"))%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(Sm.cal.abun = log(weighted.mean(Value, Weight, na.rm = TRUE)))%>%
  dplyr::rename(year=Time) %>%
  dplyr::select(year, Sm.cal.abun)

largecal <- zoodat%>% 
  dplyr::filter(Var==c("LargeCalanoida"))%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(Lg.cal.abun = log(weighted.mean(Value, Weight, na.rm = TRUE)))%>%
  dplyr::rename(year=Time) %>%
  dplyr::select(year, Lg.cal.abun)

jelly <- zoodat%>% 
  dplyr::filter(Var==c("Cnidaria"))%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(Jelly.Abun = log(weighted.mean(Value, Weight, na.rm = TRUE)))%>%
  dplyr::rename(year=Time) %>%
  dplyr::select(year, Jelly.Abun)

krill <- zoodat%>% 
  dplyr::filter(Var==c("Euphausiacea"))%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(Krill.Abun = log(weighted.mean(Value, Weight, na.rm = TRUE)))%>%
  dplyr::rename(year=Time) %>%
  dplyr::select(year, Krill.Abun)

# regional total for all taxa
ZooAbun.Reg <- zoodat%>%
  group_by(Time, EPU)%>%
  dplyr::summarise(Abund = log(sum(Value)))

# Calculate annual basin scale, area weighted average of the total
ZooAbun.An <- zoodat%>%
  group_by(Time, EPU)%>%
  dplyr::summarise(Abund = sum(Value))

ZooAbun.An <- ZooAbun.An%>%
  full_join(strat.area, by="EPU")%>%
  group_by(Time)%>%
  dplyr::summarise(ZooAbun = log(weighted.mean(Abund, Weight, na.rm = TRUE)))%>%
  dplyr::rename(year=Time)%>%
  dplyr::select(year, ZooAbun)

#           Plot Zooplankton Time Series --------------------------------------------
# Panels for taxonomic groups with colors for EPU (raw data scale)
ts.ZooReg.Taxa <- ggplot(na.omit(zoodat), aes(Time, Value, col=EPU))+
  geom_line(na.rm = T, size=1)+
  theme_bw() +  
  facet_wrap(~Var)+
  xlab("Year") +
  ylab("Zooplankton Abundance (#)")

# Total Zooplankton Abundance by Region
ts.ZooReg <- ggplot(na.omit(ZooAbun.Reg), aes(x = Time, y = Abund, col=EPU, group=EPU))  + 
  theme_bw() +  
  geom_line(na.rm = T, size=1)+
  xlab("Year") +
  ylab("Log Zooplankton Abundance (#)")

# Total Basin w3ide Zooplankton Abundance (area weighted mean)
ts.ZooAn <- ggplot(ZooAbun.An, aes(x = year, y = ZooAbun))  + 
  theme_bw() +  
  geom_line(size=1, na.rm = T)+
  xlab("Year") +
  ylab("Log Zooplankton Abundance (#)")

# Mean Annual Large Calanoida Abundance (area weighted)
ts.LgCope <- ggplot(largecal, aes(x = year, y = Lg.cal.abun))  + 
  theme_bw() + 
  geom_line(size=1, na.rm = T)+
  xlab("Year") +
  ylab("Log Abundance of Large Calanoida")

# Mean Annual Small Calanoida Abundance (area weighted)
ts.SmCope <- ggplot(smallcal, aes(x = year, y = Sm.cal.abun))  + 
  theme_bw() +  
  geom_line(size=1, na.rm = T)+
  xlab("Year") +
  ylab("Log Abundance of Small Calanoida")

ts.jelly <- ggplot(jelly, aes(x = year, y = Jelly.Abun))  + 
  theme_bw() +  
  geom_line(size=1, na.rm = T)+
  xlab("Year") +
  ylab("Log Abundance of Cnidarians")

ts.krill <- ggplot(krill, aes(x = year, y = Krill.Abun))  + 
  theme_bw() +  
  geom_line(size=1, na.rm = T)+
  xlab("Year") +
  ylab("Log Abundance of Euphausids")
# Combine and save
tiff("TS.ZooAbun.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ts.ZooReg, ts.ZooAn, ts.SmCope, ts.LgCope, ts.jelly, ts.krill, nrow=3, ncol=2)
dev.off()

# Add lags ------------------------------------------------------
# Temperature data (up to 3 years)
Temp <- merge(BT.An, OISST.Annual, by = "year", all.x = TRUE, all.y = T)%>%
  arrange(year) %>%
  mutate(OISST_1 = lag(OISST,1),
         OISST_2 = lag(OISST,2),
         OISST_3 = lag(OISST,3),
         BT_1 = lag(BT,1),
         BT_2 = lag(BT,2),
         BT_3 = lag(BT,3))

#Add lags to climate indices 
#(not including slope water, gsi and cold pool, Up to 5 years for winter NAO (see Groger))
Env <- merge(HW.An, winterNAO_PC, by = "year", all.x = TRUE, all.y = T) %>%
  arrange(year) %>%
  mutate(HW_1 = lag(HW,1),
         HW_2 = lag(HW,2),
         HW_3 = lag(HW,3),
         NAOIndex_1= lag(NAOIndex,1),
         NAOIndex_2= lag(NAOIndex,2),
         NAOIndex_3= lag(NAOIndex,3),
         NAOIndex_4= lag(NAOIndex,4),
         NAOIndex_5= lag(NAOIndex,5))

# Production
Prod <- merge(Chla.An, PPD.An, by = "year", all.x = TRUE, all.y = T) %>%
  arrange(year) %>%
  mutate(ChlA_1 = lag(ChlA,1),
         ChlA_2 = lag(ChlA,2),
         ChlA_3 = lag(ChlA,3),
         PPD_1= lag(PPD,1),
         PPD_2= lag(PPD,2),
         PPD_3= lag(PPD,3))

# Add "smart" lags 

# 1. for benthic egg predation, 1 year lag for both mackerel and haddock
predators1 <-predators
predators1$year <- predators1$year+1
colnames(predators1)[2:ncol(predators1)]<-paste(colnames(predators)[2:ncol(predators)],"_1",sep="")
predators1<-as.data.frame(predators1)
predators1<-predators1%>%
  add_row(year=1962, .before=1)%>%
  dplyr::filter((year<=2020)%>% replace_na(TRUE))

# 2. Adult food, total zoo abundance, 1 year lag --> no longer using this
#ZooAbun.An_1 <- data.frame(year=ZooAbun.An$year+1,ZooAbun_1=ZooAbun.An$ZooAbun)
#ZooAbun.An_1 <- ZooAbun.An_1[-c(nrow(ZooAbun.An_1)),] 

# 2a. Use this instead
largecal <-  largecal %>%
  arrange(year) %>%
  mutate(Lg.cal.abun_1 = lag(Lg.cal.abun, 1),
         Lg.cal.abun_2 = lag(Lg.cal.abun, 2),
         Lg.cal.abun_3 = lag(Lg.cal.abun, 3))

# 2. Early larvae food, small copepods (include both 0 and 1 year lag)
smallcal <-  smallcal %>%
  arrange(year) %>%
  mutate(Sm.cal.abun_1 = lag(Sm.cal.abun, 1))

# 3. Calanus finmarchicus stage 5 in fall in the GoM (early larval food in one area)
cfin.GoM.C5 <- cfin.GoM.C5 %>%
  dplyr::rename(year=Year) %>%
  dplyr::select(year, Value) %>%
  dplyr::rename(Cfin.C5=Value)%>%
  arrange(year) %>%
  mutate(Cfin.C5_1 = lag(Cfin.C5, 1))

# Join variables  -----------------------------------------

ecological.vars <- list(predators, predators1, Prod,
                        jelly, smallcal, largecal, cfin.GoM.C5) 
eco.vars <- ecological.vars %>% reduce(full_join, by='year')

env.vars <- list(Env, Temp, GSI.Annual, WSW, CP)
env.vars <- env.vars %>% reduce(full_join, by='year')

variables <- merge(eco.vars, env.vars,
                   by = "year", 
                   all.x = TRUE, all.y = T)

# Export and save
write.csv(variables, file = "Data/Variables.csv", row.names = F)

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

#Calculate Spawning Success Rs (lnRt/ssb-1)
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
write.csv(herring, file = "ASAPdat.csv", row.names = F)

# Plot herring recruitment indices  -------------------------------------------------------------------
rec <- ggplot(herring, aes(x = year, y = recruits))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Recruitment")

nodev <- ggplot(herring, aes(x = year, y = R.no.devs))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Recruitment (No devs)")
# is this estimated from stock recruit function

logdevs <- ggplot(herring, aes(x = year, y = logR.dev))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Log Recruitment Deviations")

srresid <- ggplot(herring, aes(x = year, y = SR.std.resid))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Stock-Recruit Residuals (Standardized)")

ssb <- ggplot(herring, aes(x = year, y = SSB))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("SSB")

SPR <-  ggplot(herring, aes(x = year, y = SPR))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Spawning potential ratio")

Rs <- ggplot(herring, aes(x = year, y = Rs))  + 
  theme_bw() +  
  geom_point(na.rm = F)+
  geom_line(size=1) +
  xlab("Year") +
  ylab("Spawning Success (Rt/SSBt-1")

# Combine and save
tiff("Recruitment Indices.tiff", width = 8, height = 6, units = 'in', res = 300)
ggarrange(ssb, rec, SPR, srresid,  nodev, Rs, logdevs)
dev.off()


# Correlatin matrix of variables
variables2 <- variables %>% 
  dplyr::select(OISST, BT, GSI, HW, WSW, CP, NAOIndex, 
                Lg.cal.abun, Sm.cal.abun, Jelly.Abun, Cfin.C5,
                MackSSB, HadSSB_GOM, HadSSB_GB)

ggstatsplot::ggcorrmat(
  data = variables2,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))
