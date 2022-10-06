##      Title:            Get, prepare, combine, and plot environmental & biological variables
##      Script Purpose:   Prepare multivariate time series for analysis, explore correlations
##      Author:           Adelle Molina
##      Created:          7/22/22
##      Updated:          10/3/22
##      To Do:            1. Area weighted annual means --> done for survey data but not for other variables, which need it?
##                        3. Organize autocorrelation section
##      Notes:            Some saved plots include unweighted averages
# Libraries ---------------------------------------------------------------

#devtools::install_github("slucey/RSurvey/Survdat")
#remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
library(ecodata)
library(Survdat)
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


library(lubridate)
library(CCA)

# Load and modify environmental data ---------------------------------------------------
# OISST: Manually downloaded daily mean SST (1982-present, missing 1981 for some reason) from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# In separate script, computed daily mean in each EPU, combined, and exported, load that back in

OISST  <- read.csv(file.choose(), header = T, blank.lines.skip = T)

# Add a column for season (Need justification for these cutoff days, I selected them myself)
OISST <- OISST %>%
  dplyr::mutate(Season = case_when(day <= 59 ~ "Winter",
                                  day >59 & day <=181~ "Spring", 
                                  day >181 & day <=243~ "Summer",
                                  day >243 ~ "Fall")) %>%
  mutate_if(is.character, as.factor)

# Calculate seasonal means 
OISST.Seasonal <- OISST %>%
  dplyr::select(Season, year, Value) %>%
  dplyr::group_by(year, Season) %>%
  dplyr::summarise(Seasonal.Mean = mean(Value))

# Calculate annual mean, not weighted
OISST.Annual <- OISST%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(OISST = mean(Value))

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

# cold pool index (from ecodata) only MAB, various indices from 1958-2021
CP <- cold_pool %>%
  dplyr::filter(Var%in%c("cold_pool_index"))%>% # select the cold pool index variable
  dplyr::select(Time, Value)%>%
  dplyr::rename(year=Time,
         CP = Value)

# gulf stream index (from ecodata) from 1954-2020 for all regions 
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

# marine heatwave intensity (from ecodata) by epu from 1982-2020, multiple variables
HW <- heatwave
HW.Annual <- HW %>%
  dplyr::filter(Var%in%c("cumulative intensity"))%>% # select cumulative intensity as the index value
  dplyr::group_by(Time) %>%
  dplyr::summarise(HW = mean(Value))%>% # calculate annual mean across all regions, should this be weighted?
  dplyr::rename(year=Time)

# annual means by region
HW.Regional <- heatwave %>%
  dplyr::filter(Var%in%c("cumulative intensity"))%>%
  dplyr::group_by(Time, EPU) %>%
  dplyr::summarise(HW = mean(Value, na.rm=T))%>%
  dplyr::rename(year=Time)%>%
  as.data.frame()

# Pivot and rename for merge
# oops was fucking around with this to get rid of nas but it was fine all along
HW.reg <- HW.Regional%>%
  tidyr::pivot_wider(names_from = EPU, values_from = HW)  %>%
  dplyr::rename(HW.GB=GB,
                HW.GOM=GOM,
                HW.MAB=MAB)

# Compute Surface and Bottom Temp and Sal from NEFSC survey data ----------------------------------------------------------

# Survey data (2 forms)
#1.  Anomalies (load in data called NEFSC Temp Anomalies)

# From ecodata github, regional bottom temperature, edited in separate script and exported as a csv, load that back in
NEFSC  <- read.csv(file.choose(), header = T, blank.lines.skip = T) # survey temperature anomalies

# calculate annual mean bottom temp by region --> add the reference to the anomaly
NMFS.BT <- NEFSC %>%
  dplyr::select(Time, EPU, Var, Value) %>%
  dplyr::filter(Var%in%c("bottom temp anomaly in situ", "reference bt in situ (1981-2010)"))%>% 
  dplyr::group_by(Time, EPU) %>%
  dplyr::mutate(ActualBT = Value + Value[Var==c("bottom temp anomaly in situ")])%>%
  dplyr::filter(!Var==c("bottom temp anomaly in situ")) %>%
  dplyr::select(Time, EPU, ActualBT)

# Repeat for all regions (this should be area weighted, or come from the annual anomaly model)
NMFS.BT.AnAll <- NMFS.BT%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(BT = mean(ActualBT))%>%
  dplyr::rename(year=Time)

# calculate annual mean SST by region  -->  add the reference to the anomaly
NMFS.SST <- NEFSC %>%
  select(Time, EPU, Var, Value) %>%
  dplyr::filter(Var%in%c("sst anomaly in situ", "reference sst in situ (1981-2010)"))%>% 
  dplyr::group_by(Time, EPU) %>%
  mutate(ActualSST = Value + Value[Var==c("sst anomaly in situ")])%>%
  dplyr::filter(!Var==c("sst anomaly in situ")) %>%
  dplyr::select(Time, EPU, ActualSST)

# Repeat for all regions
NMFS.SST.AnAll <- NMFS.SST%>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(SST = mean(ActualSST))%>%
  rename(year=Time)

#2.  RAW
# Raw surface and bottom temp, sal from Survey 
# From survdat, trimmed to EPU in another script (it's called Survdat EPU)
NEFSC.RAW <- read.csv(file.choose(), header = T, blank.lines.skip = T) 

# Summarize by season, year, and EPU
Survey <- NEFSC.RAW%>%
  dplyr::group_by(YEAR, SEASON, EPU) %>%
  dplyr::summarise(SST = mean(SURFTEMP, na.rm = TRUE),
                   SSS = mean(SURFSALIN, na.rm = TRUE),
                   BT = mean(BOTTEMP, na.rm = TRUE),
                   BS = mean(BOTSALIN, na.rm = TRUE))

# Split by season --> didn't use this yet
Fall <- Survey%>%
  dplyr::filter(SEASON==c("FALL"))

Spring <- Survey%>%
  dplyr::filter(SEASON==c("SPRING"))  

# Summarize across both seasons (should these be weighted by sampling frequency within region?)
Survey.Annual <- NEFSC.RAW%>%
  dplyr::group_by(YEAR, EPU) %>%
  dplyr::summarise(SST = mean(SURFTEMP, na.rm = TRUE),
                   SSS = mean(SURFSALIN, na.rm = TRUE),
                   BT = mean(BOTTEMP, na.rm = TRUE),
                   BS = mean(BOTSALIN, na.rm = TRUE))

# Compute annual bottom temperature, area weighted
#Generate area table
EPU_sf <- spTransform(EPU, crs)
strat.area <- Survdat::getarea(EPU_sf, 'EPU')

# calculate area weights
NEFSC.weight   <- Survdat::stratprep(as.data.table(NEFSC.RAW),   strat.area, strat.col = 'EPU', area.col = 'Area')


# Calculate area weighted annual values
Survey.Annual.weighted <- NEFSC.weight%>% 
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(Surv.SST = weighted.mean(SURFTEMP, W.h, na.rm = TRUE),
                   Surv.SSS = weighted.mean(SURFSALIN, W.h,na.rm = TRUE),
                   Surv.BT = weighted.mean(BOTTEMP, W.h,na.rm = TRUE),
                   Surv.BS = weighted.mean(BOTSALIN, W.h,na.rm = TRUE))%>%
  dplyr::rename(year=YEAR)

# Calculate annual values
#Survey.AnAll <- NEFSC.RAW%>% 
  #dplyr::group_by(YEAR) %>%
  #dplyr::summarise(Surv.SST = weighted.mean(SURFTEMP, weights$W.h, na.rm = TRUE),
                   #Surv.SSS = mean(SURFSALIN, na.rm = TRUE),
                   #Surv.BT = mean(BOTTEMP, na.rm = TRUE),
                   #Surv.BS = mean(BOTSALIN, na.rm = TRUE))%>%
  dplyr::rename(year=YEAR)

# Create object with columns for epu for sst and bt
SurvSST.Regional <-  Survey.Annual%>%
  dplyr::select(YEAR, EPU, SST)%>%
  pivot_wider(names_from = EPU, values_from = SST) %>%
  dplyr::rename(year=YEAR,
                SST.GB=GB,
                SST.GOM=GOM,
                SST.MAB=MAB,
                SST.SS=SS)

SurvBT.Regional <-  Survey.Annual%>%
  dplyr::select(YEAR, EPU, BT)%>%
  pivot_wider(names_from = EPU, values_from = BT) %>%
  dplyr::rename(year=YEAR,
                BT.GB=GB,
                BT.GOM=GOM,
                BT.MAB=MAB,
                BT.SS=SS)

# Plot environmental indices -------------------------------------------------------------------
# Cold Pool (annual), Gulf Stream (annual), Heatwave Index (annual and regional)
c <- ggplot(CP, aes(x = year, y = CP))  + 
  theme_bw() +  
  geom_point()+
  xlab("Year") +
  geom_line(na.rm=T, size = 2) +
  ylab("Cold Pool Index MAB")

g <- ggplot(GSI.Annual, aes(x = year, y = GSI))  + 
  theme_bw() +  
  geom_point()+
  xlab("Year") +
  geom_line(na.rm=T, size = 2) +
  ylab("Gulf Stream Index")

h <- ggplot(HW.Regional, aes(x = year, y = HW))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  facet_wrap(~EPU)+
  geom_line(na.rm=T, size = 2) +
  ylab("Heatwave Index")

h2 <- ggplot(HW.Annual, aes(x = year, y = HW))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  geom_line(na.rm=T, size = 2) +
  ylab("Heatwave Index")

# Save these as one combined plot
tiff("Environmental Indices.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(c, g, h2, h)
dev.off()

# Plot SST  -------------------------------------------------------
# Oisst (annual, regional --> poss not necessary, patterns look similar across epu)
# seasonal oisst (Not done)
# Survey sst annual and regional (done, from both raw and anomalies data versions)
ois <- ggplot(OISST.Regional, aes(x = year, y = Annual.Mean, col = EPU))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("OISST")

o <- ggplot(OISST.Annual, aes(x = year, y = OISST))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("OISST")

sstanom <- ggplot(NMFS.SST, aes(x = Time, y = ActualSST, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ SST (from anoms)")

sstanom2 <- ggplot(NMFS.SST.AnAll, aes(x = year, y = SST))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ SST (from anoms)")

sstraw <- ggplot(Survey, aes(x = YEAR, y = SST, col = EPU, shape = SEASON))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ SST")

# repeat the above without the season variable
sstraw2 <- ggplot(Survey.Annual, aes(x = YEAR, y = SST, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ SST")
# This plot looks different from what should be the same plot based on anomaly data, even the one with seasons has same patterns....hmmm what's going on
# maybe b/c 2019 no fall samples...?

sstraw3 <- ggplot(Survey.AnAll, aes(x = year, y = Surv.SST))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ SST")

sstraw.weighted <- ggplot(Survey.Annual.weighted, aes(x = year, y = Surv.SST))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ SST")
ggarrange(sstraw3, sstraw.weighted)

# Now combine and save
tiff("SST.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ois, o, sstanom, sstanom2, sstraw2, sstraw.weighted, sstraw, ncol=2, nrow = 4)
dev.off()
# Should remove the last year and replot b/c no fall 2019 data

# Plot Bottom Temp --------------------------------------------------------
# Survey bottom temp Annual & Regional (done, from both raw and anomalies data versions)
# seasonal bt (Not DONE, use raw survey data)

bt <- ggplot(NMFS.BT, aes(x = Time, y = ActualBT, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=3)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ BT (from anoms)")

bt1 <- ggplot(NMFS.BT.AnAll, aes(x = year, y = BT))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=3)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ BT (from anoms)")

btraw <- ggplot(Survey, aes(x = YEAR, y = BT, col = EPU, shape = SEASON))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ BT")

# repeat the above without the season variable
btraw2 <- ggplot(Survey.Annual, aes(x = YEAR, y = BT, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ BT")

btraw3 <- ggplot(Survey.AnAll, aes(x = year, y = Surv.BT))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("In Situ BT")

# Combine and save
tiff("BT.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(bt, bt1, btraw2, btraw3, btraw, ncol = 2, nrow=3 ) 
dev.off()

# Now plot just the anomalies (only for bottom temperature, prob don't need these for all metrics)
btanom <- NEFSC %>%
  dplyr::group_by(Time, EPU) %>%
  dplyr::filter(Var==c("bottom temp anomaly in situ"))

btref<- NEFSC %>%
  dplyr::group_by(Time, EPU) %>%
  dplyr::filter(Var==c("reference bt in situ (1981-2010)"))

bottom.anoms <- ggplot(btanom, aes(x = Time, y = Value, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=3)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("Bottom Temp Anomaly")

bottom.refs <- ggplot(btref, aes(x = Time, y = Value, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=3)+
  geom_line(na.rm = T, size = 2) +
  xlab("Year") +
  ylab("Bottom Temp Reference Value")

# Plot Salinity -----------------------------------------------------------
# Only from survey
surfs <- ggplot(Survey.Annual, aes(x = YEAR, y = SSS, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ Surface Sal")

surfs2 <- ggplot(Survey.AnAll, aes(x = year, y = Surv.SSS))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ Surface Sal")

bs <- ggplot(Survey.Annual, aes(x = YEAR, y = BS, col = EPU))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ Bottom Sal")

bs2 <- ggplot(Survey.AnAll, aes(x = year, y = Surv.BS))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ Bottom Sal")

# Combine and save
tiff("Salinity.tiff", width = 5, height = 5, units = 'in', res = 300)
ggarrange(bs, bs2, surfs, surfs2)
dev.off()

# Combine  Physical Variables ---------------------------------------------
# regional sst (done for raw only; objects called surv )
# regional bt (done for raw only; objects called surv)
# Don't use anomalies data b/c they are not the reanalized version
# join the selected variables together
physvar.list <- list(OISST.Annual, OISST.reg, 
                     NMFS.BT.AnAll, NMFS.SST.AnAll,
                     Survey.Annual.weighted, SurvBT.Regional, SurvSST.Regional, 
                     HW.Annual, HW.reg, GSI.Annual, CP)
physvar.list <- physvar.list %>% reduce(full_join, by='year')

# Trim down to the year range we want 
physvar.list <- physvar.list %>%
  dplyr::filter(between(year,1982,2018))

# Load and modify Biological Data ---------------------------------------------------------

# Chlorophyll has regional annual medians from 1998-2019
Primary <- chl_pp
levels(as.factor(Primary$Var))# lots of diff vars to choose from 
Primary$Year <- as.character(Primary$Time)
PP <- Primary%>%
  separate(Year, c("Extra","year"), sep = 2, remove = T, fill = "left", convert = T, extra = "merge")
PP <- PP%>%
  dplyr::filter(Var==c("ANNUAL_CHLOR_A_MEDIAN"))%>% 
  dplyr::select(year, EPU, Value)%>%
  as.data.frame()

PP.Ann <- PP%>%
  dplyr::group_by(year) %>%
  dplyr::summarise(PP=mean(Value, na.rm = TRUE)) # should these be area weighted
PP.Ann$year <- as.numeric(PP.Ann$year)

# Pivot and rename for merge
PP.reg <- pivot_wider(PP, names_from = EPU, values_from = Value) 
PP.reg <- PP.reg %>%
  dplyr::rename(PP.GB=GB,
                PP.GOM=GOM,
                PP.MAB=MAB)
# Zooplankton  ------------------------------------------------------------

# Load calanus r data, save as object, modify and summarize
copepod <- calanus_stage
copepod <- copepod%>%
  mutate_if(is.character, as.factor)%>%
  dplyr::filter(EPU%in% c("GOM", "SS", "GB"))%>% # Select the spawning regions only
  dplyr::group_by(Time, Var) %>% # what are cutoff dates for seasons (6 total var's are stage and season combined)
  dplyr::summarise(Anmean = mean(Value, na.rm = TRUE))%>%
  as.data.frame()

# split up by season and stage, get regional annual averages for both stages
cope <- copepod %>%
  separate(Var, c("Stage","Season"), sep = "-", remove = T,  convert = T, extra = "merge")

# now summarize across all seasons for both adults,  stage c5, and both (later in the process might want to keep seasons)
cope.adult <- cope%>%
  dplyr::filter(Stage==c("adt"))%>%
  group_by(Time)%>%
  dplyr::summarise(CAD = mean(Anmean))%>%
  dplyr::rename(year=Time)

cope.c5 <- cope%>%
  dplyr::filter(Stage==c("c5"))%>%
  group_by(Time)%>%
  dplyr::summarise(CC5 = mean(Anmean))%>%
  dplyr::rename(year=Time)

cope.Ann <- copepod%>%
  group_by(Time)%>%
  dplyr::summarise(Cope = mean(Anmean))%>%
  dplyr::rename(year=Time)

cope.reg <- calanus_stage%>%
  dplyr::filter(EPU%in% c("GOM", "SS", "GB"))%>% # Select the spawning regions only
  group_by(Time, EPU)%>%
  dplyr::summarise(Cope = mean(Value))%>%
  dplyr::rename(year=Time)

# Pivot and rename for merge
c.reg <- pivot_wider(cope.reg, names_from = EPU, values_from = Cope) 
c.reg <- c.reg %>%
  dplyr::rename(cal.GB=GB,
                cal.GOM=GOM,
                cal.SS=SS)

# 5 other Rdata files for zooplankton
# 1. abundance
head(zoo_abund) # regional abundance
range(zoo_abund$Value) # weird this one has negative values, don't use
# 2. diversity (diversity index, nah)
# 3. oi
head(zoo_oi) # regional density by calanus family and season (centropages, pseudocalanus, temora) with sd
levels(as.factor(zoo_oi$Var))
# 4. anom (anomaly by region, nah)
# 5. strat abund by region and season
head(zoo_strat_abun) # regional abundance by taxa (euphausids, cnidaria, large and small calanoida)

# Load the abundance data and edit
zoo.numb <- zoo_strat_abun%>%
  group_by(Time, EPU, Var)%>%
  dplyr::summarise(Abund = mean(Value))%>%
  dplyr::rename(year=Time)

# Create another object with annual average across taxa
zoonum_reg <- zoo.numb%>%
  group_by(year, EPU)%>%
  dplyr::summarise(Abund = mean(Abund))

# Pivot and rename for merge
zn.reg <- pivot_wider(zoonum_reg, names_from = EPU, values_from = Abund) 
zn.reg <- zn.reg %>%
  dplyr::rename(Num.GB=GB,
                Num.GOM=GOM,
                Num.MAB=MAB)

# Calculate annual average
zoonum_an <- zoo_strat_abun%>%
  group_by(Time)%>%
  dplyr::summarise(Abund = mean(Value))%>%
  dplyr::rename(year=Time)

# Load the oi (regional density) data and edit
# First split up var column by species and season
zoo.den <- zoo_oi %>%
  separate(Var, c("Taxon","Season"), sep = "zoo", remove = T,  convert = T, extra = "merge")

zoo.den <- zoo.den%>%
  mutate_if(is.character, as.factor)%>%
  group_by(Time, EPU, Taxon, Season)%>%
  dplyr::summarise(Dens = mean(Value))%>%
  dplyr::rename(year=Time)

# Remove the sd factor levels in the taxon
zoo.den <- zoo.den[zoo.den$Season %in% c(" fall", " spring"), ]
zoo.den <- droplevels(zoo.den)

# Grouping 
# First get regional annual average for all Calanus taxa
zooden_reg <-zoo.den%>%
  group_by(year, EPU)%>%
  dplyr::summarise(Density = mean(Dens))

# repeat across all regions
zooden_an <-zoo.den%>%
  group_by(year)%>%
  dplyr::summarise(ZooDens = mean(Dens))

# Then compute average density by region by species (average across season))
zooden <- zoo.den%>%
  group_by(year, EPU, Taxon)%>%
  dplyr::summarise(Dens = mean(Dens))

# Plot biological variables -----------------------------------------------
prim  <- ggplot(PP, aes(x = year, y = Value, col = EPU))  + 
  theme_bw() +  
  geom_point() +
  geom_line(size = 2) +
  xlab("Year") +
  ylab("Primary Production (gC/m2/day")

# Zooplankton abundance
zooabun_an <- ggplot(zoo.num.an, aes(x = year, y = log(Abund)))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Abundance (log(Number)")

zooabun <- ggplot(zoo.num.reg, aes(x = year, y = log(Abund), col = EPU))  + # use log so scale is easily visible
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Abundance (log(Number)")

zooabun2 <- ggplot(zoo.numb, aes(x = year, y = log(Abund), col = Var))  + # use log so scale is easily visible
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Abundance (log(Number)")+
  facet_wrap(~EPU)

ZD <- ggplot(zoo.den2, aes(x = year, y = Dens, col = Taxon))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  facet_wrap(vars(EPU))+
  ylab("Zooplankton Density (#/m2/day)")

ZDr <- ggplot(zooden_reg, aes(x = year, y = Density, col = EPU))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Density(#/m2/day)")

ZDa <- ggplot(zooden_an, aes(x = year, y = ZooDens))  + 
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Density(#/m2/day)")

# For copepods, plot adults and c5 stage separately
ad <- ggplot(cope.adult, aes(x = year, y = CAD))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year")+
  ylab("Adult Copepod Abundance")

c5 <- ggplot(cope.c5, aes(x = year, y = CC5))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("C5 Copepod Abundance")

# Copepods together (facet by stage with colors by season)
cal <- ggplot(cope, aes(x = Time, y = Anmean, col=Season))  + 
  theme_bw() +  
  geom_point(na.rm = T)+
  geom_line(na.rm = T) +
  xlab("Year") +
  ylab("Copepod Abundance")+
  facet_wrap(~Stage)

# Mean annual copepod abundance by region
cal2 <- ggplot(cope.reg, aes(x = year, y = Cope, col = EPU))  + 
  theme_bw() +  
  geom_point(na.rm = T)+
  geom_line(na.rm = T) +
  xlab("Year") +
  ylab("Copepod Abundance")

# Mean annual copepod abundance

cal3 <- ggplot(cope.Ann, aes(x = year, y = Cope))  + 
  theme_bw() +  
  geom_point(na.rm = T)+
  geom_line(na.rm = T) +
  xlab("Year") +
  ylab("Copepod Abundance")

# Combine and save
tiff("Biological Variables.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(zooabun_an, zooabun, zooabun2, ZDa, ZDr, ZD, cal3, c5, ad, cal, cal2, prim, nrow=4, ncol=3)
dev.off()

# Combine Biological Variables ------------------------------------------------------

# Annual and regional primary productivity 
# Annual and regional zooplankton abundance 
# Annual zooplankton density
# Annual copepod (total, adult, and c5) & regional (total)
# Taxon specific zoo abundance/density (not done)
# Seasonal copepod (not done)

biovar.list <- list(zoonum_an, zn.reg, zooden_an,
                    cope.adult, cope.c5, cope.Ann, c.reg,
                    PP.Ann, PP.reg)

biovar.list <- biovar.list %>% reduce(full_join, by='year')

# Trim down to the year range we want 
biovar.list <- biovar.list  %>% 
  dplyr::filter(between(year,1982,2018))

# Load Recruitment Data (from ASAP model run) ---------------------------------------------------

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

# save as csv
write.csv(herring, file = "ASAPdat.csv", row.names = F)

# Plot herring recruitment indices  -------------------------------------------------------------------

rec <- ggplot(herring, aes(x = year, y = recruits))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Recruitment")

nodev <- ggplot(herring, aes(x = year, y = R.no.devs))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Recruitment (No devs)")
# is this estimated from stock recruit function

devs <- ggplot(herring, aes(x = year, y = logR.dev))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Log Recruitment Deviations")

srresid <- ggplot(herring, aes(x = year, y = SR.std.resid))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Stock-Recruit Residuals (Standardized)")

ssb <- ggplot(herring, aes(x = year, y = SSB))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("SSB")

SPR <-  ggplot(herring, aes(x = year, y = SPR))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Spawner per Recruit")

Rs <- ggplot(herring, aes(x = year, y = Rs))  + 
  theme_bw() +  
  geom_point(na.rm = F, size=3)+
  geom_line() +
  xlab("Year") +
  ylab("Spawning Success")

# Combine and save
tiff("Recruitment Indices.tiff", width = 8, height = 6, units = 'in', res = 300)
ggarrange(ssb, rec, SPR, srresid, devs, nodev, Rs)
dev.off()

# Join biological and environmental covariates to herring data -----------------------------------------
vars <- merge(physvar.list, biovar.list ,
                      by = "year", 
                      all.x = TRUE, all.y = T)
multivariate <- merge(herring, vars,
        by = "year", 
        all.x = TRUE, all.y = T)

# Trim off some years 
multivariate <- multivariate %>% 
  filter(between(year,1982,2018))

# Save as a csv
write.csv(multivariate, file = "combined.data.csv", row.names = F)

# Pull out primary production variables and salinity (need a shorter time series for those)
# also could creater a longer series with just survey temps
multi.short <- multivariate%>% 
  select(year, OISST, BT, SST, Surv.SST, Surv.BT, 
         BT.GOM, BT.GB, BT.MAB , BT.SS, 
         SST.GOM, SST.GB, SST.MAB , SST.SS, 
         GSI, CP, HW, HW.GB, HW.GOM, HW.MAB,
         Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,
         recruits, SPR, logR.dev, Rs)

# Create a smaller version of variable list (no redundant survey data (anomalies vs raw, just use raw))
var.list <- list(OISST.Annual, OISST.reg, Survey.Annual.weighted, SurvBT.Regional, SurvSST.Regional, 
                 HW.Annual, HW.reg, GSI.Annual, CP, PP.Ann, PP.reg,
                 zoonum_an, zn.reg, zooden_an,cope.adult, cope.c5, cope.Ann, c.reg) 
var.list <- var.list %>% reduce(full_join, by='year')
var.list <- var.list%>% 
  #dplyr::select(-Surv.SSS, -Surv.BS) %>% 
  dplyr::filter(between(year,1982,2018))

write.csv(var.list, file = "Data/variables.csv", row.names = F)

# Auto Correlation ---------------------------------------------------------

# convert full dataset to long form, 
multi.long <- multivariate %>% 
  pivot_longer(cols = c(OISST, BT, SST, Surv.SST, Surv.BT, 
                        BT.GOM, BT.GB, BT.MAB , BT.SS, 
                        SST.GOM, SST.GB, SST.MAB , SST.SS, 
                        GSI, CP, HW, HW.GB, HW.GOM, HW.MAB,
                        Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,
                        recruits, SPR, logR.dev, Rs),
               names_to = "names", values_to = "value")%>% 
  select(year, names, value)

#multi.edit <- multivariate %>% 
  #pivot_longer(cols = c(c(names(vars))[2:41], recruits, SPR, logR.dev, Rs), names_to = "names", values_to = "value")%>% 
  #select(year, names, value)

# Plot time series (looks terrible too many variables)
#multi.edit %>% 
  #plot_time_series(.date_var = year,
                   #.value = value,
                   #.facet_vars = names)
# acf/pacf plots
#multi.edit %>%
  #group_by(names) %>% 
  #plot_acf_diagnostics(.date_var = year,
                       #.value = value,
                       #.show_white_noise_bars = T)

acf(herring$SPR, type = "covariance")
acf(herring$SPR)
pacf(herring$SPR)
# looks like an AR(1) type of pattern
ccf(multivariate$SPR, multivariate$HW, plot = T, na.action = na.pass)

# acf/pacf plots
bio.melt %>%
  group_by(series)%>%
  plot_acf_diagnostics(.date_var = year,
                       .value = value,
                       .show_white_noise_bars = T, .facet_ncol = 2)

# Combined correlation matrix ---------------------------------------------

# Create correlation matrix
cormat <- melt(cor(as.data.frame(multivariate), use="na.or.complete"))

# Get lower triangle of the correlation matrix
get_lower_tri<-function(co){
  co[upper.tri(co)] <- NA
  return(co)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(co){
  co[lower.tri(co)]<- NA
  return(co)
}
do <- get_upper_tri(cor(as.data.frame(multivariate), use="na.or.complete"))
cormat <- melt(do, na.rm=T)
cor1 <- ggplot(data = cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


cor2pcor(m, tol)
pcor2cor(m, tol)

# Grouped correlation matrices --------------------------------------------

## Repeat for just SST variables and recruitment INDICES
temps <- multivariate%>% 
  select(year, OISST, SST, Surv.SST,  
         SST.GOM, SST.GB, SST.MAB , SST.SS, 
         recruits, SPR, logR.dev, Rs)

#do2 <- get_upper_tri(cor(as.data.frame(temps), use="na.or.complete", method="spearman"))
#cormat2 <- melt(do2, na.rm=T)
#cor2 <- ggplot(data = cormat2, aes(Var2, Var1, fill = value))+
  #geom_tile(color = "white")+
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
  #theme_minimal()+ 
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  #coord_fixed()

# correlogram 
cor.sst <-ggstatsplot::ggcorrmat(
  data = temps,
  type = "nonparametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)

# SPR 3 sig correlations with *SST (from anom), and GB/GOM SST (from survey), Rs one sig with GB sst

# From those results, look at autorrelation in an even more reduced set
temps.long <- temps %>% 
  pivot_longer(cols = c(OISST, SST, Surv.SST,  SST.GOM, SST.GB, SST.MAB , SST.SS, SPR),
               names_to = "names", values_to = "value")%>% 
  select(year, names, value)

# pacfs
temps.long %>% 
  plot_time_series(.date_var = year,.value = value,.facet_vars = names, .facet_ncol = 2, .smooth = F)
temps.long %>%
  group_by(names) %>% 
  plot_acf_diagnostics(.date_var = year, .value = value)

# Use built in acf to look at autocorrelation within the TS (prob don't need)
acf(herring$SPR, type = "covariance")
acf(lh, type = "covariance")
pacf(herring$SPR)
ccf(multivariate$logR.dev, multivariate$HW, plot = T, na.action = na.pass)

# BT variables and one recruitment index
BTS <- multivariate%>% 
  select(year, BT, Surv.BT, BT.GOM, BT.GB, BT.MAB , BT.SS, recruits, SPR, logR.dev, Rs)
#do3 <- get_upper_tri(cor(as.data.frame(BTS), use="na.or.complete"))
#cormat3 <- melt(do3, na.rm=T)

cor.bt <-ggstatsplot::ggcorrmat(
  data = BTS,
  type = "nonparametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)
# SPR 2 sig correlations with *GB/GOM BT

# Repeat for environmental indices
ENVS<- multivariate%>% 
  select(year, GSI, CP, HW, HW.GB, HW.GOM, HW.MAB,recruits, SPR, logR.dev, Rs)
#do4 <- get_upper_tri(cor(as.data.frame(ENVS), use="na.or.complete"))
#cormat4 <- melt(do4, na.rm=T)
cor.ind <-ggstatsplot::ggcorrmat(
  data = ENVS,
  type = "nonparametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)
# SPR 2 sig correlations with HW & GB HW*

## Repeat for food indices
FOODS <- multivariate%>% 
  select(year,   Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,
         recruits, SPR, logR.dev, Rs)
#do5 <- get_upper_tri(cor(as.data.frame(FOODS), use="na.or.complete"))
#cormat5 <- melt(do5, na.rm=T)
#cor5 <- ggplot(data = cormat5, aes(Var2, Var1, fill = value))+
  #geom_tile(color = "white")+
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       #midpoint = 0, limit = c(-1,1), space = "Lab", #name="Pearson\nCorrelation") +
  #theme_minimal()+ 
  #theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1))+
  #coord_fixed()

cor.food <-ggstatsplot::ggcorrmat(
  data = FOODS,
  type = "nonparametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)
# none

# Combine and save, use this to decide what to include in the version to be used in trees
ggarrange(cor.ind, cor.bt, cor.sst, cor.food)


###########     CODE BELOW IS OLD VERSION     ######
# Add lags ----------------------------------------------------------------

# Create object for each lag
multi.simp <- multivariate%>%
  dplyr::select(year, SPR, logR.dev, OISST, BT, SST, HW, GSI, CP, CC5, CAD)

df1 <-multi.simp 
df1$year <- df1$year+1
colnames(df1)[2:ncol(df1)]<-paste(colnames(multi.simp)[2:ncol(multi.simp)],"_1",sep="")
df1<-as.data.frame(df1)
df1<-df1%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
df1<-df1[,-c(1,2,3)] # remove extra columns


df2 <-multi.simp 
df2$year <- df2$year+2
colnames(df2)[2:ncol(df2)]<-paste(colnames(multi.simp)[2:ncol(multi.simp)],"_2",sep="")
df2<-as.data.frame(df2)
df2<-df2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
df2<-df2[,-c(1,2,3)] # remove extra columns

df3 <-multi.simp 
df3$year <- df3$year+3
colnames(df3)[2:ncol(df3)]<-paste(colnames(multi.simp)[2:ncol(multi.simp)],"_3",sep="")
df3<-as.data.frame(df3)
df3<-df3%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
df3<-df3[,-c(1,2,3)] # remove extra columns

# Join these back together
dflag<-cbind(multi.simp,df1,df2,df3)
regdat <- data.frame(dflag)
ncol(dflag)

# boosted regression trees
plot(regdat$year, regdat$SPR)

# Select optimal learning rate and number of trees
mod<-gbm.step(data=regdat,
              gbm.x=c(4:27), # Exclude 3 year lag for now
              gbm.y=2, # Spawner per recruit
              family="gaussian",
              tree.complexity=1, # (1 = no interactions)
              learning.rate=0.0001,
              bag.fraction=0.5) # fraction of data used in test set (.7)
# Results for tree complexity of 1
# bag fraction .7 and learning rate 0.001 --> 1050 trees for SPR, now try to reduce the lr
# bag fraction .7 and learning rate 0.005 --> 350 trees for SPR,
# But then tried again with same settings and it didn't work...
# and then again with .7, lr = 0.0001 --> 9150

# Tree complexity of 2
# bag fraction .7 and learning rate 0.001 --> 550 trees for SPR

summary(mod) # all the top variables have a two year lag (prob b/c this is spawner per recruit)

# Repeat for log deviations
mod2<-gbm.step(data=regdat,
               gbm.x=c(4:27), # Exclude 3 year lag for now
               gbm.y=3, # Log deviations
               family="gaussian",
               tree.complexity=1, # (1 = no interactions)
               learning.rate=0.01,
               bag.fraction=0.9)
# Results for tc 1
# bag fraction .7 and learning rate 0.001 --> 2500 trees 
# bag fraction .7 and learning rate 0.01 --> 300 trees 
summary(mod2) # top variables include hw2, sst, bt2, sst1, but ofc they change with each run
# Second run with larger learning rate includes sst (&1) hw2, cp2, bt2

# percent deviance explained ( (null dev - resid dev) / null dev ) * 100
null.dev<-mod.nolag$self.statistics$mean.null
resid.dev<-mod.nolag$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100

dev.expl

# relative influence
ri<-summary(mod2)

# plot
ggplot(data=ri,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()



# Remove OISST it's highly correlated with the others
head(regdat)
regdat2<-regdat%>%
  dplyr::select(!c("OISST", "OISST_1", "OISST_2"))
str(regdat2)

mod3<-gbm.step(data=regdat2,
               gbm.x=c(4:24), # Exclude 3 year lag for now
               gbm.y=3, # Log deviations
               family="gaussian",
               tree.complexity=1, # (1 = no interactions)
               learning.rate=0.001,
               bag.fraction=0.9)
ri3<-summary(mod3)
ggplot(data=ri3,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# SST(2&1), bt2, gsi2, cad2
mod.simp <- gbm.simplify(mod3)
