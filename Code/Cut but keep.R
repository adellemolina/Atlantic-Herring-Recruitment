
###########    OLD VERSION     ######
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

dev.expl1
# hmm ok yeah it keeps changing

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


basin.vars <- merge(physvar.basin, biovar.basin,
                    by = "year", 
                    all.x = TRUE, all.y = T)
reg.vars <- merge(physvar.reg, biovar.reg,
                  by = "year", 
                  all.x = TRUE, all.y = T)
write.csv(basin.vars, file = "Basin wide biophys variables.csv", row.names = F)
write.csv(reg.vars, file = "Regional biophys variables.csv", row.names = F)


# Old version with just basin scale variables
biovar.basin <- list(Chla, ppd, predators,
                     Cal.Stage.An, Cal.An,
                     Smcope.dens.sea, Smcope.dens.An,
                     smallcal, largecal,  ZooAbun.An.Wt)
biovar.basin <- biovar.basin %>% reduce(full_join, by='year')

biovar.reg <- list(PPD.reg, Smcope.dens.reg, zn.reg)
biovar.reg <- biovar.reg %>% reduce(full_join, by='year')

# Trim down to the year range we want --> not yet
#biovar.list <- biovar.list  %>% 
#dplyr::filter(between(year,1982,2020))
names(biovar.basin)


# Basin wide annual averages (Oisst, survey (sst, bt, sss, bs), indices)
#physvar.basin <- list(OISST.Annual, Surv.An.Wt, Fall, Spring,  HW.An, GSI.Annual, CP, WSW)
#physvar.basin <- physvar.basin %>% reduce(full_join, by='year')

# Regional (Oisst, survey (raw sst&bt), 
#physvar.reg <- list(OISST.reg, SurvBT.Reg, SurvSST.Reg, HW.reg)
#physvar.reg <- physvar.reg %>% reduce(full_join, by='year')
basin.vars <- read.csv(file = 'Data/Basin wide biophys variables.csv', header=T)
reg.vars <- read.csv(file = 'Data/Regional biophys variables.csv', header=T)

# join these
comb.vars <- merge(basin.vars, reg.vars,by = "year")

# Also bring in the weight at age ssb matrix
WAA = structure(c(
  1.300000e-002, 1.600000e-002, 1.600000e-002, 1.100000e-002, 1.100000e-002, 1.100000e-002, 1.400000e-002, 3.100000e-002, 1.100000e-002, 8.000000e-003, 1.500000e-002, 1.500000e-002, 1.300000e-002, 3.200000e-002, 1.500000e-002, 7.000000e-003, 1.500000e-002, 1.700000e-002, 2.400000e-002, 7.000000e-003, 6.000000e-003, 3.200000e-002, 1.000000e-002, 2.700000e-002, 2.700000e-002, 
  2.400000e-002, 2.400000e-002, 2.400000e-002, 2.400000e-002, 2.400000e-002, 2.700000e-002, 2.800000e-002, 1.400000e-002, 2.700000e-002, 2.600000e-002, 2.700000e-002, 3.300000e-002, 3.000000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 2.700000e-002, 3.200000e-002, 2.700000e-002, 2.700000e-002, 
  2.700000e-002, 2.700000e-002, 2.700000e-002, 2.200000e-002, 2.000000e-002, 3.000000e-002, 2.000000e-002, 
  3.800000e-002, 4.700000e-002, 4.300000e-002, 3.800000e-002, 4.100000e-002, 6.100000e-002, 6.800000e-002, 6.900000e-002, 5.100000e-002, 4.500000e-002, 5.500000e-002, 8.800000e-002, 4.500000e-002, 5.100000e-002, 7.300000e-002, 5.400000e-002, 3.900000e-002, 5.000000e-002, 6.900000e-002, 6.400000e-002, 4.700000e-002, 5.700000e-002, 6.800000e-002, 6.600000e-002, 6.800000e-002, 
  6.200000e-002, 6.300000e-002, 6.000000e-002, 4.700000e-002, 5.400000e-002, 5.100000e-002, 5.500000e-002, 5.900000e-002, 5.200000e-002, 6.000000e-002, 6.500000e-002, 5.600000e-002, 5.900000e-002, 5.900000e-002, 4.700000e-002, 5.400000e-002, 6.200000e-002, 6.400000e-002, 6.800000e-002, 5.700000e-002, 4.300000e-002, 4.900000e-002, 4.900000e-002, 6.100000e-002, 6.600000e-002, 
  5.700000e-002, 6.500000e-002, 5.800000e-002, 5.000000e-002, 4.700000e-002, 5.600000e-002, 4.200000e-002, 
  9.500000e-002, 9.600000e-002, 1.070000e-001, 6.900000e-002, 1.020000e-001, 1.260000e-001, 1.440000e-001, 1.540000e-001, 1.330000e-001, 1.240000e-001, 1.330000e-001, 1.320000e-001, 1.310000e-001, 1.190000e-001, 1.330000e-001, 1.040000e-001, 1.350000e-001, 1.390000e-001, 1.440000e-001, 1.400000e-001, 1.460000e-001, 1.160000e-001, 1.080000e-001, 1.170000e-001, 1.160000e-001, 
  1.060000e-001, 9.600000e-002, 1.020000e-001, 9.600000e-002, 8.600000e-002, 9.500000e-002, 8.800000e-002, 9.100000e-002, 9.200000e-002, 9.100000e-002, 1.110000e-001, 9.900000e-002, 9.900000e-002, 9.900000e-002, 9.100000e-002, 8.700000e-002, 8.900000e-002, 1.060000e-001, 1.060000e-001, 9.500000e-002, 8.900000e-002, 7.600000e-002, 9.000000e-002, 9.000000e-002, 1.060000e-001, 
  1.030000e-001, 8.000000e-002, 9.300000e-002, 9.700000e-002, 1.120000e-001, 1.090000e-001, 9.600000e-002, 
  1.130000e-001, 1.700000e-001, 1.720000e-001, 1.780000e-001, 1.340000e-001, 1.630000e-001, 1.700000e-001, 1.970000e-001, 1.700000e-001, 1.690000e-001, 1.880000e-001, 1.840000e-001, 1.750000e-001, 1.780000e-001, 1.870000e-001, 1.850000e-001, 1.920000e-001, 2.000000e-001, 2.140000e-001, 1.930000e-001, 2.080000e-001, 1.760000e-001, 1.590000e-001, 1.540000e-001, 1.720000e-001, 
  1.560000e-001, 1.420000e-001, 1.350000e-001, 1.370000e-001, 1.200000e-001, 1.230000e-001, 1.250000e-001, 1.240000e-001, 1.170000e-001, 1.230000e-001, 1.370000e-001, 1.340000e-001, 1.260000e-001, 1.370000e-001, 1.290000e-001, 1.310000e-001, 1.330000e-001, 1.400000e-001, 1.350000e-001, 1.380000e-001, 1.210000e-001, 1.100000e-001, 1.070000e-001, 1.240000e-001, 1.190000e-001, 
  1.360000e-001, 1.140000e-001, 1.210000e-001, 1.420000e-001, 1.560000e-001, 1.370000e-001, 1.240000e-001, 
  2.020000e-001, 2.240000e-001, 2.060000e-001, 2.230000e-001, 2.220000e-001, 1.910000e-001, 2.020000e-001, 2.350000e-001, 2.380000e-001, 1.960000e-001, 2.110000e-001, 2.100000e-001, 2.150000e-001, 2.080000e-001, 2.290000e-001, 2.500000e-001, 2.360000e-001, 2.400000e-001, 2.650000e-001, 2.390000e-001, 2.370000e-001, 2.270000e-001, 2.020000e-001, 1.920000e-001, 2.010000e-001, 
  1.890000e-001, 1.710000e-001, 1.640000e-001, 1.560000e-001, 1.380000e-001, 1.450000e-001, 1.500000e-001, 1.500000e-001, 1.380000e-001, 1.400000e-001, 1.560000e-001, 1.530000e-001, 1.430000e-001, 1.530000e-001, 1.550000e-001, 1.590000e-001, 1.630000e-001, 1.640000e-001, 1.620000e-001, 1.590000e-001, 1.460000e-001, 1.410000e-001, 1.230000e-001, 1.320000e-001, 1.550000e-001, 
  1.480000e-001, 1.510000e-001, 1.480000e-001, 1.600000e-001, 1.680000e-001, 1.550000e-001, 1.310000e-001, 
  2.650000e-001, 2.790000e-001, 2.270000e-001, 2.650000e-001, 2.650000e-001, 2.390000e-001, 2.480000e-001, 2.680000e-001, 2.950000e-001, 2.700000e-001, 2.480000e-001, 2.360000e-001, 2.430000e-001, 2.390000e-001, 2.530000e-001, 2.940000e-001, 3.010000e-001, 2.720000e-001, 2.970000e-001, 2.860000e-001, 2.680000e-001, 2.520000e-001, 2.380000e-001, 2.290000e-001, 2.340000e-001, 
  2.160000e-001, 2.050000e-001, 1.900000e-001, 1.800000e-001, 1.590000e-001, 1.620000e-001, 1.710000e-001, 1.740000e-001, 1.640000e-001, 1.570000e-001, 1.720000e-001, 1.660000e-001, 1.670000e-001, 1.710000e-001, 1.730000e-001, 1.830000e-001, 1.840000e-001, 1.840000e-001, 1.750000e-001, 1.790000e-001, 1.690000e-001, 1.680000e-001, 1.550000e-001, 1.440000e-001, 1.580000e-001, 
  1.690000e-001, 1.580000e-001, 1.690000e-001, 1.850000e-001, 1.820000e-001, 1.930000e-001, 1.770000e-001, 
  2.980000e-001, 3.020000e-001, 2.420000e-001, 2.980000e-001, 2.980000e-001, 2.760000e-001, 2.960000e-001, 2.890000e-001, 3.520000e-001, 2.900000e-001, 2.950000e-001, 2.780000e-001, 2.490000e-001, 2.520000e-001, 3.020000e-001, 3.190000e-001, 3.390000e-001, 3.280000e-001, 3.320000e-001, 3.130000e-001, 3.180000e-001, 2.710000e-001, 2.560000e-001, 2.640000e-001, 2.600000e-001, 
  2.330000e-001, 2.250000e-001, 2.200000e-001, 2.090000e-001, 1.800000e-001, 1.750000e-001, 1.880000e-001, 1.940000e-001, 1.870000e-001, 1.860000e-001, 1.980000e-001, 1.810000e-001, 1.830000e-001, 1.920000e-001, 1.940000e-001, 1.990000e-001, 2.030000e-001, 2.030000e-001, 1.880000e-001, 1.910000e-001, 1.830000e-001, 1.830000e-001, 1.880000e-001, 1.800000e-001, 1.650000e-001, 
  1.700000e-001, 1.710000e-001, 1.860000e-001, 1.960000e-001, 2.010000e-001, 2.080000e-001, 1.910000e-001, 
  3.550000e-001, 3.550000e-001, 3.710000e-001, 3.550000e-001, 3.110000e-001, 4.190000e-001, 3.530000e-001, 3.440000e-001, 3.790000e-001, 3.520000e-001, 3.620000e-001, 3.710000e-001, 3.420000e-001, 3.210000e-001, 3.890000e-001, 3.660000e-001, 3.790000e-001, 3.680000e-001, 4.130000e-001, 3.790000e-001, 2.690000e-001, 3.190000e-001, 3.150000e-001, 3.160000e-001, 3.290000e-001, 
  3.120000e-001, 3.060000e-001, 3.050000e-001, 3.090000e-001, 3.070000e-001, 2.750000e-001, 2.280000e-001, 2.220000e-001, 2.160000e-001, 2.050000e-001, 2.210000e-001, 2.010000e-001, 1.950000e-001, 1.980000e-001, 2.030000e-001, 1.980000e-001, 2.040000e-001, 2.070000e-001, 2.010000e-001, 2.090000e-001, 2.030000e-001, 1.980000e-001, 1.980000e-001, 1.990000e-001, 1.960000e-001, 
  1.950000e-001, 1.900000e-001, 1.850000e-001, 2.320000e-001, 2.430000e-001, 2.130000e-001, 2.130000e-001),
  .Dim = c(57,8),
  .Dimnames = list(
    c("1965", "1966", "1967", "1968", "1969", "1970", "1971", "1972", "1973", "1974", "1975", "1976", "1977", "1978", "1979", "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989",
      "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014",
      "2015", "2016", "2017", "2018", "2019", "2020", "2021"), 
    c("1", "2", "3", "4", "5", "6", "7", "8")))

# Figure out how to turn this into a time series....maybe average weight of reproductive adult
# I think maybe pull out ages 3+ and calculate average weight to age ratio --> nah the ssb already accounts for this, but what I need is the parameter that explains why


#Quickly compare EPU data from different packages (weights)
strat.area # from techdoc rdata
Weights # from survdat, computed by poststratifying
epu_sf # from ecodata, its the converted shapefile


# Trim down to the year range we want --> don't run yet
#physvar.list <- physvar.list %>%
#dplyr::filter(between(year,1982,2020))

# Combine and save
tiff("BT.tiff", width = 8, height = 8, units = 'in', res = 300)
ggarrange(ts.anombt.reg,ts.anombt.an,ts.rawbt.season, ts.rawbt.reg, ts.rawbt.an, ts.rawT.fall, ts.rawT.spring, ncol=2, nrow=4) 
dev.off()

# Raw survey sst & bt (seasonal, basin wide)
ts.rawT.fall <- ggplot(Fall, aes(x = year, y = FaSST))  + 
  theme_bw() +   
  geom_line(aes(col="red"), size=2)+
  geom_line(aes(y = FaBT, col="black"), na.rm = T, size = 2) +
  xlab("Year") +
  ylab("Fall T") +
  scale_colour_manual(name='Data',
                      values=c('black'="black", 'red'="red"),
                      labels=c('SST', 'BT'))

ts.rawT.spring <- ggplot(Spring, aes(x = year, y = SpSST))  + 
  theme_bw() +   
  geom_point()+
  geom_line(aes(col="red"), size=2)+
  geom_line(aes(y = SpBT, col="black"), na.rm = T, size = 2) +
  xlab("Year") +
  ylab("Spring T")+
  scale_colour_manual(name='Data',
                      values=c('black'="black", 'red'="red"),
                      labels=c('SST', 'BT'))

# Survey bottom temp Seasonal, Annual & Regional(raw)
ts.rawbt.season <- ggplot(Survey, aes(x = YEAR, y = BT, col = EPU, shape = SEASON))  + 
  theme_bw() +   
  geom_point(na.rm = F, size=2)+
  geom_line(na.rm = T, size = 1) +
  xlab("Year") +
  ylab("In Situ BT") # this is messy/noisy, don't use

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

# Summarize across both seasons 
Surv.An <- NEFSC.RAW%>%
  dplyr::group_by(YEAR, EPU) %>%
  dplyr::summarise(SST = mean(SURFTEMP, na.rm = TRUE),
                   SSS = mean(SURFSALIN, na.rm = TRUE),
                   BT = mean(BOTTEMP, na.rm = TRUE),
                   BS = mean(BOTSALIN, na.rm = TRUE))


# Summarize by season, year, and EPU
Surv <- NEFSC.RAW%>%
  dplyr::group_by(YEAR, SEASON, EPU) %>%
  dplyr::summarise(SST = mean(SURFTEMP, na.rm = TRUE),
                   SSS = mean(SURFSALIN, na.rm = TRUE),
                   BT = mean(BOTTEMP, na.rm = TRUE),
                   BS = mean(BOTSALIN, na.rm = TRUE))

###########################################################################################################################
# Cut from qBRT (old code)
###########################################################################################################################
# 2. Using a smaller set of possible variables (no regional and no PP) (NOT UPDATED) --------
# not right b/c years wrbng
# remove unimportant variables (primary prod, salinity) & some regional indices
ts.vars2 <- ts.vars%>% 
  dplyr::select(-cal.GB, -cal.GOM, -Num.GB, - Num.GOM, -Num.MAB, -SST.GB, -SST.GOM, -SST.MAB, -SST.SS, -BT.GB, -BT.GOM, -BT.MAB, -BT.SS)
lag1.1 <-ts.vars2 
lag1.1$year <- lag1.1$year+1
colnames(lag1.1)[2:ncol(lag1.1)]<-paste(colnames(ts.vars2)[2:ncol(ts.vars2)],"_1",sep="")
lag1.1<-as.data.frame(lag1.1)
lag1.1<-lag1.1%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag1.1<-lag1.1[,-c(1)] # remove extra year column

lag2.1 <-ts.vars2 
lag2.1$year <- lag2.1$year+2
colnames(lag2.1)[2:ncol(lag2.1)]<-paste(colnames(ts.vars2)[2:ncol(ts.vars2)],"_2",sep="")
lag2.1<-as.data.frame(lag2.1)
lag2.1<-lag2.1%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) 
lag2.1<-lag2.1[,-c(1)] 

lag3.1 <-ts.vars2 
lag3.1$year <- lag3.1$year+3
colnames(lag3.1)[2:ncol(lag3.1)]<-paste(colnames(ts.vars2)[2:ncol(ts.vars2)],"_3",sep="")
lag3.1<-as.data.frame(lag3.1)
lag3.1<-lag3.1%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) 
lag3.1<-lag3.1[,-c(1)] 

dat.noreg <- merge(herr, ts.vars2,by = "year")
dat.noreglag <- cbind(dat.noreg,lag1.1,lag2.1, lag3.1)
brtdat.noreg <- data.frame(dat.noreglag)

# 3. Add lags to the dataset without region (NOT UPDATED) -------------------------------
lag1.2 <-dat.all
lag1.2$year <- lag1.2$year+1
colnames(lag1.2)[2:ncol(lag1.2)]<-paste(colnames(dat.all)[2:ncol(dat.all)],"_1",sep="")
lag1.2<-as.data.frame(lag1.2)
lag1.2<-lag1.2%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag1.2<-lag1.2[,-c(1)] # remove extra year column

lag2.2 <-dat.all 
lag2.2$year <- lag2.2$year+2
colnames(lag2.2)[2:ncol(lag2.2)]<-paste(colnames(dat.all)[2:ncol(dat.all)],"_2",sep="")
lag2.2<-as.data.frame(lag2.2)
lag2.2<-lag2.2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) 
lag2.2<-lag2.2[,-c(1)] 

lag3.2 <-dat.all
lag3.2$year <- lag3.2$year+3
colnames(lag3.2)[2:ncol(lag3.2)]<-paste(colnames(dat.all)[2:ncol(dat.all)],"_3",sep="")
lag3.2<-as.data.frame(lag3.2)
lag3.2<-lag3.2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) 
lag3.2<-lag3.2[,-c(1)] 

# Combine 
short.dat <- merge(her, dat.all,by = "year")
lag.alldat <- cbind(short.dat,lag1.2,lag2.2, lag3.2)
lagged.alldat <- data.frame(lag.alldat)

#vars <- read.csv(file = 'variables.csv', header=T)
#ncol(vars)
#vars <- vars[,-1]

# Create one with just rs, devs, ssb, and year 
#her <- herring%>%
#dplyr::select(year, Rs, logR.dev, SSB)%>%
#dplyr::filter(year%in%c(seq(1982,2020)))

# join full herring dataset to the variables --> not yet
#dat <- merge(herring, vars,by = "year")


# Reduce to exclude the subregion variables
#dat.all <- dat%>%
#dplyr::select(year, OISST, SST, SSS, BT, BS, FaSST, FaBT, SpSST, SpBT, HW, GSI, CP, ZooAbun, CalAbun, ZooDens, HadSSB)

# 2. Using larger dataset -------------------------------------------------
#this vars variable doesn't exist anymore'
# Filter to years
herr <- herring%>% 
  dplyr::select(year, SPR, logR.dev, Rs, SSB)%>%
  dplyr::filter(year%in%c(seq(1982,2020)))

# remove highly correlated variables (haha nope, this doesn't work b/c includes vars that arent in object)
names(vars)
#ts.vars <- vars%>% 
#dplyr::select(-ZooDens,-OISST.GB,-OISST.GOM,-OISST.MAB,-OISST.SS,-SSS, -CAD, -CC5, -cal.SS, -CP, -GSI, -HW.GB, -HW.GOM, -HW.MAB, -PP.GB, -PP.GOM, -PP.MAB,-PP.GOM, -Abund)

# did this for now but fix this b/c diff list of pars
ts.vars <- vars

# Create a separate object for each lag (1-3 years), shift, append lag to column name
lag1 <-ts.vars 
lag1$year <- lag1$year+1
colnames(lag1)[2:ncol(lag1)]<-paste(colnames(ts.vars)[2:ncol(ts.vars)],"_1",sep="")
lag1<-as.data.frame(lag1)
lag1<-lag1%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag1<-lag1[,-c(1)] # remove extra year column

lag2 <-ts.vars 
lag2$year <- lag2$year+2
colnames(lag2)[2:ncol(lag2)]<-paste(colnames(ts.vars)[2:ncol(ts.vars)],"_2",sep="")
lag2<-as.data.frame(lag2)
lag2<-lag2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) 
lag2<-lag2[,-c(1)] 

lag3 <-ts.vars 
lag3$year <- lag3$year+3
colnames(lag3)[2:ncol(lag3)]<-paste(colnames(ts.vars)[2:ncol(ts.vars)],"_3",sep="")
lag3<-as.data.frame(lag3)
lag3<-lag3%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2021)%>% replace_na(TRUE)) 
lag3<-lag3[,-c(1)] 

# Combine 
dat2 <- merge(herr, ts.vars,by = "year")
lagdat<- cbind(dat2, lag1, lag2, lag3)
brtdat <- data.frame(lagdat)
###########################################################################################################################
# Cut from quarto do
temp.dat <- dat%>% 
  dplyr::select(year, c(names(dat)[9:29]),SPR, logR.dev, Rs, SSB)
temps.melt <- reshape2::melt(temp.dat, id.vars="year", variable.name="series")

tempts <- ggplot(temps.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)

tempcor <- ggstatsplot::ggcorrmat(
  data = temp.dat,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))

## # Biological 
bio.dat <- dat%>% 
  dplyr::select(year,  c(names(dat)[36:51]),SPR, logR.dev, Rs, SSB)
bio.melt <- reshape2::melt(bio.dat, id.vars="year", variable.name="series")
ggplot(bio.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)
ggstatsplot::ggcorrmat(
  data = bio.dat,
  type = "nonparametric",
  colors = c("darkred", "white", "steelblue"))


# old plotting and other code from correlations/time series and BRT scripts

# old versions of plots
ts.ZooReg.old <- ggplot(zoonum_reg, aes(x = year, y = Abund, col = EPU))  + # use log so scale is easily visible
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T)+
  xlab("Year") +
  ylab("Log Zooplankton Abundance (# of ind)")

zooabun2 <- ggplot(zoo.numb, aes(x = year, y = log(Abund), col = Var))  + # use log so scale is easily visible
  theme_bw() +  
  geom_point()+
  geom_line(na.rm = T, size = 2)+
  xlab("Year") +
  ylab("Zooplankton Abundance (log(Number)")+
  facet_wrap(~EPU)

# old versions of averaging with old data (only 2 stages)
# non regional regional
copepod <- copepod%>%
  mutate_if(is.character, as.factor)%>%
  dplyr::filter(EPU%in% c("GOM", "SS", "GB"))%>% # Select the spawning regions only
  dplyr::group_by(Time, Var) %>% # (6 total var's are stage and season combined)
  dplyr::summarise(Anmean = mean(Value, na.rm = TRUE))%>%
  as.data.frame()

#get regional annual averages for both stages
cope.reg <- calanus_stage%>%
  dplyr::filter(EPU%in% c("GOM", "SS", "GB"))%>% # Select the spawning regions only
  group_by(Time, EPU)%>%
  dplyr::summarise(CalAbun = mean(Value))%>%
  dplyr::rename(year=Time)

# Pivot and rename for merge
c.reg <- pivot_wider(cope.reg, names_from = EPU, values_from = CalAbun) 
c.reg <- c.reg %>%
  dplyr::rename(Cal.GB=GB,
                Cal.GOM=GOM,
                Cal.SS=SS)



################################        MESSY RANDOM CODE       ###############
bio.dat <- multivariate%>% 
  select(year,  Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,SPR, logR.dev, Rs)

# Environmental Indices
env.dat<- dat%>% 
  select(year, GSI, CP, HW, HW.GB, HW.GOM, HW.MAB,SPR, logR.dev, Rs)
env.melt <- melt(env.dat, id.vars="year", variable.name="series")
ggplot(env.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)
cor.env <-ggstatsplot::ggcorrmat(
  data = env.dat,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))
temp.dat <- basin%>% 
  dplyr::select(year, BT.GOM, BT.GB, BT.MAB , BT.SS,
                SST.GOM, SST.GB, SST.MAB , SST.SS, 
                Surv.SST,  Surv.BT) 
temps.melt <- melt(temp.dat, id.vars="year", variable.name="series")

tempts <- ggplot(tempvar.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)

# Old variable lists ------------------------------------------------------
# Join to herring data, don't use this b/c add lags first, then join
#multivariate <- merge(herring, vars, by = "year", all.x = TRUE, all.y = T)

# Trim off some years 
#multivariate <- multivariate %>% 
#filter(between(year,1982,2018))

# Save as a csv
#write.csv(multivariate, file = "combined.data.csv", row.names = F)

# Pull out primary production variables and salinity (need a shorter time series for those)
#multi.short <- multivariate%>% 
#select(year, OISST, BT, SST, Surv.SST, Surv.BT, 
#BT.GOM, BT.GB, BT.MAB , BT.SS, 
#SST.GOM, SST.GB, SST.MAB , SST.SS, 
#GSI, CP, HW, HW.GB, HW.GOM, HW.MAB,
#Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,
#recruits, SPR, logR.dev, Rs)

# Create a smaller version of variable list (no redundant survey data (anomalies vs raw, just use raw))
#var.list <- list(OISST.Annual, OISST.reg, Survey.Annual.weighted, SurvBT.Regional, SurvSST.Regional, 
#Fall, Spring,
#HW.Annual, HW.reg, GSI.Annual, CP, PP.Ann, PP.reg,
#zoonum_an, zn.reg, zooden_an,cope.adult, cope.c5, cope.Ann, c.reg, had) 
#var.list <- var.list %>% reduce(full_join, by='year')
#var.list <- var.list%>% 
#dplyr::select(-Surv.SSS, -Surv.BS) %>% 
#dplyr::filter(between(year,1982,2018))
#write.csv(vars, file = "Data/variables.csv", row.names = F)

#had  <- read.csv(file.choose(), header = T, blank.lines.skip = T)


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