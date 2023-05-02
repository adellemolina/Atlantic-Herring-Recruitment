##      Title:            Evaluate possible predictors
##      Script Purpose:   Run boosted regression trees
##      Author:           Adelle Molina
##      Created:          9/13/22
##      Updated:          4/28/23
##      Notes:            Old version where I add lags here, Clean up, organize, run more trees

# Packages ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(reshape2)
library(ggstatsplot)
library(gbm)
library(dismo)
library(gridExtra)

# Load data ---------------------------------------------------
# load variables lists (from csv exported from "other"Corr & ts" script where raw data is processed)

basin.vars <- read.csv(file = 'Basin wide biophys variables.csv', header=T)
reg.vars <- read.csv(file = 'Regional biophys variables.csv', header=T)
# create herring data
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

# Plot time series and correlations --------------------------------------------------------

# Basin scale temp, sal, and environnmental indices
basin.env <- basin.vars%>% 
  dplyr::select(year, WSW, CP, GSI, HW, SpBT, SpSST, FaBT, FaSST, BS, BT, SSS, SST, OISST) 
basin.melt <- melt(basin.env, id.vars="year", variable.name="series")

ggplot(basin.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)

cor.env <-ggstatsplot::ggcorrmat(
  data = basin.env,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))

## # Basin scale biological variables
bio.dat <- basin.vars%>% 
  dplyr::select(year,ChlA, PPD, HadSSB_GoM, MackSSB, CalAD, CalC5, CalC4, CalC3, Cal.AnAbun, Smcope.Falldens, Smcope.Springdens, Smcope.AnDens, Sm.cal.abun, Lg.cal.abun, ZooAbun)
bio.melt <- melt(bio.dat, id.vars="year", variable.name="series")
ggplot(bio.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)
ggstatsplot::ggcorrmat(
  data = bio.dat,
  type = "nonparametric",
  colors = c("darkred", "white", "steelblue"))

# LAGS & BRT's ----------------------------------------------------------------

# 1. Add lags to full basin scale variable list -------------------------------------------------
lag_1 <-basin.vars %>%
  dplyr::filter(year%in%c(seq(1983,2020)))
lag_1$year <- lag_1$year+1
colnames(lag_1)[2:ncol(lag_1)]<-paste(colnames(basin.vars)[2:ncol(basin.vars)],"_1",sep="")
lag_1<-as.data.frame(lag_1)
lag_1<-lag_1%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag_1<-lag_1[,-c(1)] # remove extra year column

lag_2 <-basin.vars %>%
  dplyr::filter(year%in%c(seq(1983,2020)))
lag_2$year <- lag_2$year+2
colnames(lag_2)[2:ncol(lag_2)]<-paste(colnames(basin.vars)[2:ncol(basin.vars)],"_2",sep="")
lag_2<-as.data.frame(lag_2)
lag_2<-lag_2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) 
lag_2<-lag_2[,-c(1)] 

lag_3 <-basin.vars %>%
  dplyr::filter(year%in%c(seq(1983,2020)))
lag_3$year <- lag_3$year+3
colnames(lag_3)[2:ncol(lag_3)]<-paste(colnames(basin.vars)[2:ncol(basin.vars)],"_3",sep="")
lag_3<-as.data.frame(lag_3)
lag_3<-lag_3%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2020)%>% replace_na(TRUE)) 
lag_3<-lag_3[,-c(1)] 

# Combine 
basin.dat <- merge(herring, basin.vars,by = "year")%>%
  dplyr::filter(year%in%c(seq(1983,2019)))
lag.basindat<- cbind(basin.dat, lag_1, lag_2, lag_3)
lagbasin <- data.frame(lag.basindat)


# 2. Small tree without regional variables or lags ---------
names(basin.dat)
simp.ldev <-gbm.step(data=basin.dat,
                 gbm.x=c(6,9:36), 
                 gbm.y=4, # log recruitment deviations
                 family="gaussian", tree.complexity=1, # (1 = no interactions)
                 learning.rate=0.01, bag.fraction=0.7)
(((simp.ldev$self.statistics$mean.null)-(simp.ldev$cv.statistics$deviance.mean))/(simp.ldev$self.statistics$mean.null))*100
#31 % dev explained 
summary(simp.ldev) # Had, zooabun, hw, wsw, mack 

# not important --> ppd, chla
# other low importance and correlated with important vars (seasonal sst & bt, seasonal small copepod density)

# 3. Smaller tree (remove unimportant/correlated variables) ---------
basin.reduced <- basin.dat%>%
  dplyr::select(-PPD, -ChlA, -SpBT, -FaSST, -FaBT, -SpSST,-Smcope.Springdens, -Smcope.Falldens)%>%
  dplyr::filter(year%in%c(seq(1983,2018)))
names(basin.reduced)
simp2.ldev <-gbm.step(data=basin.reduced,
                     gbm.x=c(6,9:28), 
                     gbm.y=4, # log recruitment deviations
                     family="gaussian", tree.complexity=1, # (1 = no interactions)
                     learning.rate=0.01, bag.fraction=0.7)
(((simp2.ldev$self.statistics$mean.null)-(simp2.ldev$cv.statistics$deviance.mean))/(simp2.ldev$self.statistics$mean.null))*100
#34 % dev explained --> less vars but went up

summary(simp2.ldev) # Had, wsw, smcope dens, zooabun, cal ad, calc4
# still too many confounding food variables
# remove small and large cal abun corr with zoo abun
# each stage correlated to annual, try removing stages first, vs removing the annual b/c density corr with adults
basin.reduced2 <- basin.reduced%>%
  dplyr::select(-SSS, -Lg.cal.abun, -Sm.cal.abun, -CalAD, -CalC3, -CalC4,-CalC5)
simp3.ldev <-gbm.step(data=basin.reduced2,
                      gbm.x=c(6,9:21), 
                      gbm.y=4, # log recruitment deviations
                      family="gaussian", tree.complexity=1, # (1 = no interactions)
                      learning.rate=0.01, bag.fraction=0.7)
(((simp3.ldev$self.statistics$mean.null)-(simp3.ldev$cv.statistics$deviance.mean))/(simp3.ldev$self.statistics$mean.null))*100
#35 % dev explained --> less vars but went up, goes down if remove more env
summary(simp3.ldev) # had, zooabun, hw, mack, cal.anabun

basin.reduced3 <- basin.reduced%>%
  dplyr::select(-SSS, -Lg.cal.abun, -Sm.cal.abun,  - SST)
simp4.ldev <-gbm.step(data=basin.reduced3,
                      gbm.x=c(6,9:24), 
                      gbm.y=4, # log recruitment deviations
                      family="gaussian", tree.complexity=1, # (1 = no interactions)
                      learning.rate=0.01, bag.fraction=0.7)
(((simp4.ldev$self.statistics$mean.null)-(simp4.ldev$cv.statistics$deviance.mean))/(simp4.ldev$self.statistics$mean.null))*100
#41 % dev explained --> hmm more vars but went up, suggests keeping the stage jawns, went down once removed GSI and CP
summary(simp4.ldev) # had, zooabun, hw, mack, calc4, no zeroes 
# add lags to this

# 4. Add lags to reduced dataset -------------------------------------------------

basin.reduced3.1 <-basin.reduced3 %>%
  dplyr::select(-recruits, -R.no.devs, -logR.dev, -SR.std.resid, -SPR, - Rs) , , , , , , , , ,

lag_1.1 <-basin.reduced3.1 %>%
  dplyr::filter(year%in%c(seq(1983,2018)))
lag_1.1$year <- lag_1.1$year+1
colnames(lag_1.1)[2:ncol(lag_1.1)]<-paste(colnames(basin.reduced3.1)[2:ncol(basin.reduced3.1)],"_1",sep="")
lag_1.1<-as.data.frame(lag_1.1)
lag_1.1<-lag_1.1%>%
  add_row(year=1983, .before=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag_1.1<-lag_1.1[,-c(1)] # remove extra year column

lag_2.1 <-basin.reduced3.1 %>%
  dplyr::filter(year%in%c(seq(1983,2018)))
lag_2.1$year <- lag_2.1$year+2
colnames(lag_2.1)[2:ncol(lag_2.1)]<-paste(colnames(basin.reduced3.1)[2:ncol(basin.reduced3.1)],"_2",sep="")
lag_2.1<-as.data.frame(lag_2.1)
lag_2.1<-lag_2.1%>%
  add_row(year=1983, .before=1)%>%
  add_row(year=1984, .after=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) 
lag_2.1<-lag_2.1[,-c(1)] 

lag_3.1 <-basin.reduced3.1 %>%
  dplyr::filter(year%in%c(seq(1983,2018)))
lag_3.1$year <- lag_3.1$year+3
colnames(lag_3.1)[2:ncol(lag_3.1)]<-paste(colnames(basin.reduced3.1)[2:ncol(basin.reduced3.1)],"_3",sep="")
lag_3.1<-as.data.frame(lag_3.1)
lag_3.1<-lag_3.1%>%
  add_row(year=1983, .before=1)%>%
  add_row(year=1984, .after=1)%>%
  add_row(year=1985, .after=2)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) 
lag_3.1<-lag_3.1[,-c(1)] 

# Combine 
herr <- herring%>% 
  dplyr::select(year, logR.dev)
basin.comb <- merge(herr, basin.reduced3.1,by = "year")%>%
  dplyr::filter(year%in%c(seq(1983,2018)))
lag.basindat2<- cbind(basin.comb, lag_1.1, lag_2.1, lag_3.1)
lagbasin2 <- data.frame(lag.basindat2)


# 5. Large tree with 3 year lags and only important vars  ---------
names(lagbasin2)
lag.ldev <-gbm.step(data=lagbasin2,
                      gbm.x=c(3:70), 
                      gbm.y=2, # log recruitment deviations
                      family="gaussian", tree.complexity=1, # (1 = no interactions)
                      learning.rate=0.01, bag.fraction=0.7)
(((lag.ldev$self.statistics$mean.null)-(lag.ldev$cv.statistics$deviance.mean))/(lag.ldev$self.statistics$mean.null))*100
#27 % dev explained --> less than in model without lags

summary(lag.ldev) # lots of zeroes, 3 year lag on most food variables doesn't matter, 


# Old trees (based on outdated datasets) ----------------------------------

# NO LAGS
basic <-gbm.step(data=dat.all,
                   gbm.x=c(4:20), 
                   gbm.y=2, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)
basic.dev <- (((basic$self.statistics$mean.null)-(basic$cv.statistics$deviance.mean))/(basic$self.statistics$mean.null))*100
basic.dev # 46%
basic.sum <-summary(basic)
basic.sum # top vars are had, zooabun, ssb, bs, hw and all have high relative influence

 # Repeat for logdevs
basic2 <-gbm.step(data=dat.all,
                 gbm.x=c(4:20), 
                 gbm.y=3, # logrdev
                 family="gaussian", tree.complexity=1, # (1 = no interactions)
                 learning.rate=0.01, bag.fraction=0.7)
basic2.dev <- (((basic2$self.statistics$mean.null)-(basic2$cv.statistics$deviance.mean))/(basic2$self.statistics$mean.null))*100
basic2.dev # 36%
basic2.sum <-summary(basic2)
basic2.sum # top 5 (had, zooabun, zooden, hw, fabt)

# NO REGIONALS
ncol(lagged.alldat)
lag.all <-gbm.step(data=lagged.alldat,
                   gbm.x=c(4:68), 
                   gbm.y=2, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)
lag.all.dev <- (((lag.all$self.statistics$mean.null)-(lag.all$cv.statistics$deviance.mean))/(lag.all$self.statistics$mean.null))*100
lag.all.dev # was 31%, now 27 with fixed vars/lags, still much lower than in the models above without lags
summary(lag.all) # top 5 (had, sss, bs, sss1, spsst3)

# Repeat for deviations
lag.all2 <-gbm.step(data=lagged.alldat,
                    gbm.x=c(4:68), 
                    gbm.y=3, # logRdev
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
lag.all2.dev <- (((lag.all2$self.statistics$mean.null)-(lag.all2$cv.statistics$deviance.mean))/(lag.all2$self.statistics$mean.null))*100
lag.all2.dev # 35%,  higher than same model for Rs
summary(lag.all2) # top 5 (had, hw2, spsst3, oisst1, hw)

# 2. No region with lags --------------------------------------------------
simplag.ldev <-gbm.step(data=lagbasin,
                    gbm.x=c(6,9:116), 
                    gbm.y=4, # logRdev
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
(((simplag.ldev$self.statistics$mean.null)-(simplag.ldev$cv.statistics$deviance.mean))/(simplag.ldev$self.statistics$mean.null))*100
#26 % dev explained
summary(simplag.ldev) # Had, zoo, hw, calc4, mack (chla & ppd 0)
# 2. Remove seasonal temps -----------------------------------------------
dat.all.sm <- dat.all%>%
  dplyr::select(-FaSST, -FaBT, -SpSST, -SpBT)
# this won't work now...need to add back herring vars

basic1 <-gbm.step(data=dat.all.sm,
                  gbm.x=c(4:16), 
                  gbm.y=2, # Rs
                  family="gaussian", tree.complexity=1, # (1 = no interactions)
                  learning.rate=0.01, bag.fraction=0.7)
basic1.dev <- (((basic1$self.statistics$mean.null)-(basic1$cv.statistics$deviance.mean))/(basic1$self.statistics$mean.null))*100
basic1.dev # 57 when year was included, down to 46 without year aka something represented in year is missing but same as model above
basic.sum1 <- summary(basic1) # top 5 (had, zooabun, ssb, hw, zoodens)


# diff model for deviations based on larger dataset
# 1. No regionals
names(brtdat)
devsmod1 <-gbm.step(data=brtdat,
                 gbm.x=c(5,6, 11,13, 23, 27:28,33, 38:40,44,48),
                 gbm.y=3, # log R devs
                 family="gaussian", tree.complexity=1, # (1 = no interactions)
                 learning.rate=0.01, bag.fraction=0.7)
devs.dev <- (((devsmod1$self.statistics$mean.null)-(devsmod1$cv.statistics$deviance.mean))/(devsmod1$self.statistics$mean.null))*100
devs.dev
devs.sum <-summary(devsmod1)
plot(devsmod1,13)
# Similar results as for Rs but some more sensible ones...hmm need to decide on maybe one food variable
# Early BRT models using Rs ----------------------------------------------------------------------

# 1. Using a full dataset
ncol(brtdat) # how many columns
which(colnames(brtdat)=="Rs") # which column is Rs, the dependent variable
names(brtdat)
Rs.mod <-gbm.step(data=brtdat,
                   gbm.x=c(1, 5, 6:85), # Select all columns after herring r indices, ssb, year
                   gbm.y=4, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod$self.statistics$mean.null
resid.dev<-Rs.mod$cv.statistics$deviance.mean
dev_a <-((null.dev-resid.dev)/null.dev)*100
dev_a # dev explained is higher than for SPR (using both reduced and full possible vars) 

Rs.rel <-summary(Rs.mod)
ggplot(data=Rs.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# still lots of unimportant variables, zoo and copepods in GB had an extremely high relative influence (maybe too much)

# 2. Repeat w slightly smaller dataset (remove variables with 0 or near 0 influence)
brtdat.sm <- brtdat%>% 
  dplyr::select(-PP, -PP_1, -PP_2, -PP_3, -Cope, -Cope_1, -Cope_2, -Cope_3, -Surv.SST, -Surv.SST_1, -Surv.SST_2, -Surv.SST_3, -Surv.BT, -Surv.BT_1, -Surv.BT_2, -Surv.BT_3)
ncol(brtdat.sm) # how many columns
names(brtdat.sm) # 
Rs.mod1 <-gbm.step(data=brtdat.sm,
                   gbm.x=c(1, 5, 6:77), 
                   gbm.y=4, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)
# This was the "best" model
null.dev<-Rs.mod1$self.statistics$mean.null
resid.dev<-Rs.mod1$cv.statistics$deviance.mean
dev.expl1<-((null.dev-resid.dev)/null.dev)*100
dev.expl1 

# 3. Similar to 2, but was using a slightly different version of brtdat.sm (about 8 less parameters)
Rs.mod1.1 <-gbm.step(data=brtdat.sm,
                   gbm.x=c(1, 5, 6:69),
                   gbm.y=4, # Rs
                   family="gaussian", tree.complexity=2, 
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod1.1$self.statistics$mean.null
resid.dev<-Rs.mod1.1$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # this model has highest deviance explained thus far

Rs.rel1 <-summary(Rs.mod1.1)
ggplot(data=Rs.rel1,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# same pattern as model above
# Can probably remove even more 0 influence sets (like SST.MAB) and in other cases can remove lags (like sst.ss, and surv bs)

# 4. Repeat for Rs using even smaller dataset
brtdat.noreg2 <- brtdat.noreg%>% 
  dplyr::select(-PP, -PP_1, -PP_2, -PP_3)
which(colnames(brtdat.noreg)=="Rs") 
Rs.mod2 <-gbm.step(data=brtdat.noreg2,
                  gbm.x=c(1, 5, 6:29), # Select all columns after herring r indices
                  gbm.y=4, # Rs
                  family="gaussian", tree.complexity=1, # (1 = no interactions)
                  learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod2$self.statistics$mean.null
resid.dev<-Rs.mod2$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is much lower (prob b/c removed regional food indices)

Rs.rel2 <-summary(Rs.mod2)
ggplot(data=Rs.rel2,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

# Early BRT models using SPR or Deviations ---------------------------------------------------------------------
# with fuller version of data
ncol(brtdat) # how many columns
which(colnames(brtdat)=="SPR") # which column is SPR, the dependent variable
SPR.mod <-gbm.step(data=brtdat,
                   gbm.x=c(1, 5, 6:85), # Select all columns after herring r indices, and also ssb and year?
                   gbm.y=2, # Spawner per recruit
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-SPR.mod$self.statistics$mean.null
resid.dev<-SPR.mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is lower than past models for same variable (14%) w 1750 trees
SPR.rel <-summary(SPR.mod)
ggplot(data=SPR.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# lots of unimportant variables, highest relative influence is 10% for GOM zoo abun

# Repeat with no regional indices
ncol(brtdat.noreg) # how many columns
which(colnames(brtdat.noreg)=="SPR") # which column is SPR, the dependent variable
names(brtdat.noreg)
SPR.mod.noreg <-gbm.step(data=brtdat.noreg,
                         gbm.x=c(1, 5, 6:33), # Select all columns after herring r indices, ssb and year
                         gbm.y=2, # Spawner per recruit
                         family="gaussian", tree.complexity=1, # (1 = no interactions)
                         learning.rate=0.001, bag.fraction=0.7)

null.dev<-SPR.mod.noreg$self.statistics$mean.null
resid.dev<-SPR.mod.noreg$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is again lower (actually negative, is that a problem?), so taking away variables was worse, but oisst has high relative influence
SPR.noreg.rel <-summary(SPR.mod.noreg)
ggplot(data=SPR.noreg.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

# Run a model with full data for the recruitment deviations
devs.mod <-gbm.step(data=brtdat,
                    gbm.x=c(1, 5, 6:85),
                    gbm.y=3, # logRdevs
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)

null.dev<-devs.mod$self.statistics$mean.null
resid.dev<-devs.mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is very low for devs (using full possible vars) and pretty low (using reduced vars)
devs.rel <-summary(devs.mod)
ggplot(data=devs.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()


# Partial Dependence Plots ---------------------------------------------------------------

# For the model with highest deviance explained 
names(dat.all.sm)
x<-c()
for(i in 1:(ncol(dat.all.sm) -1)){
  pdp<-plot.gbm(basic1,i,return.grid=T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(basic.sum1$rel.inf[which(row.names(basic.sum1)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","Rs","var","ri")
  x<-rbind(x,pdp)
}

# order by relative influence
x$var<-factor(x$var,levels=row.names(basic.sum1))

# top 6
top6<-as.character(basic.sum1$var[1:6])
relinf<-round(basic.sum1$rel.inf[1:6],1)

x%>%
  filter(var==top6[1])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[1],"(",relinf[1],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p1

x%>%
  filter(var==top6[2])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[2],"(",relinf[2],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p2

x%>%
  filter(var==top6[3])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[3],"(",relinf[3],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p3

x%>%
  filter(var==top6[4])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[4],"(",relinf[4],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p4

x%>%
  filter(var==top6[5])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[5],"(",relinf[5],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p5

x%>%
  filter(var==top6[6])%>%
  ggplot(aes(x=val,y=Rs))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[6],"(",relinf[6],"%)",sep=""),y="Marginal effect on Rs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p6

pdpz <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)

# More messy BRTs ---------------------------------------------------
# not organized below her

# Check to see if the reduced dataset still has any sig correlations
brtdat.sm2.nolag <- brtdat.sm2[,1:22]
cors <- ggstatsplot::ggcorrmat(
  data = brtdat.sm2.nolag,
  type = "nonparametric",
  colors = c("darkred", "white", "steelblue"))
cors
# Fit with different settings
mod.nolag<-gbm.step(data=brtdat.sm,
                    gbm.x=c(1,5:23), # Select all columns after herring r indices with  no lag
                    gbm.y=4, # Spawner per recruit
                    family="gaussian",
                    tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01,
                    bag.fraction=0.7)# fraction of data used in test set (.7)
summary(mod.nolag)
mod.nolag$n.trees
brd.mod$n.trees
null.dev<-mod.nolag$self.statistics$mean.null
resid.dev<-mod.nolag$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl
# this model with no lags had almost just as high of deviance explained (40, oope nope that's for a spr model), nope it's a bit lower at 36
ggplot(data=summary(mod.nolag),aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

mod.nolag<-gbm.step(data=regdat,
                    gbm.x=c(5:21), # Select all columns after herring r indices with  no lag
                    gbm.y=2, # Spawner per recruit
                    family="gaussian",
                    tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01,
                    bag.fraction=0.7)# fraction of data used in test set (.7)
summary(mod.nolag)

# Original tree exploration (settings and play around) -----------------------------------------------

# Create object for each lag using just the physical variables 
lag1 <-vars 
lag1$year <- lag1$year+1
colnames(lag1)[2:ncol(lag1)]<-paste(colnames(vars)[2:ncol(vars)],"_1",sep="")
lag1<-as.data.frame(lag1)
lag1<-lag1%>%
  add_row(year=1982, .before=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag1<-lag1[,-c(1)] # remove extra year column

lag2 <-vars 
lag2$year <- lag2$year+2
colnames(lag2)[2:ncol(lag2)]<-paste(colnames(vars)[2:ncol(vars)],"_2",sep="")
lag2<-as.data.frame(lag2)
lag2<-lag2%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
lag2<-lag2[,-c(1)] # remove extra columns

dflag<-cbind(herr,vars[,-1], lag1,lag2)
regdat <- data.frame(dflag)

# Select optimal learning rate and number of trees --> with reduced possible pars
ncol(regdat) # how many columns
which(colnames(regdat)=="SPR") # which column is SPR, the dependent variable
mod<-gbm.step(data=regdat,
              gbm.x=c(5:55), # Select all columns after herring r indices
              gbm.y=2, # Spawner per recruit
              family="gaussian",
              tree.complexity=1, # (1 = no interactions)
              learning.rate=0.01,
              bag.fraction=0.7)# fraction of data used in test set (.7)
#lr=.1, bf = .7, nt = 2050

# percent deviance explained ( (null dev - resid dev) / null dev ) * 100
null.dev<-mod$self.statistics$mean.null
resid.dev<-mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl
ri<-summary(mod)
ggplot(data=ri,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

mod2<-gbm.step(data=regdat,
               gbm.x=c(5:38), # Exclude 3 year lag for now
               gbm.y=2, # SPR
               family="gaussian",
               tree.complexity=1, # (1 = no interactions)
               learning.rate=0.01,
               bag.fraction=0.5,
               max.trees=20000,
               step.size=50)

mod3<-gbm.step(data=regdat2,
               gbm.x=c(4:24), # Exclude 3 year lag for now
               gbm.y=3, # Log deviations
               family="gaussian",
               tree.complexity=1, # (1 = no interactions)
               learning.rate=0.001,
               bag.fraction=0.9,
               max.trees=20000,
               step.size=50,
               n.trees=50)
ri3<-summary(mod2)
ggplot(data=ri3,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# SST(2&1), bt2, gsi2, cad2

mod.simp <- gbm.simplify(mod3)
mod<-gbm.step(data=regdat,
              gbm.x=c(2:25, 26:103), # Exclude 3 year lag for now
              gbm.y=26, # Spawner per recruit
              family="gaussian",
              tree.complexity=1, # (1 = no interactions)
              learning.rate=0.07,
              bag.fraction=0.5) # fraction of data used in test set (.7)