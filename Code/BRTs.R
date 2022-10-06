##      Title:            Evaluate possible predictors
##      Script Purpose:   Run boosted regression trees
##      Author:           Adelle Molina
##      Created:          9/13/22
##      Updated:          9/13/22
##      Notes:            1. this is a mess...all diff tree stuff all over; model with only 2 year lags explained 30% of the deviance at first...now it's no higher than 6...the more I reduce variables
##      To Do:            2. Save figs --> using reduced lists so easier to see
##                        3. Run more trees --> use fuller dataset
##                        4. Plot tree results 
    
# Packages
library(dplyr)
library(ggplot2)
library(corrplot)
library(reshape2)
library(ggstatsplot)
library(gbm)
library(dismo)
library(gridExtra)

# Load data ---------------------------------------------------
# load possible variables list (from csv exported from other script where data is processed)
# this is a reduced list...prob should have just kept everything and thrown it all in
vars <- read.csv(file = 'variables.csv', header=T)
more <- read.csv(file = 'multivariatedata.csv', header=T)
ncol(vars)
vars <- vars[,-1]
names(vars)

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

#Calculate Spawning Success Rs (lnRt/ssb-1)
herring$Rs <- NA
for(i in 2:nrow(herring)){
  herring$Rs[i] <- log(herring$recruits[i]/herring$SSB[i-1])
}

# join 
dat <- merge(herring, vars,by = "year")
names(dat)
# Plot time series and correlations --------------------------------------------------------
# Temperature
temp.dat <- dat%>% 
  dplyr::select(year, BT.GOM, BT.GB, BT.MAB , BT.SS,
         SST.GOM, SST.GB, SST.MAB , SST.SS, 
         Surv.SST,  Surv.BT, 
         SPR, logR.dev, Rs) # woops doesn't work now b/c I took some of these out of vars
temps.melt <- melt(temp.dat, id.vars="year", variable.name="series")
ggplot(temps.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)

cor.temps <-ggstatsplot::ggcorrmat(
  data = temp.dat,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))

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

## # Biological 
bio.dat <- multivariate%>% 
  select(year,  Abund, Num.GB, Num.GOM, Num.MAB, ZooDens, CAD, CC5, Cope,SPR, logR.dev, Rs)
bio.melt <- melt(bio.dat, id.vars="year", variable.name="series")
ggplot(bio.melt, aes(year,value)) +
  geom_line()+ 
  theme_bw()+
  facet_wrap(~ series, scales = "free", ncol=4)
ggstatsplot::ggcorrmat(
  data = bio.dat,
  type = "nonparametric",
  colors = c("darkred", "white", "steelblue"))

# these are now in quarto....with slight updates



# Add lags ----------------------------------------------------------------
# remove vars that didn't seem as important --> where is best step to do this, prob not here wanna keep
vars <- vars%>% 
  dplyr::select(-ZooDens, -CAD, -CC5, -Cope, -CP, -GSI)

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

df3 <-dat 
df3$year <- df3$year+3
colnames(df3)[2:ncol(df3)]<-paste(colnames(dat)[2:ncol(dat)],"_3",sep="")
df3<-as.data.frame(df3)
df3<-df3%>%
  add_row(year=1982, .before=1)%>%
  add_row(year=1983, .after=1)%>%
  add_row(year=1984, .after=2)%>%
  dplyr::filter((year<2019)%>% replace_na(TRUE)) # add na to earlier years and chop off extra years
df3<-df3[,-c(1)] # remove extra columns

# Join lagged data with the herring data (select only 3 possible variables) and filter to years
full.dat <- herring%>% 
  dplyr::select(year, SPR, logR.dev, Rs)%>%
  dplyr::filter(year%in%c(seq(1982,2018)))

dflag<-cbind(herr,vars[,-1], lag1,lag2)
regdat <- data.frame(dflag)
names(regdat)


# Add lags to a smaller set of possible variables (no regional indices and no pp)
ts.vars2 <- ts.vars%>% 
  dplyr::select(-cal.GB, -cal.GOM, -Num.GB, - Num.GOM, -Num.MAB, -SST.GB, -SST.GOM, -SST.MAB, -SST.SS, -BT.GB, -BT.GOM, -BT.MAB, -BT.SS)
# Create a separate object for each lag (1-3 for now), shift, append lag to column name
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
dat.noreg <- merge(ts.recr, ts.vars2,by = "year")
dat.noreglag <- cbind(dat.noreg,lag1.1,lag2.1, lag3.1)
brtdat.noreg <- data.frame(dat.noreglag)
# Boosted Regression Trees ------------------------------------------------

# Run full trees (all possible vars with up to 3 year lags) for all three possible recruitment index
# 1. SPR (spawner per recruit vector...all I've presented thus far is for this one, but based on last chat w/ John perhaps this isn't the best)
# 2. Rs (this is the one to use)
# 3. Log recruit devs
#########   1. SPR with fuller version of data
ncol(brtdat) # how many columns
which(colnames(brtdat)=="SPR") # which column is SPR, the dependent variable
names(brtdat)
SPR.mod <-gbm.step(data=brtdat,
              gbm.x=c(1, 6, 9:92), # Select all columns after herring r indices, and also ssb and year?
              gbm.y=7, # Spawner per recruit
              family="gaussian", tree.complexity=1, # (1 = no interactions)
              learning.rate=0.01, bag.fraction=0.7)

null.dev<-SPR.mod$self.statistics$mean.null
resid.dev<-SPR.mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is lower than past models for same variable
SPR.rel <-summary(SPR.mod)
ggplot(data=SPR.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()


# Repeat with no regional indices
ncol(brtdat.noreg) # how many columns
which(colnames(brtdat.noreg)=="SPR") # which column is SPR, the dependent variable
names(brtdat.noreg)
SPR.mod.noreg <-gbm.step(data=brtdat.noreg,
                   gbm.x=c(1, 6, 9:40), # Select all columns after herring r indices, and also ssb and year?
                   gbm.y=4, # Spawner per recruit
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-SPR.mod.noreg$self.statistics$mean.null
resid.dev<-SPR.mod.noreg$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is again lower than past models for same variable, so taking away variables was worse
SPR.noreg.rel <-summary(SPR.mod.noreg)
ggplot(data=SPR.noreg.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

# now run a model for Rs
#########   2. Rs
ncol(brtdat) # how many columns
which(colnames(brtdat)=="Rs") # which column is Rs, the dependent variable
names(brtdat)
Rs.mod <-gbm.step(data=brtdat,
                   gbm.x=c(1, 6, 9:92), # Select all columns after herring r indices, and also ssb and year?
                   gbm.y=8, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod$self.statistics$mean.null
resid.dev<-Rs.mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is higher than for SPR (using both reduced and full possible vars)
Rs.rel <-summary(Rs.mod)
ggplot(data=Rs.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# again here the copepods in GB had an extremely (maybe too high relative influence)

# Repeat for Rs using slightly smaller dataset (all ones with 0 influence)

brtdat.sm <- brtdat%>% 
  dplyr::select(-PP, -PP_1, -PP_2, -PP_3, -Cope, -Cope_1, -Cope_2, -Cope_3, -Surv.SST, -Surv.SST_1, -Surv.SST_2, -Surv.SST_3, -Surv.BT, -Surv.BT_1, -Surv.BT_2, -Surv.BT_3)
ncol(brtdat.sm) # how many columns
names(brtdat.sm) # how many columns
Rs.mod1 <-gbm.step(data=brtdat.sm,
                   gbm.x=c(1, 6, 9:76), # Select all columns after herring r indices, and also ssb and year?
                   gbm.y=8, # Rs
                   family="gaussian", tree.complexity=1, # (1 = no interactions)
                   learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod1$self.statistics$mean.null
resid.dev<-Rs.mod1$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # this model has highest deviance explained thus far

Rs.rel1 <-summary(Rs.mod1)
ggplot(data=Rs.rel1,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

# Repeat for Rs using even smaller dataset
ncol(brtdat.noreg) # how many columns
which(colnames(brtdat.noreg)=="Rs") 
names(brtdat.noreg)
brtdat.noreg <- brtdat.noreg%>% 
  dplyr::select(-PP, -PP_1, -PP_2, -PP_3)
Rs.mod2 <-gbm.step(data=brtdat.noreg,
                  gbm.x=c(1, 6, 9:36), # Select all columns after herring r indices, and also ssb and year?
                  gbm.y=8, # Rs
                  family="gaussian", tree.complexity=1, # (1 = no interactions)
                  learning.rate=0.01, bag.fraction=0.7)

null.dev<-Rs.mod2$self.statistics$mean.null
resid.dev<-Rs.mod2$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is higher than for SPR (using both reduced and full possible vars)

Rs.rel2 <-summary(Rs.mod2)
ggplot(data=Rs.rel2,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()

#########   3. Recruitment deviations
devs.mod <-gbm.step(data=brtdat,
                  gbm.x=c(1, 6, 9:92), # Select all columns after herring r indices, and also ssb and year?
                  gbm.y=4, # logRdevs
                  family="gaussian", tree.complexity=1, # (1 = no interactions)
                  learning.rate=0.01, bag.fraction=0.7)

null.dev<-devs.mod$self.statistics$mean.null
resid.dev<-devs.mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl # dev explained is very low for devs (using full possible vars) and pretty low (using reduced vars), but interestingly one variable explains a huge amount of variance in full (cal GB), but perhaps spurious

# relative influence
devs.rel <-summary(devs.mod)

# plot
ggplot(data=devs.rel,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()




# old trees for only up to 2 year lags
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

summary(mod)
# percent deviance explained ( (null dev - resid dev) / null dev ) * 100
null.dev<-mod$self.statistics$mean.null
resid.dev<-mod$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl

# relative influence
ri<-summary(mod)

# plot
ggplot(data=ri,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()


# BRT Plots ---------------------------------------------------------------

plot.gbm(Rs.mod1,71)
ncol(brtdat.sm)
x<-c()
for(i in 1:(ncol(brtdat.sm) -6)){
  pdp<-plot.gbm(Rs.mod1,i,return.grid=T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(Rs.rel1$rel.inf[which(row.names(Rs.rel1)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","Rs","var","ri")
  x<-rbind(x,pdp)
}

# order by relative influence
x$var<-factor(x$var,levels=row.names(Rs.rel1))

# top 6
top6<-as.character(Rs.rel1$var[1:6])
relinf<-round(Rs.rel1$rel.inf[1:6],1)

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

# ok so for the variable I should be using I found a model that has high deviance explained and one variable explains tons of the variance



# BRT grouped variables ---------------------------------------------------
# only use nonlagged data
names(regdat)
# Fit with different settings
mod.nolag<-gbm.step(data=regdat,
              gbm.x=c(5:21), # Select all columns after herring r indices with  no lag
              gbm.y=2, # Spawner per recruit
              family="gaussian",
              tree.complexity=1, # (1 = no interactions)
              learning.rate=0.01,
              bag.fraction=0.7)# fraction of data used in test set (.7)

mod.nolag$n.trees
brd.mod$n.trees
null.dev<-mod.nolag$self.statistics$mean.null
resid.dev<-mod.nolag$cv.statistics$deviance.mean
dev.expl<-((null.dev-resid.dev)/null.dev)*100
dev.expl
ggplot(data=summary(mod.nolag),aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# GB zooplankton abundance $ MAB HW
help("plot.gbm")
plot(mod.nolag, )

names(regdat)
mod.nolag<-gbm.step(data=regdat,
                    gbm.x=c(5:21), # Select all columns after herring r indices with  no lag
                    gbm.y=2, # Spawner per recruit
                    family="gaussian",
                    tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01,
                    bag.fraction=0.7)# fraction of data used in test set (.7)

# messy disorganized
# old summaries of nt results
summary(mod.nolag)

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

mod

ggplot(data=ri3,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()
# SST(2&1), bt2, gsi2, cad2

# lr .1, bf .7, nt = 450
# lr .01, bf .7, nt = 1200
# lr .01, bf .7, nt = 9300
# lr .01, bf .7, nt = 1350, 1700
# bag fraction .7 and learning rate 0.0001 --> 7500 trees 
# bag fraction .7 and learning rate 0.01 --> 3500 trees 
# removed some predictors, bag fraction .7 and learning rate 0.001 --> 8500 trees 
# But predicted deviance is very low, doesn't seem right

mod.simp <- gbm.simplify(mod3)


mod<-gbm.step(data=regdat,
              gbm.x=c(2:25, 26:103), # Exclude 3 year lag for now
              gbm.y=26, # Spawner per recruit
              family="gaussian",
              tree.complexity=1, # (1 = no interactions)
              learning.rate=0.07,
              bag.fraction=0.5) # fraction of data used in test set (.7)