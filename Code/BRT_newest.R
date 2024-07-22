##      Title:            Herring Recruitment Boosted Regression Trees
##      Author:           Adelle Molina
##      Created:          4/28/23
##      Updated:          7/18/24
##      Notes:            This script has more recent model runs

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

# Load data ---------------------------------------------------------------

###############WHICH ONE


# BRT with WHAM output and Sarah zoo indices ------------------------------



# BRT for log recruitment deviations ---------
names(combined.dat)
ncol(combined.dat)

# old best
#(((BRT.ldev$self.statistics$mean.null)-(BRT.ldev$cv.statistics$deviance.mean))/(BRT.ldev$self.statistics$mean.null))*100 #46 % dev explained 
#summary(BRT.ldev) # Had, zooabun, sst_1, hw, had_1 --> only ppd and chla had zero relative abundance

# old best with 0 vars removed
#(((BRT.ldev_2$self.statistics$mean.null)-(BRT.ldev_2$cv.statistics$deviance.mean))/(BRT.ldev_2$self.statistics$mean.null))*100 #52 % dev explained 
#ri <- summary(BRT.ldev_2) # Had, zooabun, sst_1, hw, had_1 --> none with zero relative abundance

# Have run this a few times with diff sized start datasets 
# Previous had less lags, now run again with more
BRT.ldev_wlags <-gbm.step(data=combined.dat,
                     gbm.x=c(9:ncol(combined.dat)), 
                     gbm.y=4, # log recruitment deviations
                     family="gaussian", tree.complexity=1, # (1 = no interactions)
                     learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev_wlags$self.statistics$mean.null)-(BRT.ldev_wlags$cv.statistics$deviance.mean))/(BRT.ldev_wlags$self.statistics$mean.null))*100 #35 % dev explained (SAME WITH MORE LAGS)
summary(BRT.ldev_wlags)

# BRT for log recruitment deviations with 0 relative influence removed ---------
#combined.dat <- combined.dat%>%dplyr::select(-PPD, -ChlA, -GSI,  -NAOIndex, -NAOIndex_1, -NAOIndex_2, -NAOIndex_3)
combined.dat2 <- combined.dat%>%dplyr::select(-PPD, -PPD_1, -PPD_2,-PPD_3,
                                              -ChlA,  -ChlA_1, -ChlA_2,-ChlA_3, -GSI)

BRT.ldev_2_wlags <-gbm.step(data=combined.dat2,
                    gbm.x=c(9:ncol(combined.dat2)), 
                    gbm.y=4, # log recruitment deviations
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev_2_wlags$self.statistics$mean.null)-(BRT.ldev_2_wlags$cv.statistics$deviance.mean))/(BRT.ldev_2_wlags$self.statistics$mean.null))*100 #32 % dev explained (Weird that went down with less (all 0) vars in it)
ri <- summary(BRT.ldev_2_wlags)

# BRT for log recruitment deviations with correlated variables removed ---------
combined.dat3 <- combined.dat2%>%
  dplyr::select(-HW, -HW_1, -HW_2, -HW_3, -CP)

BRT.ldev_3_wlags <-gbm.step(data=combined.dat3,
                            gbm.x=c(9:ncol(combined.dat3)), 
                            gbm.y=4, # log recruitment deviations
                            family="gaussian", tree.complexity=1, # (1 = no interactions)
                            learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev_3_wlags$self.statistics$mean.null)-(BRT.ldev_3_wlags$cv.statistics$deviance.mean))/(BRT.ldev_3_wlags$self.statistics$mean.null))*100 #38 % dev explained (new is 36)
summary(BRT.ldev_3_wlags)
ri <- summary(BRT.ldev_3)
# Best model for now

# BRT for log recruitment deviations with more zero influences removed ---------
# drop ones that don't make sense 
# bt in year of recruit...theyre already off bottom, and 3 years before also nonsense
# why had I previously dropped cfin in year of recruitment (I think theyre still eating those at that stage but it had and still has 0 rel infl)
# Drop mackerel, its correlated and doesn't make sesne
combined.dat4 <- combined.dat3%>%
  dplyr::select(-BT, -BT_3, -MackSSB, -MackSSB_1)

BRT.ldev_4_wlags <-gbm.step(data=combined.dat4,
                      gbm.x=c(9:ncol(combined.dat4)), 
                      gbm.y=4, # log recruitment deviations
                      family="gaussian", tree.complexity=1, # (1 = no interactions)
                      learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev_4_wlags$self.statistics$mean.null)-(BRT.ldev_4_wlags$cv.statistics$deviance.mean))/(BRT.ldev_4_wlags$self.statistics$mean.null))*100 #34 % dev explained 
summary(BRT.ldev_4_wlags)

# Plot relative influence
ggplot(data=ri,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()+
  theme_bw()

#write_csv(ri,"Relativeinf_newest.csv")
# Export 
tiff("Figures/Relative Influence.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=ri,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()+
  theme_bw()
dev.off()


# Partial Dependence Plots ---------------------------------------------------------------
# For the model with highest deviance explained --> BRT.ldev_2
x<-c()
for(i in 1:(ncol(reduced.dat) -8)){
  pdp<-plot.gbm(BRT.ldev_3,i,return.grid=T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(ri$rel.inf[which(row.names(ri)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","logR.dev","var","ri")
  x<-rbind(x,pdp)
}

# order by relative influence
x$var<-factor(x$var,levels=row.names(ri))

# top 6
top6<-as.character(ri$var[1:6])
relinf<-round(ri$rel.inf[1:6],1)

x%>%
  filter(var==top6[1])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[1],"(",relinf[1],"%)",sep=""),y="Marginal effect on Recruit Devs")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p1

x%>%
  filter(var==top6[2])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[2],"(",relinf[2],"%)",sep=""),y="Marginal effect on logR.dev")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p2

x%>%
  filter(var==top6[3])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[3],"(",relinf[3],"%)",sep=""),y="Marginal effect on logR.dev")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p3

x%>%
  filter(var==top6[4])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[4],"(",relinf[4],"%)",sep=""),y="Marginal effect on logR.dev")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p4

x%>%
  filter(var==top6[5])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[5],"(",relinf[5],"%)",sep=""),y="Marginal effect on logR.dev")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p5

x%>%
  filter(var==top6[6])%>%
  ggplot(aes(x=val,y=logR.dev))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[6],"(",relinf[6],"%)",sep=""),y="Marginal effect on logR.dev")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p6

pdpz <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)

# Export 
tiff("Figures/Partial Dependence.tiff", width = 8, height = 8, units = 'in', res = 300)
pdpz <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()

# Other plots
par(mfrow=c(3,4))
gbm.plot(BRT.ldev_2_wlags, n.plots=12, write.title = F)

gbm.plot.fits(BRT.ldev_2_wlags)


# BRT for Rs with plots ---------
names(reduced.dat) # oops for these shoud prob take out ssb1...that may explain why its so much hieher
BRT.Rs_new <-gbm.step(data=reduced.dat,
                    gbm.x=c(9:ncol(reduced.dat)), 
                    gbm.y=8, # Recruitment Success
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
(((BRT.Rs_new$self.statistics$mean.null)-(BRT.Rs_new$cv.statistics$deviance.mean))/(BRT.Rs_new$self.statistics$mean.null))*100 #61 % dev explained (63 without ppd/chla) (53 with smaller dataset)
ri.Rs <- summary(BRT.Rs_new) # Had, jelly, zoo, ssb_1, mack_1 --> ppd and chla had zero relative abundance

tiff("Figures/Relative Influence (Rs).tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(data=ri.Rs,aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity")+
  labs(x="",y="relative influence")+
  coord_flip()+
  theme_bw()
dev.off()

# Partial dependence For the Rs model with highest deviance explained --> BRT.Rs
ncol(combined.dat)
x<-c()
for(i in 1:(ncol(combined.dat) -9)){
  pdp<-plot.gbm(BRT.Rs,i,return.grid=T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(ri.Rs$rel.inf[which(row.names(ri.Rs)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","Rs","var","ri")
  x<-rbind(x,pdp)
}

# order by relative influence
x$var<-factor(x$var,levels=row.names(ri.Rs))

# top 6
top6<-as.character(ri.Rs$var[1:6])
relinf<-round(ri.Rs$rel.inf[1:6],1)

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

# Save this
tiff("Partial Dependence (Rs).tiff", width = 8, height = 8, units = 'in', res = 300)
pdpz <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()
