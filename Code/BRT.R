##      Title:            Herring Recruitment Boosted Regression Trees
##      Author:           Adelle Molina
##      Created:          4/28/23
##      Updated:          5/2/23
##      Notes:             

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

# Combine the data called herring and variables into one dataframe
combined.dat <-merge(herring, variables,
                     by = "year", 
                     all.x = TRUE, all.y = T)
str(combined.dat) # lots of NA's still

# Export to reload easily into quarto
write.csv(combined.dat, file = "Data/Combined.Data.csv", row.names = F)

# chop it down
combined.dat <- combined.dat %>% 
  filter(between(year,1982,2019))

# BRT for log recruitment deviations ---------
names(combined.dat)
BRT.ldev <-gbm.step(data=combined.dat,
                     gbm.x=c(9:29), 
                     gbm.y=4, # log recruitment deviations
                     family="gaussian", tree.complexity=1, # (1 = no interactions)
                     learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev$self.statistics$mean.null)-(BRT.ldev$cv.statistics$deviance.mean))/(BRT.ldev$self.statistics$mean.null))*100 #46 % dev explained 
summary(BRT.ldev) # Had, zooabun, sst_1, hw, had_1 --> only ppd and chla had zero relative abundance

# BRT for log recruitment deviations with ppd and chla removed ---------
names(combined.dat)

combined.dat <- combined.dat%>%
  dplyr::select(-PPD, -ChlA)

BRT.ldev_2 <-gbm.step(data=combined.dat,
                    gbm.x=c(9:27), 
                    gbm.y=4, # log recruitment deviations
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
(((BRT.ldev_2$self.statistics$mean.null)-(BRT.ldev_2$cv.statistics$deviance.mean))/(BRT.ldev_2$self.statistics$mean.null))*100 #52 % dev explained 
ri <- summary(BRT.ldev_2) # Had, zooabun, sst_1, hw, had_1 --> none with zero relative abundance

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
ncol(combined.dat)
x<-c()
for(i in 1:(ncol(combined.dat) -9)){
  pdp<-plot.gbm(BRT.ldev_2,i,return.grid=T)
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
tiff("Partial Dependence.tiff", width = 8, height = 8, units = 'in', res = 300)
pdpz <- grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()

# BRT for Rs with plots ---------
names(combined.dat)
BRT.Rs <-gbm.step(data=combined.dat,
                    gbm.x=c(9:27), 
                    gbm.y=8, # Recruitment Success
                    family="gaussian", tree.complexity=1, # (1 = no interactions)
                    learning.rate=0.01, bag.fraction=0.7)
(((BRT.Rs$self.statistics$mean.null)-(BRT.Rs$cv.statistics$deviance.mean))/(BRT.Rs$self.statistics$mean.null))*100 #61 % dev explained (63 without ppd/chla)
ri.Rs <- summary(BRT.Rs) # Had, jelly, zoo, ssb_1, mack_1 --> ppd and chla had zero relative abundance

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
