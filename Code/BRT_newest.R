##      Title:            Herring Recruitment Boosted Regression Trees
##      Author:           Adelle Molina
##      Created:          4/28/23
##      Updated:          8/19/24
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
#remotes::install_github("gbm-developers/gbm3", build_vignettes = TRUE, force = TRUE)
library(gbm3)

# Load data ---------------------------------------------------------------
names(newdat)

# Try with other gbm functions (instead of step) ---------------------------

# here using function gbmt
train_params <- training_params(num_trees = 1000, interaction_depth = 1, min_num_obs_in_node = 5,num_features = 3)
gauss_fit <- gbmt(recrt~BT.GB_1+BT.GOM_1+c5fall.gom_1+fall.sm_1+spring.lg+HW.GOM_1,
                  data=newdat, cv_folds =5, keep_gbm_data = TRUE, train_params = train_params, distribution = "gaussian")
#summary(gauss_fit)

# here using function gbm.fit
gbm2 <- gbm(logRdev~BT.GB_1+BT.GOM_1+
              OISST.GB + OISST.GOM+
              spring.lg + fall.sm_1 + HadSSB_GB+ HadSSB_GoM + HadSSB_GB_1 + HadSSB_GoM_1 + c5fall.gom_1+
              HW.GB + HW.GOM,         # formula (all columns with . otherwise type out formula)
            data=newdat,                   # dataset
            distribution="gaussian",     
            n.trees=50,                # number of trees
            shrinkage=0.001,             # shrinkage/learning rate, 0.001 to 0.1 usually work
            interaction.depth=1,         # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
            #train.fraction = 0,        # fraction of data for training, 
            #mFeatures = 3,               # Number of features to consider at each node.
            n.minobsinnode = 5,         # minimum number of obs needed in each node (10 is too high for dataset size)
            keep.data=TRUE,
            cv.folds=10,                 # do 10-fold cross-validation
            verbose = F) 

gbm2 # Evaluate model fit
gbm.perf(gbm2, method = "cv") # this shows loss function
summary(gbm2)
# "better" models when I use logRdev as opposed to raw recruitment values

# how do I compute deviance explained? I think one of the two methods below works (one uses all trees the other uses just one)
abs((last(gbm2$train.error)-last(gbm2$cv_error))/last(gbm2$train.error))*100
(sum(gbm1$train.error)-sum(gbm1$cv_error))/sum(gbm1$train.error)*100

# pretty plot of relative influence
effects <- tibble::as_tibble(summary(gbm3))
effects %>% 
  # arrange descending to get the top influencers
  dplyr::arrange(desc(rel_inf)) %>%
  # plot these data using columns
  ggplot(aes(x = forcats::fct_reorder(.f = var, 
                                      .x = rel_inf), 
             y = rel_inf, 
             fill = rel_inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title = element_text()) + 
  xlab('Features') +
  theme_bw()+
  ylab('Relative Influence') +
  ggtitle("Top Drivers of Recruitment")

# Play with gbm3 functions to eval model ----------------------------------

# explore the model and other functions in the pacakge
plot(gbm2, var_index = 5)
#vip::vip(gbm2) 
str(gbm2)
#permutation_relative_influence(gbm3, num_trees = 50) # i think this is to add trees --> do if you need more iteration

# predict to the data, 
gbm_pred <- predict(object = gbm2,  newdata = newdat,n.trees = 50, type="link")
# trying to plot the fit relative to the actual data but the scales are wrong (multiply by 100 got us there even though that prob isn't right)
ggplot(newdat, aes(x=year))+
  geom_line(aes(y=100*gbm2$fit))+
  geom_line(aes(y=logRdev), col="red")


best.iter <- gbm.perf(gbm2, method = "cv")

Yhat <- predict(gbm2, newdata = newdat, n.trees = best.iter, type = "link")
print(sum((newdat$logRdev - Yhat)^2))# least squares error

# plot marginal effectfor one variable at a time
plot(gbm3)



# GBM3 without correlated vars --------------------------------------------
# added gbm4 with interaction depth set to 2....

# Try with correlated variables removed
gbm4 <- gbm(logRdev~BT.GB_1+BT.GOM_1+ OISST.GB + OISST.GOM+
              spring.lg + fall.sm_1 +  HadSSB_GB_1 + HadSSB_GoM_1 +
              HW.GB + HW.GOM,        
            data=newdat,                   # dataset
            distribution="gaussian",     
            n.trees=5000,                
            shrinkage=0.0005,             
            interaction.depth=1,         
            bag.fraction = 0.5,                
            n.minobsinnode = 5,         
            keep.data=TRUE,
            cv.folds=10,                 
            verbose = F) 
gbm.perf(gbm4, method = "cv") # this shows loss function
sqrt((gbm4$train.error[1]-gbm4$cv_error[1])/gbm4$train.error[1])
hist(newdat$logRdev)
sqrt(gbm4$cv_error[1])
sqrt(gbm4$train.error[1])
gbm4 # Evaluate model fit (lower pseudo R than gbm2)
gbm.perf(gbm4, method = "cv") # this shows loss function
summary(gbm4)
abs((last(gbm3$train.error)-last(gbm3$cv_error))/last(gbm3$train.error))*100 # very low fake dev explained

# Plot that
effects_gbm4 <- tibble::as_tibble(summary(gbm4))
effects_plot <- effects_gbm4 %>% 
  # arrange descending to get the top influencers
  dplyr::arrange(desc(rel_inf)) %>%
  # plot these data using columns
  ggplot(aes(x = forcats::fct_reorder(.f = var, 
                                      .x = rel_inf), 
             y = rel_inf, 
             fill = rel_inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title = element_text()) + 
  xlab('Features') +
  theme_bw()+
  ylab('Relative Influence') +
  ggtitle("Top Drivers of Recruitment")

effects_plot

ggsave("newBRT.png", width = 8, height = 6)

# using the parameters from gbm.fit, rerun using gbm.step with those values
BRT.mm192_opt <-gbm.step(data=newdat,
                         gbm.x=c(5:16), minobsinnode=5,
                         gbm.y=4, # mm192 log recruitment deviations
                         family="gaussian", tree.complexity=1, # (1 = no interactions)
                         learning.rate=0.0001, bag.fraction=0.5, n.trees = 1000, n.folds = 10) # doesn't run? 

# BRT with WHAM output and Sarah zoo indices ------------------------------
BRT.mm192_raw <-gbm.step(data=newdat,
                         gbm.x=c(5:ncol(newdat)), 
                         gbm.y=2, # mm192 recruits
                         family="gaussian", tree.complexity=1, # (1 = no interactions)
                         learning.rate=0.05, bag.fraction=0.7)
dev.exp(BRT.mm192_raw)

BRT.mm192 <-gbm.step(data=newdat,
                          gbm.x=c(5:ncol(newdat)), 
                          gbm.y=4, # mm192 log recruits
                          family="gaussian", tree.complexity=1, # (1 = no interactions)
                          learning.rate=0.00001, bag.fraction=0.7) # why is LR so diff from raw model
dev.exp(BRT.mm192)
#WHAM$logR <- as.integer(WHAM$logR)

BRT.mm192_bio <-gbm.step(data=newdat,
                     gbm.x=c(5:ncol(newdat)), 
                     gbm.y=4, # mm192 log deviations
                     family="gaussian", tree.complexity=1, # (1 = no interactions)
                     learning.rate=0.000005, bag.fraction=0.8)


#(((BRT.mm192_raw$self.statistics$mean.null)-(BRT.mm192_raw$cv.statistics$deviance.mean))/(BRT.mm192_bio$self.statistics$mean.null))*100 
(1.114876e+13 -1.08101e+13)/1.114876e+13 # very low calculated deviance explained for this model

# BRT for log recruitment deviations ---------
names(combined.dat)
ncol(combined.dat)
summary(BRT.mm192_raw)
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
