##      Title:            Herring Recruitment Boosted Regression Trees
##      Author:           Adelle Molina
##      Created:          4/28/23
##      Updated:          12/1/24
##      Notes:            This script has most recent model runs
##      To Do:            Clean this up, a tiny bit more tweaking/tuning, not sure where I had left off
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
# use data called newdat from data preparation script
names(newdat)

# Use gbmt function (from gbm3, instead of gbm.step) ---------------------------
names(newdat)
train_params <- training_params(num_trees = 50, interaction_depth = 2, min_num_obs_in_node = 5, shrinkage = 0.001, bag_fraction = 0.5)
gauss_fit <- gbmt(logRdev ~ BT.GB_1+BT.GOM_1+ NAOIndex_3+OISST.GB + OISST.GOM+
                  spring.lg + fall.sm_1 +c5fall.gom_1+ Had_1 + jelly_GB_1 + jelly_GOM_1,   
                  data=newdat, cv_folds = 10, keep_gbm_data = TRUE, train_params = train_params, distribution = gbm_dist("Laplace"))
summary(gauss_fit)
best_iter_cv <- gbmt_performance(gauss_fit, method='cv')
plot(best_iter_cv)

# using gbm A (from gbm3) --------------------------------------
# several runs with WHAM and VAST and a range of other vars
gbm2 <- gbm3::gbm(logRdev~BT.GB_1+BT.GOM_1+
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

# how do I compute deviance explained? 
#I think one of the methods below works (one uses all trees the other uses just the last tree)
#abs((last(gbm2$train.error)-last(gbm2$cv_error))/last(gbm2$train.error))*100
#(sum(gbm2$train.error)-sum(gbm2$cv_error))/sum(gbm2$train.error)*100
#(mean(gbm2$train.error)-mean(gbm2$cv_error))/mean(gbm2$train.error)*100
abs(mean((gbm2$train.error-gbm2$cv_error)/gbm2$train.error))*100 # --> turned this one into a function

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

# Using gbm B (from gbm3) with WHAM, VAST & other "smart" lags ---------------------------------------------------------------------

# Using new (as of 9/26) "smart" lags
brt_smart1 <- gbm3::gbm(logR~BT.GB_1+BT.GOM_1+ NAOIndex_5+
                    OISST.GB + OISST.GOM+
                    spring.lg + fall.sm_1 + HadSSB_GB_1 + HadSSB_GoM_1 + c5fall.gom_1+
                    jelly_GB_1 + jelly_GOM_1,         # formula (all columns with . otherwise type out formula)
                  data=newdat, distribution="gaussian",     
                  n.trees=50,                # number of trees
                  shrinkage=0.001,             # shrinkage/learning rate, 0.001 to 0.1 usually work
                  interaction.depth=1,         # 1: additive model, 2: two-way interactions, etc
                  bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.8,        # fraction of data for training, 
                  #mFeatures = 3,               # Number of features to consider at each node.
                  n.minobsinnode = 5,         # minimum number of obs needed in each node (10 is too high for dataset size)
                  keep.data=TRUE,
                  cv.folds=5,                 # do 10-fold cross-validation
                  verbose = F) 

brt_smart1 # Evaluate model fit
gbm.perf(brt_smart1, method = "cv") # this shows loss function
summary(brt_smart1)
dev.exp2(brt_smart1)
# calculate deviance explained
(sum(brt_smart1$train.error)-sum(brt_smart1$cv_error))/sum(brt_smart1$train.error)*100

# Using newest (as of 10/31) "smart" lags and Micah haddock index instead of SSB
brt_smart2 <- gbm3::gbm(logRdev~BT.GB_1+BT.GOM_1+ NAOIndex_5+OISST.GB + OISST.GOM+
                        spring.lg + fall.sm_1 + Had_1 +jelly_GB_1 + jelly_GOM_1, 
                        data=newdat, distribution="gaussian",     
                        n.trees=300,  shrinkage=0.005,             
                        interaction.depth=2,  n.minobsinnode = 5,      
                        bag.fraction = 0.7,   cv.folds = 2, train.fraction = 0.8,      
                        keep.data=TRUE, verbose = F) 
gbm.perf(brt_smart2, method = "cv") 
gbmt_performance(brt_smart2, method='test')
summary(brt_smart2)
brt_smart2 
dev.exp2(brt_smart2) # this doesn't seem right it was too low now its too high....

# Using gbm C (from gbm3) with WHAM, VAST, Micah, thermal indicators, strata means (DEVEL)---------------------------------------------------------------------
# rerun with wint/fall index, one of my duration jawns and/or mean sst and bt
names(newdat)
# What is the best shrinkage rate? I have had it as low as 0.001, had determined previously to use 0.05, but am now using 0.01

# Using newest (as of 11/18) indicator list...also switching to using logR since I don't think I did devs right
# 11/19 added back in no lag SST
# Run this one last time and call it a damn day ok
herr_logR_gbm <- gbm3::gbm(logR ~  avg.fall.BT_1  + avg.fall.SST + days.optimal.SST_1 + 
                        fall.wint.sm_1  + #NAOIndex_3 +fall.sm_1 +avg.fall.SST_1 +
                        spring.lg +  Had_1 +jelly_GB_1 + jelly_GOM_1 +
                        WSW + NAOIndex_5 + GSI , 
                        data=newdat, distribution="gaussian",     
                        n.trees=250,  shrinkage=0.01,             
                        interaction.depth=2,  n.minobsinnode = 5,      
                        bag.fraction = .8,   cv.folds = 10, train.fraction = 0.8,      
                        keep.data=TRUE, verbose = F) 

gbm3::gbm.perf(herr_logR_gbm, method = "cv") # plot error functions
gbm3::gbm.perf(herr_logR_gbm, method='test')
#gbmt_performance(herr_logR_gbm, method='OOB')

# Save optimal number of trees for fitting
opt.trees <- gbm.perf(herr_logR_gbm, method = "cv")
opt.trees.test <- gbm.perf(herr_logR_gbm, method = "test")

# Calculate mean errors and look at model
mean(herr_logR_gbm$cv_error)
dev.exp2(herr_logR_gbm)
dev.exp3(herr_logR_gbm)
print(herr_logR_gbm)
summary(herr_logR_gbm)

#sqrt(mean(herr_logR_gbm$cv_error))
#mean(herr_logR_gbm$train.error)
#mean(herr_logR_gbm$valid.error)
#(abs(herr_logR_gbm$train.error-herr_logR_gbm$cv_error)/herr_logR_gbm$train.error)*100
# Determine optimal number of trees (loop) --------------------------------
# Write a loop to summarize error loss as trees are added once shrinkage rate has been set

cv.results.n10 <- c()
for(i in seq(100, 4000, 100)){
  mod <- gbm3::gbm(logR ~  avg.fall.BT_1  + avg.fall.SST_1 + spring.lg + fall.sm_1 + Had_1 +jelly_GB_1 + jelly_GOM_1 +
                  WSW + NAOIndex_3 + NAOIndex_5 + GSI , data=newdat, distribution="gaussian",     
                  n.trees=i,  shrinkage=0.01,  interaction.depth=2,  n.minobsinnode = 5,      
                  bag.fraction = 0.7,   cv.folds = 10, train.fraction = 0.8, keep.data=TRUE, verbose = F) 
  fit <- data.frame(year=newdat$year, fit=mod$fit,  y=newdat$logR)
  results <- fit%>%
    dplyr::summarize(SSE = sum((fit - y)^2),
                     SSE.SD=sd((fit - y)^2),
                     RMSE = sqrt(mean(mod$cv_error)))%>% 
    mutate(N.trees=i)  
  cv.results.n10 <-  rbind(cv.results.n10, results)
}
# results for cv.folds = 10
goal.error <- min(cv.results.n10$SSE)+mean(cv.results.n10$SSE.SD)
ggplot(cv.results.n10, aes(x = n.trees.fix, y = SSE))+
  geom_point()+
  theme_bw()+ 
  geom_errorbar(aes(ymin=SSE-SSE.SD, ymax=SSE+SSE.SD), width=.2, 
                position=position_dodge(0.05))+
  geom_hline(yintercept = goal.error)

# results for when cv.folds = 2
cv.results.n10$n.trees.fix <- seq(100, 4000, 100)
ggplot(cv.results, aes(x = n.trees.fix, y = SSE))+
  geom_point()+
  theme_bw()+ 
  geom_errorbar(aes(ymin=SSE-SSE.SD, ymax=SSE+SSE.SD), width=.2, 
                position=position_dodge(0.05))+
  geom_hline(yintercept = goal.error)
goal.error <- min(cv.results$SSE)+mean(cv.results$SSE.SD)
# this exercise (based in part on leathwick complexity appendix section) suggests that 800-1000 trees is best
 #using more folds got similar result but with slightly less than 1000
# diff ways of estimating optimal number of trees
# but the method above suggests 120-200 is the optimal number of trees

# Plot and Model Diagnostics --------------------------------------------------------------

# predict back onto the data using the optimum number of trees to visualize the fit 
Yhat <- predict(herr_logR_gbm, newdata = newdat, n.trees = opt.trees)
Yhat.test <- predict(herr_logR_gbm, newdata = newdat, n.trees = opt.trees.test)
fitted  <- data.frame(year=newdat$year, brt.pred=Yhat, expect = newdat$logR)
fitted.test  <- data.frame(year=newdat$year, brt.pred=Yhat.test, expect = newdat$logR)
#fitted  <- fitted %>%mutate(devs = abs(expect-brt.pred))

print(sum((newdat$logR - Yhat)^2))# least squares error
#gbm_roc_area(newdat$logR, Yhat) # hmmm this is always .5 even when I use diff data...does tha tmake sense its diff nrows

# Plot obs vs fitted
ggplot(fitted, aes(x=year, y=expect))+ 
  theme_bw()+  
  geom_line(na.rm=T, col="blue")+ #observed data
  #geom_point(data=fitted, mapping=aes(x=year, y=brt.pred)) +
  geom_line(aes(y=brt.pred), lwd=1) + # fitted data (all years)
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +  ylab("log Recruitment")
# add a legend

# Model fit to model data
# Extract fit --> this is same same but diff b/c it's using final number of trees I think
mod.fit <- data.frame(year=newdat$year, fit=herr_logR_gbm$fit,  y=newdat$logR)
#mod.fit <- data.frame(year=newdat$year, fit=herr_logR_gbm$fit,  y=newdat$logR)
print(sum((mod.fit$fit - mod.fit$y)^2))# least squares error of the last model
print(sd((mod.fit$fit - mod.fit$y)^2))# least squares error 

# now add in fits from the model output in red
ggplot(fitted, aes(x=year, y=expect))+ 
  theme_bw()+  
  geom_line(na.rm=T, col="blue")+ #observed data
  geom_line(aes(y=brt.pred), lwd=1) + # fit (optimal trees)
  #geom_point(data=mod.fit, mapping=aes(x=year, y=fit), col="red") + # fit (model)
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  xlab("Year") +  ylab("log Recruitment")
# how do I add error barss here or even a cloud...I'd have to repeat or look at residuals

# Compare with short data -------------------------------------------------
# Repeat as above but use a shortened dataset --> ope add back in the non lagged sst
newdat.short <- newdat %>% filter(between(year,1987,2020))
logR_gbm_short <- gbm3::gbm(logR ~  avg.fall.BT_1  + avg.fall.SST_1 +days.optimal.SST_1 + avg.fall.SST+ 
                             spring.lg + fall.sm_1 + Had_1 +jelly_GB_1 + jelly_GOM_1 +
                             WSW + NAOIndex_5 + GSI , 
                           data=newdat.short, distribution="gaussian",     
                           n.trees=250,  shrinkage=0.01,             
                           interaction.depth=2,  n.minobsinnode = 5,      
                           bag.fraction = .8,   cv.folds = 10, train.fraction = 0.9,      
                           keep.data=TRUE, verbose = F) 

gbm.perf(logR_gbm_short, method = "cv") 
opt.trees.short <- gbm.perf(logR_gbm_short, method = "cv")

# predict onto the full data using the optimum number of trees  
Yhat.short <- predict(logR_gbm_short, newdata = newdat, n.trees = opt.trees.short)
fitted.short  <- data.frame(year=newdat$year, brt.pred=Yhat.short, expect = newdat$logR)

# Plot the fit to the short data and make predictions for the full data

ggplot(fitted.short, aes(x=year, y=expect))+ # plot obs vs pred
  theme_bw()+  
  geom_line(na.rm=T, aes(color=""),lwd=1)+ #observed data
  geom_line(aes(color=" ", y=brt.pred), lwd=1, show.legend=T) + # fitted data (short)
  #geom_point(data=mod.fit, mapping=aes(x=year, y=fit), col="red") + # fitted data (short)
  scale_x_continuous(breaks = seq(0,2022,5),minor_breaks = seq(0,2022,1))+
  guides(color=guide_legend(position='inside'))+
  scale_color_manual(labels=c("Observed", "Predicted"),
                    values = c( "blue", "black"))+
  theme(legend.position.inside = c(.2,.2),
        legend.background = element_blank(),
        legend.title = element_blank())+
  xlab("Year") +  ylab("log Recruitment")
ggsave("Fitted.11.18.png", width = 8, height = 6)

# Other Figures & Diagnostics -----------------------------------------------------------

permutation_relative_influence(herr_logR_gbm, rescale = T,  sort_it = T, num_trees = opt.trees)

# pretty plot of relative influence
ri <- summary(herr_logR_gbm)
effects <- ri %>% 
  dplyr::arrange(desc(rel_inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var, 
                                      .x = rel_inf), 
             y = rel_inf, 
             fill = rel_inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.title = element_text()) + 
  xlab('Variables') +
  theme_bw()+
  ylab('Relative Influence') +
  ggtitle("Top Predictors of log Recruitment")
ggsave("PartialDep_11.18.png", width = 8, height = 6)

# calibration plot...whatever that is
calibrate_plot(newdat$logR, Yhat, "Gaussian")

# Explore interactions
interact(herr_logR_gbm, newdat, var_indices=c(1,5), opt.trees)
# is there a way to only plot particular interactions....like how do we choose the ones that had high var importance....oh reorder
# can I write a loop to find any sig interactions
for(i in 1:nrow(ri)){
  pdp<-plot(herr_logR_gbm,c(i,),return_grid = T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(ri$rel_inf[which(row.names(ri)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","y","var","ri")
  x<-rbind(x,pdp)
}

# Plot Partial Dependence for chosen model --------------------------------
# Figure out syntax
plot(herr_logR_gbm, var_index=c(1)) # can do one at a time...but they're in order of the model input
plot(herr_logR_gbm, var_index=c(4,6)) # can also do two at a time, which shows an interesting grid type figure...
#plot(herr_logR_gbm, var_index=c(1,2,3))#if you do three at a time you get a grid of grids

# Loop to plot only the top 6 variables, maybe make it more like 8?
ri <- summary(herr_logR_gbm)
x<-c()
for(i in 1:nrow(ri)){
  pdp<-plot(herr_logR_gbm,i,return_grid = T)
  pdp$var<-rep(colnames(pdp)[1],nrow(pdp))
  pdp$ri<-rep(ri$rel_inf[which(row.names(ri)==colnames(pdp)[1])],nrow(pdp))
  colnames(pdp)<-c("val","y","var","ri")
  x<-rbind(x,pdp)
}

# order by relative influence
x$var<-factor(x$var,levels=row.names(ri))

# top 6
top6<-as.character(ri$var[1:6])
relinf<-round(ri$rel_inf[1:6],1)

x%>%
  filter(var==top6[1])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  theme_bw()+
  labs(x=paste(top6[1],"(",relinf[1],"%)",sep=""),y="Marginal effect")+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p1

x%>%
  filter(var==top6[2])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[2],"(",relinf[2],"%)",sep=""),y="")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p2

x%>%
  filter(var==top6[3])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[3],"(",relinf[3],"%)",sep=""),y="")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p3
x%>%
  filter(var==top6[4])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[4],"(",relinf[4],"%)",sep=""),y="Marginal effect ")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p4

x%>%
  filter(var==top6[5])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[5],"(",relinf[5],"%)",sep=""),y="")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p5

x%>%
  filter(var==top6[6])%>%
  ggplot(aes(x=val,y=y))+
  geom_smooth(method="loess",se=T,color="darkgray")+
  labs(x=paste(top6[6],"(",relinf[6],"%)",sep=""),y="")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=12))->p6
pdpz <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
ggsave("PartialDep_11.18.png", width = 8, height = 6)


tiff("Figures/ts.Competitors.tiff", width = 8, height = 5, units = 'in', res = 300)
ggarrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()

# Using gbm A2(from gbm3) with reduced dataset --------------------------------------------
# added gbm4 with interaction depth set to 2....
# this is where I had previously left off and was playing around and tuning 

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
gbm.perf(gbm4, method = "cv") # this shows loss function, but clearly haven't gotten to that point of max error loss
sqrt((gbm4$train.error[1]-gbm4$cv_error[1])/gbm4$train.error[1])
hist(newdat$logRdev)
sqrt(gbm4$cv_error[1])
sqrt(gbm4$train.error[1])
gbm4 # Evaluate model fit (lower pseudo R than gbm2)
summary(gbm4)
abs((last(gbm4$train.error)-last(gbm4$cv_error))/last(gbm4$train.error))*100 # very low fake dev explained, very high for model with interaction depth 2

# Plot relative influence
effects_gbm4 <- tibble::as_tibble(summary(gbm4))
effects_plot <- effects_gbm4 %>% 
  dplyr::arrange(desc(rel_inf)) %>%
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

# Play with gbm3 functions to eval model ----------------------------------

# explore the model and other functions in the pacakge
plot(gbm2, var_index = 6)

# predict to the data, 
gbm_pred <- predict(object = gbm2,  newdata = newdat,n.trees = 50, type="link")
# trying to plot the fit relative to the actual data but the scales are wrong
#(multiply by 100 got us there even though that prob isn't right)
ggplot(newdat, aes(x=year))+
  geom_line(aes(y=100*gbm2$fit))+
  geom_line(aes(y=logRdev), col="red")

best.iter <- gbm.perf(gbm2, method = "cv")

Yhat <- predict(gbm2, newdata = newdat, n.trees = best.iter, type = "link")
print(sum((newdat$logRdev - Yhat)^2))# least squares error

# plot marginal effectfor one variable at a time
plot(gbm3)

# BRT with WHAM and VAST (gbm.step) ------------------------------
BRT.mm192_step <-gbm.step(data=newdat,gbm.x=c(5:ncol(newdat)),  gbm.y=3, # mm192 recruits (log)
                           family="gaussian", tree.complexity=1,
                          learning.rate=0.001, bag.fraction=0.8, train.fraction = 0.8, n.minobsinnode=5, step.size = 1) 

dev.exp(BRT.mm192_smart) #
# error that dataset is too small or subsampling rate is too large
# other error to restart model with a smaller learning rate or smaller step size...


# ok so these dismo gbm.step models just are not working at all...
BRT.mm192_smart <-gbm.step(data=newdat,
                         gbm.x=c(5:ncol(newdat)), 
                         gbm.y=3, # mm192 recruits (log)
                         family="gaussian", tree.complexity=2, # (1 = no interactions)
                         learning.rate=0.00001, bag.fraction=0.5) 

dev.exp(BRT.mm192_smart) #
names(newdat)
effects <- tibble::as_tibble(summary(BRT.mm192_smart))
effects %>% 
  dplyr::arrange(desc(rel.inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var, 
                                      .x = rel.inf), 
             y = rel.inf, 
             fill = rel.inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  xlab('Features') +
  theme_bw()+
  ylab('Relative Influence') 


# "older"
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

# BRT for log recruitment deviations ASAP (gbm.step) ---------
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
