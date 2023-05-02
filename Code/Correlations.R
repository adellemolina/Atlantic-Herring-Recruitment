##      Title:            Correlations
##      Script Purpose:   Explore correlations in multivariate time series (auto, multicollinearity)
##      Author:           Adelle Molina
##      Created:          5/2/22
##      Updated:          5/2/23
##      To Do:            This is a mess of code from the end of script



# Correlation ---------------------------------------------------------
# Environmental Data
temp.vars <- basin.vars[, c(1:3, 5,7:10)]
temp.vars2 <- reg.vars%>% 
  dplyr::select(year,OISST.GB,OISST.GOM,OISST.MAB,OISST.SS,
                BT.GB, BT.GOM, BT.MAB, BT.SS,
                SST.GB, SST.GOM, SST.MAB, SST.SS)
temp.vars <- temp.vars%>% full_join(temp.vars2, by="year")
tempvar.melt <- reshape2::melt(temp.vars, id.vars="year", variable.name="series")


tempcor <- ggstatsplot::ggcorrmat(
  data = temp.vars,
  type = "nonparametric", 
  colors = c("darkred", "white", "steelblue"))

# old autocorrelation exploration (not updated)
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

# Combined correlation matrix (not updated) ---------------------------------------------

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

# Grouped correlation matrices (not updated) --------------------------------------------

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


