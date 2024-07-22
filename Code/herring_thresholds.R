dev.exp <- function(model){
  (((model$self.statistics$mean.null)-(model$cv.statistics$deviance.mean))/(model$self.statistics$mean.null))*100 
}


# make a fx for this



# Explore seasonal oisst using the data i already have
OISST%>%
  dplyr::filter(EPU==c("GOM"))%>%
  ggplot(aes(x=day, y = Value, col = as.factor(year)))+
  geom_line(na.rm = T) +
  #facet_wrap(~EPU)+
  scale_x_continuous(breaks = seq(0,365,30),minor_breaks = seq(0,365,5)) +
  ylab("OISST")+
  theme_bw()+
  geom_hline(yintercept=c(16,20))+
  geom_vline(xintercept=c(90, 181, 273))

# will need to maybe redo this with spatial stuff but for now real quick calculate number of days above 16 for each year
OISST%>%
  dplyr::filter(EPU==c("GOM"))%>%
  group_by(year)%>%
  summarize(number = sum(Value>=16))

testing <- OISST%>%
  dplyr::filter(EPU==c("GOM"))%>%
  group_by(year)%>%
  summarize(number.opt = sum(dayValue>=16),
            number.limit = sum(Value>=20))

ggplot(testing, aes(x=year, y = number.opt))+
  geom_point()+
  geom_line()
  
# having projection issues using joes method
# joe and abby have a similar but slightly different way of doing a similar thing
# I'm getting stuck on lack of overlap in the data for some reason
###Abby stuff for number of days function
# need a GoM strata geometry thing

GOM_strat <- ecodata::epu_sf%>%
  dplyr::select(EPU, geometry) %>%
  filter(EPU == "GOM")



