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