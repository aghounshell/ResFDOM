### Getting data together for 2020 WVWA Report back
### DOC, Abs, and FL data for DBP formation
### Focus on: long-term DOC data; SUVA and Humic C data for 2019

setwd("C:/Users/ahoun/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load long-term DOC data
fcrchem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DN_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  filter(Site == 50 | Site == 100 | Site == 200) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  #filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1|Depth_m == 9.0)

fcrchem <- fcrchem %>% select(-c("Reservoir")) %>% 
  rename(Date = DateTime, Depth = Depth_m, Station = Site)

fcrchem$Location <- paste(fcrchem$Station,fcrchem$Depth)

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

fcrchem <- completeFun(fcrchem,"DOC_mgL")

fcrchem_long <- fcrchem %>% filter(Station == 50 | Station == 100)

# Plot long-term DOC (add in oxygenation schedule?)
ggplot(fcrchem_long,mapping=aes(x=Date,y=DOC_mgL,color=as.factor(Location)))+
  annotate(geom="rect",xmin = as.POSIXct("2014-05-06"), xmax = as.POSIXct("2014-06-03"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2014-06-29"), xmax = as.POSIXct("2014-07-29"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2014-08-18"), xmax = as.POSIXct("2014-11-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-12-01"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  geom_line(size=1)+
  geom_point(size=3)+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  scale_color_manual(breaks=c('50 0.1','50 9','100 0.1'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Site 50 0.1m","Site 50 9m","Inflow"))+
  xlim(as.POSIXct("2014-01-01"),as.POSIXct("2019-12-31"))+
  #ylim(0,8)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggplot(fcrchem_long,mapping=aes(x=Date,y=DOC_mgL,color=as.factor(Location)))+
  geom_line(size=1)+
  geom_point(size=3)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  scale_color_manual(breaks=c('50 0.1','50 9','100 0.1'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Site 50 0.1m","Site 50 9m","Inflow"))+
  xlim(as.POSIXct("2019-01-01"),as.POSIXct("2019-12-31"))+
  ylim(0,7.5)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=7.5,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

### Load in fluorescence data (as fluorescent indicators, Peak C)
eems <- read_csv("./Data/20201103_ResultsFiles_ResEEMs2019.csv")

# Select FCR data for station 50, 100, 200
eems <- eems %>% filter(Reservoir == "FCR") %>% filter(Station == "50"|Station == "100"|Station == "200")
eems$Date <- as.POSIXct(strptime(eems$Date, "%m/%d/%Y", tz = "EST"))
m_eems <- eems %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_eems <- eems %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))
sd_eems_full <- eems %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

# Combine columns station and depth for plotting
# Combine columns (station, depth) for plotting
m_eems$Location <- paste(m_eems$Station,m_eems$Depth)
sd_eems$Location <- paste(sd_eems$Station,sd_eems$Depth)

# Plots
m_eems_50 <- m_eems %>% 
  filter(Station == "50" | Station == "100") %>% 
  filter(Depth == "0.1"|Depth == "9") %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21"))

sd_eems_50 <- sd_eems %>% 
  filter(Station == "50" | Station == "100") %>% 
  filter(Depth == "0.1"|Depth == "9") %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21"))

peakc <- ggplot(m_eems_50,mapping=aes(x=Date,y=C,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_eems_50,mapping=aes(ymin=m_eems_50$C-C,ymax=m_eems_50$C+C))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,0.36)+
  xlab("Date")+
  ylab("Peak C (RFU)")+
  scale_color_manual(breaks=c('50 0.1','50 9','100 0.1'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Site 50 0.1m","Site 50 9m","Inflow"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=0.35,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

### Calculate SUVA using abs values and DOC
# Load in abs data
abs <- read_csv("./Data/20200930_ResultsFiles_Abs2019.csv")

abs <- abs %>% filter(Reservoir == "FCR") %>% 
  filter(Station == "50"|Station == "100"|Station == "200") %>% 
  filter(Dilution == "1") %>% 
  filter(Depth == "0.1" | Depth == "9")
abs$Date <- as.POSIXct(strptime(abs$Date, "%m/%d/%Y", tz = "EST"))

m_abs <- abs %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_abs <- abs %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

m_abs$Location <- paste(m_abs$Station,m_abs$Depth)
sd_abs$Location <- paste(sd_abs$Station,sd_abs$Depth)

# Plot abs 254
### NOTE: NEED TO QA/QC THE DATA!
ggplot(m_abs,mapping=aes(x=Date,y=a254,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(sd_abs,mapping=aes(ymin=m_abs$a254-a254,ymax=m_abs$a254+a254))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("254 nm")+
  xlim(as.POSIXct("2019-04-01"),as.POSIXct("2019-11-30"))+
  scale_color_manual(breaks=c('50 0.1','50 9','100 0.1','200 0.1'),values=c("#7EBDC2","#393E41","#F0B670","#FE5F55"),labels=c("Site 50 0.1m","Site 50 9m","Inflow","Wetlands"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=70,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

#### Calculate SUVA
m_suva <- left_join(m_abs,fcrchem,by=c("Date","Station","Depth"))
m_suva <- m_suva %>% 
  select(Date,Station,Depth,a254,Location.x,DOC_mgL) %>% 
  rename(Location = Location.x)

m_suva <- m_suva %>% 
  mutate(suva = a254/2.303/DOC_mgL)

m_suva <- completeFun(m_suva,"suva")

m_suva_lim <- m_suva %>% filter(Station == "50"|Station == "100")

suva <- ggplot(m_suva_lim,mapping=aes(x=Date,y=suva,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("SUVA")+
  ylim(0,10)+
  xlim(as.POSIXct("2019-04-01"),as.POSIXct("2019-11-30"))+
  scale_color_manual(breaks=c('50 0.1','50 9','100 0.1'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Site 50 0.1m","Site 50 9m","Inflow"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=10,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(peakc,suva,ncol=1,nrow=2,common.legend = TRUE)
