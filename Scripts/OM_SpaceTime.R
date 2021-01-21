### Getting data together for 2020 WVWA Report back
### DOC, Abs, and FL data for DBP formation
### Focus on: long-term DOC data; SUVA and Humic C data for 2019

setwd("C:/Users/ahoun/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

### Load in fluorescence data (as fluorescent indicators, Peak C)
eems <- read_csv("./Data/20201103_ResultsFiles_ResEEMs2019.csv")

# Select FCR data for station 50, 100, 200
eems <- eems %>% filter(Reservoir == "FCR") %>% filter(Station == "20"|Station == "30"|Station == "45"|Station == "50"|Station == "100"|Station == "200")
eems$Date <- as.POSIXct(strptime(eems$Date, "%m/%d/%Y", tz = "EST"))
m_eems <- eems %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_eems <- eems %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

# Combine columns station and depth for plotting
# Combine columns (station, depth) for plotting
m_eems$Location <- paste(m_eems$Station,m_eems$Depth)
sd_eems$Location <- paste(sd_eems$Station,sd_eems$Depth)

# Plots
ggplot(m_eems,mapping=aes(x=Date,y=C,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_eems,mapping=aes(ymin=m_eems$C-C,ymax=m_eems$C+C))+
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
  scale_color_manual(breaks=c('20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9','100 0.1','200 0.1'),values=c("#F4D35E","#7FC6A4","#BFACC8","#7EBDC2","#5C7E82","#393E41","#EEAA55","#E7804B"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=0.35,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())


values=c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#E7804B","#DA2C38","#7DAF4B")

### Calculate SUVA using abs values and DO)C
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
