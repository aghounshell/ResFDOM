### Getting data together for 2020 WVWA Report back
### DOC, Abs, and FL data for DBP formation
### Focus on: long-term DOC data; SUVA and Humic C data for 2019

setwd("C:/Users/ahoun/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

### Load in fluorescence data (as fluorescent indicators, Peak C, Peak T)
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
peakc <- ggplot(m_eems,mapping=aes(x=Date,y=C,color=(as.factor(Location))))+
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
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=0.35,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

peakt <- ggplot(m_eems,mapping=aes(x=Date,y=T,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_eems,mapping=aes(ymin=m_eems$T-T,ymax=m_eems$T+T))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,0.46)+
  xlab("Date")+
  ylab("Peak T (RFU)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=0.45,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(peakc,peakt,ncol=1,nrow=2)

hix <- ggplot(m_eems,mapping=aes(x=Date,y=HIX,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_eems,mapping=aes(ymin=m_eems$HIX-HIX,ymax=m_eems$HIX+HIX))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,10)+
  xlab("Date")+
  ylab("HIX")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=9.8,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bix <- ggplot(m_eems,mapping=aes(x=Date,y=BIX,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_eems,mapping=aes(ymin=m_eems$BIX-BIX,ymax=m_eems$BIX+BIX))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,0.9)+
  xlab("Date")+
  ylab("BIX")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=0.89,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(hix,bix,ncol=1,nrow=2)

### Load in PARAFAC model results
parafac <- read_csv("./Data/20201103_Site50_Mod4.csv")
parafac$Date <- as.POSIXct(strptime(parafac$Date, "%m/%d/%Y", tz = "EST"))
m_parafac <- parafac %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_parafac <- parafac %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

# Combine columns (station, depth) for plotting
m_parafac$Location <- paste(m_parafac$Station,m_parafac$Depth)
sd_parafac$Location <- paste(sd_parafac$Station,sd_parafac$Depth)

fmax1 <- ggplot(m_parafac,mapping=aes(x=Date,y=Fmax1,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_parafac,mapping=aes(ymin=m_parafac$Fmax1-Fmax1,ymax=m_parafac$Fmax1+Fmax1))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("Fmax1 (RFU)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1.5,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

fmax2 <- ggplot(m_parafac,mapping=aes(x=Date,y=Fmax2,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_parafac,mapping=aes(ymin=m_parafac$Fmax2-Fmax2,ymax=m_parafac$Fmax2+Fmax2))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("Fmax2 (RFU)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1.7,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(fmax1,fmax2,ncol=1,nrow=2)

fmax3 <- ggplot(m_parafac,mapping=aes(x=Date,y=Fmax3,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_parafac,mapping=aes(ymin=m_parafac$Fmax3-Fmax3,ymax=m_parafac$Fmax3+Fmax3))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("Fmax3 (RFU)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=3.2,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

fmax4 <- ggplot(m_parafac,mapping=aes(x=Date,y=Fmax4,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
  geom_errorbar(sd_parafac,mapping=aes(ymin=m_parafac$Fmax4-Fmax4,ymax=m_parafac$Fmax4+Fmax4))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab("Fmax4 (RFU)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(fmax3,fmax4,ncol=1,nrow=2)

### Load in DOC data
fcrchem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DN_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  filter(Site == 20| Site ==30| Site ==45 | Site == 50 | Site == 100 | Site ==200) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1|Depth_m == 5.0|Depth_m == 9.0)

fcrchem <- fcrchem %>% select(-c("Reservoir")) %>% 
  rename(Date = DateTime, Depth = Depth_m, Station = Site)

# Combine columns (station, depth) for plotting
fcrchem$Location <- paste(fcrchem$Site,fcrchem$Depth)

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

fcrchem <- completeFun(fcrchem,"DOC_mgL")

# Plot DOC data
ggplot(fcrchem,mapping=aes(x=Date,y=DOC_mgL,color=(as.factor(Location))))+
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
  ylab("DOC (mg/L)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=9,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

### Calculate SUVA using abs values and DO)C
# Load in abs data
abs <- read_csv("./Data/20200930_ResultsFiles_Abs2019.csv")

abs <- abs %>% filter(Reservoir == "FCR") %>% 
  filter(Station == "20"|Station =="30"|Station == "45"|Station == "50"|Station == "100"|Station == "200") %>% 
  filter(Dilution == "1") %>% 
  filter(Depth == "0.1" | Depth == "5" | Depth == "9")
abs$Date <- as.POSIXct(strptime(abs$Date, "%m/%d/%Y", tz = "EST"))

m_abs <- abs %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_abs <- abs %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

m_abs$Location <- paste(m_abs$Station,m_abs$Depth)
sd_abs$Location <- paste(sd_abs$Station,sd_abs$Depth)

# Plot abs 254
### NOTE: NEED TO QA/QC THE DATA!
a254 <- ggplot(m_abs,mapping=aes(x=Date,y=a254,color=(as.factor(Location))))+
  geom_line(size=1)+
  geom_point(size=4)+
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
  ylab("a254 (1/m)")+
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=65,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
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

suva <- ggplot(m_suva,mapping=aes(x=Date,y=suva,color=(as.factor(Location))))+
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
  scale_color_manual(breaks=c('100 0.1','200 0.1','20 0.1','30 0.1','45 0.1','50 0.1','50 5','50 9'),values=c("#DA2C38","#E7804B","#EEAA55","#F4D35E","#7FC6A4","#7EBDC2","#5C7E82","#393E41"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=28,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(a254,suva,ncol=1,nrow=2)
