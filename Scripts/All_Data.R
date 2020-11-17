## Script to look at data and collate for PCA/Time-series analysis
## A Hounshell, 03 Nov 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,PerformanceAnalytics)

# Load PARAFAC data
parafac <- read_csv("C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20201103_Site50_Mod4.csv")
parafac$Date <- as.POSIXct(strptime(parafac$Date, "%m/%d/%Y", tz = "EST"))
m_parafac <- parafac %>% group_by(Date,Station,Depth) %>% summarise_all(funs(mean))
sd_parafac <- parafac %>% group_by(Date,Station,Depth) %>% summarise_all(funs(sd))

# Plot SD for each sample?
plot(sd_parafac$Date,sd_parafac$Fmax1)

# Combine columns (station, depth) for plotting
m_parafac$Location <- paste(m_parafac$Station,m_parafac$Depth)
sd_parafac$Location <- paste(sd_parafac$Station,sd_parafac$Depth)

# Plot by location
ggplot(m_parafac,mapping=aes(x=Date,y=Fmax1,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_parafac,mapping=aes(x=Date,y=Fmax1,ymin=Fmax1-sd_parafac$Fmax1,ymax=Fmax1+sd_parafac$Fmax1,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

ggplot(m_parafac,mapping=aes(x=Date,y=Fmax2,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_parafac,mapping=aes(x=Date,y=Fmax2,ymin=Fmax2-sd_parafac$Fmax2,ymax=Fmax2+sd_parafac$Fmax2,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

ggplot(m_parafac,mapping=aes(x=Date,y=Fmax3,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_parafac,mapping=aes(x=Date,y=Fmax3,ymin=Fmax3-sd_parafac$Fmax3,ymax=Fmax3+sd_parafac$Fmax3,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

ggplot(m_parafac,mapping=aes(x=Date,y=Fmax4,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_parafac,mapping=aes(x=Date,y=Fmax4,ymin=Fmax4-sd_parafac$Fmax4,ymax=Fmax4+sd_parafac$Fmax4,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

##########################
# Get data together to merge with others
m_parafac_50 <- m_parafac %>% select(Date,Station,Depth,Rep,Fmax1,Fmax2,Fmax3,Fmax4,Location) %>% 
  filter(Station == 50)
m_parafac_100 <- m_parafac %>% select(Date,Station,Depth,Rep,Fmax1,Fmax2,Fmax3,Fmax4,Location) %>% 
  filter(Station == 100)
m_parafac_200 <- m_parafac %>% select(Date,Station,Depth,Rep,Fmax1,Fmax2,Fmax3,Fmax4,Location) %>% 
  filter(Station == 200)
############################

# Load EEMs data (HIX, BIX, etc.)
eems <- read_csv("C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20201103_ResultsFiles_ResEEMs2019.csv")

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

# Plot by location
ggplot(m_eems,mapping=aes(x=Date,y=HIX,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_eems,mapping=aes(x=Date,y=HIX,ymin=HIX-sd_eems$HIX,ymax=HIX+sd_eems$HIX,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

ggplot(m_eems,mapping=aes(x=Date,y=BIX,color=Location))+
  geom_line()+
  geom_point()+
  geom_errorbar(m_eems,mapping=aes(x=Date,y=BIX,ymin=BIX-sd_eems$BIX,ymax=BIX+sd_eems$BIX,color=Location))+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-21"))+
  theme_classic(base_size=15)

# Plots for AGU
m_eems_50 <- m_eems %>% 
  filter(Station == "50") %>% 
  filter(Depth == "0.1"|Depth == "9") %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21"))

sd_eems_50 <- sd_eems %>% 
  filter(Station == "50") %>% 
  filter(Depth == "0.1"|Depth == "9") %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21"))

bix <- ggplot(m_eems_50,mapping=aes(x=Date,y=BIX,color=(as.factor(Depth))))+
  geom_hline(yintercept = 0.7, color = "grey")+
  geom_hline(yintercept = 0.8, color = "grey")+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(sd_eems_50,mapping=aes(ymin=m_eems_50$BIX-BIX,ymax=m_eems_50$BIX+BIX))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,1)+
  xlab("Date")+
  scale_color_manual(breaks=c('0.1','9'),values=c("#7EBDC2","#393E41"),labels=c("Epi","Hypo"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

#################### Combine PARAFAC and EEMs data by location ##############################
fl_data <- full_join(m_parafac,m_eems,by=c("Date","Location","Station","Depth"))
fl_data <- fl_data %>% select(-c("SampleName","Reservoir.x","i","Rep.x","Rep.y","Sample Name","Reservoir.y","Date analyzed","Max Fl Ex","Max Fl Em","Max Fl"))
fl_data <- fl_data %>% relocate("Location",.after = "Depth")
fl_data <- fl_data %>% relocate("Dilution",.after= "Location")
fl_data <- fl_data %>% mutate(Date = as.Date(Date,"%Y-%m-%d")) %>% arrange(Date,Station,Depth)

# Check autocorrelation of all fluorescence parameters
# Select variables only
fl_data <- ungroup(fl_data)
auto_fl_data <- fl_data %>% select(-c("Date","Station","Depth","Location","Dilution"))

chart.Correlation(auto_fl_data,histogram=TRUE,method=c("spearman"))

# Remove autocorrelated parameters: R2 > 0.8
# A, C, M, N, T, B, HIX, FI, T/M, T/N, T/C, A/T, C/N, A/N
auto_fl_data <- auto_fl_data %>% select("Fmax1","Fmax2","Fmax3","Fmax4","BIX","HIX","T/B","A/C","M/C")

chart.Correlation(auto_fl_data,histogram=TRUE,method=c("pearson"))

# Select fluorescence data to be used in analyses
fl_data_2 <- fl_data %>% select("Date","Station","Depth","Dilution","Fmax1","Fmax2","Fmax3","Fmax4","BIX","HIX","T/B","A/C","M/C")

# Also need to combine SD data for eems and parafac
sd_parafac <- sd_parafac %>% select(-c("i","SampleName","Reservoir","Rep","Location"))
sd_parafac <- sd_parafac[complete.cases(sd_parafac),]

sd_eems <- sd_eems %>% select(-c("Sample Name","Reservoir","Date analyzed","Rep","Dilution","Location"))
sd_eems <- sd_eems[complete.cases(sd_eems),]

sd_fl_data <- full_join(sd_parafac,sd_eems,by=c("Date","Station","Depth"))
sd_fl_data <- sd_fl_data %>% select("Date","Station","Depth","Fmax1","Fmax2","Fmax3","Fmax4","B","BIX","HIX","A/T","A/C","M/C")

# Export out Fl data
write_csv(fl_data_2,'C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20201105_FluorescentData_QAQC.csv')

# Select 2019 data (aka: remove data from March 2020 RC day)
fl_data_2 <- fl_data_2 %>% filter(Date > as.POSIXct("2018-12-31") & Date < as.POSIXct("2020-01-01")) %>% 
  select(-c("Dilution"))

################################ Load in GHG data ##################################
fcrghg <- read.csv("./Data/ghg.csv", header=T) %>%
  select(DateTime:co2_umolL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1|Depth_m == 5.0|Depth_m == 9.0|Depth_m == 100.0|Depth_m == 200.0)

m_ghg <- fcrghg %>% group_by(DateTime,Depth_m) %>% summarise_all(funs(mean(.,na.rm=TRUE)))
sd_ghg <- fcrghg %>% group_by(DateTime,Depth_m) %>% summarise_all(funs(sd(.,na.rm=TRUE)))

# Plot for AGU: Station 50 epi and hypo only
m_ghg_50 <- m_ghg %>% 
  filter(Depth_m == "0.1"|Depth_m == "9") %>% 
  filter(DateTime > as.POSIXct("2019-05-20") & DateTime < as.POSIXct("2019-11-21"))

sd_ghg_50 <- sd_ghg %>% 
  filter(Depth_m == "0.1"|Depth_m == "9") %>% 
  filter(DateTime > as.POSIXct("2019-05-20") & DateTime < as.POSIXct("2019-11-21"))


ch4 <- ggplot(m_ghg_50,mapping=aes(x=DateTime,y=ch4_umolL,color=(as.factor(Depth_m))))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(sd_ghg_50,mapping=aes(ymin=m_ghg_50$ch4_umolL-ch4_umolL,ymax=m_ghg_50$ch4_umolL+ch4_umolL))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab(expression(paste("CH"[4]*" (", mu,"mol L"^-1*")")))+
  scale_color_manual(breaks=c('0.1','9'),values=c("#7EBDC2","#393E41"),labels=c("Epi","Hypo"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=125,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

co2 <- ggplot(m_ghg_50,mapping=aes(x=DateTime,y=co2_umolL,color=(as.factor(Depth_m))))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(sd_ghg_50,mapping=aes(ymin=m_ghg_50$co2_umolL-co2_umolL,ymax=m_ghg_50$co2_umolL+co2_umolL))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab(expression(paste("CO"[2]*" (", mu,"mol L"^-1*")")))+
  scale_color_manual(breaks=c('0.1','9'),values=c("#7EBDC2","#393E41"),labels=c("Epi","Hypo"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1000,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

m_ghg_50_epi <- m_ghg_50 %>% 
  filter(Depth_m == "0.1")
sd_ghg_50_epi <- sd_ghg_50 %>% 
  filter(Depth_m == "0.1")

# Plot CO2 epi ONLY
co2_surf <- ggplot(m_ghg_50_epi,mapping=aes(x=DateTime,y=co2_umolL,color=(as.factor(Depth_m))))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(sd_ghg_50_epi,mapping=aes(ymin=m_ghg_50_epi$co2_umolL-co2_umolL,ymax=m_ghg_50_epi$co2_umolL+co2_umolL))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab(expression(paste("CO"[2]*" (", mu,"mol L"^-1*")")))+
  scale_color_manual(breaks=c('0.1'),values=c("#7EBDC2"),labels=c("Epi"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=300,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Plot Ch4 epi ONLY
ch4_surf <- ggplot(m_ghg_50_epi,mapping=aes(x=DateTime,y=ch4_umolL,color=(as.factor(Depth_m))))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(sd_ghg_50_epi,mapping=aes(ymin=m_ghg_50_epi$ch4_umolL-ch4_umolL,ymax=m_ghg_50_epi$ch4_umolL+ch4_umolL))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab(expression(paste("CH"[4]*" (", mu,"mol L"^-1*")")))+
  scale_color_manual(breaks=c('0.1'),values=c("#7EBDC2"),labels=c("Epi"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=1.5,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(co2,co2_surf,ch4,ch4_surf,common.legend = TRUE)

# Collate data for future use in AR model

m_ghg <- m_ghg %>% select(-c("Reservoir"))
m_ghg <- m_ghg %>% mutate(Station = ifelse(Depth_m < 10.0, 50, round(Depth_m)),
                          Depth = ifelse(Depth_m > 10.0, 0.1, Depth_m))

m_ghg <- m_ghg %>% select(-c("Depth_m","Rep"))
m_ghg <- m_ghg %>% relocate("Station",.after="DateTime") %>% relocate("Depth",.after="Station") %>% 
  rename(Date = DateTime)

m_fl_ghg <- full_join(fl_data_2,m_ghg,by=c("Date","Station","Depth"))
m_fl_ghg <- m_fl_ghg %>% arrange(Date,Station,Depth)

# Export out GHG data for the full time period
write_csv(m_ghg,'C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20201103_GHG_Mean.csv')

# Collate SD data
sd_ghg <- sd_ghg %>% select(-c("Reservoir","Rep"))
sd_ghg <- sd_ghg[complete.cases(sd_ghg),]

sd_ghg <- sd_ghg %>% mutate(Station = ifelse(Depth_m < 10.0, 50, round(Depth_m)),
                          Depth = ifelse(Depth_m > 10.0, 0.1, Depth_m))
sd_ghg <- sd_ghg %>% select(-c("Depth_m"))

sd_ghg <- sd_ghg %>% relocate("Station",.after="DateTime") %>% relocate("Depth",.after="Station") %>% 
  rename(Date = DateTime)

# Combine FL and GHG SD data
sd_fl_ghg <- full_join(sd_fl_data,sd_ghg,by=c("Date","Station","Depth"))
sd_fl_ghg <- sd_fl_ghg %>% arrange(Date,Station,Depth)

########################### Load in DOC/Nutrient data ##########################
fcrchem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DN_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  filter(Site == 50 | Site == 100 | Site ==200) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1|Depth_m == 5.0|Depth_m == 9.0)

fcrchem <- fcrchem %>% select(-c("Reservoir")) %>% 
  rename(Date = DateTime, Depth = Depth_m, Station = Site)

# Plot for AGU
fcrchem_50 <- fcrchem %>% 
  filter(Station == "50") %>% 
  filter(Depth == "0.1"|Depth == "9") %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21"))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

fcrchem_50 <- completeFun(fcrchem_50,"DOC_mgL")

# Plot
doc <- ggplot(fcrchem_50,mapping=aes(x=Date,y=DOC_mgL,color=(as.factor(Depth))))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  scale_color_manual(breaks=c('0.1','9'),values=c("#7EBDC2","#393E41"),labels=c("Epi","Hypo"))+
  annotate("text",x=c(as.POSIXct("2019-05-23"),as.POSIXct("2019-06-10"),as.POSIXct("2019-06-27"),
                      as.POSIXct("2019-07-14"),as.POSIXct("2019-07-27"),as.POSIXct("2019-08-12"),
                      as.POSIXct("2019-08-27"),as.POSIXct("2019-09-12"),as.POSIXct("2019-11-05")),
           y=8,label=c("Off","On","Off","On","Off","On","Off","On","Turnover"),size=5.5)+
  ylim(0,8)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(ch4,doc,co2,bix,ncol=2,nrow=2,common.legend = TRUE)

ggarrange(ch4,co2,nrow=2,ncol=1,common.legend = TRUE)

ggarrange(doc,bix,nrow=2,ncol=1,common.legend = TRUE)

m_fl_ghg_chem <- full_join(m_fl_ghg,fcrchem,by=c("Date","Station","Depth"))
m_fl_ghg_chem <- m_fl_ghg_chem %>% arrange(Date,Station,Depth)

############################## Temp and DO data ##############################
# Load in merged CTD and YSI data (see Cast_Data)
casts <- read_csv("./Data/YSICTD_Merge.csv")

# Then average around 0.1 m, 5.0 m, and 9.0 m
ctd_1 <- casts %>%  
  filter(depth>=0 & depth<0.2) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(Depth=0.1)

ctd_2 <- casts %>% 
  filter(depth>=4.9 & depth<5.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(Depth=5.0)

ctd_3 <- casts %>% 
  filter(depth>=8.9 & depth<9.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(Depth=9.0)

ctd_all <- rbind(ctd_1,ctd_2,ctd_3)
ctd_all <- ctd_all %>% select(-c(depth)) %>% mutate(Station = 50) %>% rename(Date = time) %>% 
  arrange(Date,Station,Depth)

ctd_all <- ctd_all %>% relocate("Station",.after = "Date") %>% relocate("Depth",.after = "Station")

data_all <- full_join(m_fl_ghg_chem,ctd_all,by=c("Date","Station","Depth"))
data_all <- data_all %>% arrange(Date,Station,Depth)

################################ Add Flora data ##########################################
flora <- read.csv("./Data/flora.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Site == 50)

flora_1 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=0 & Depth_m<0.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=0.1)

flora_2 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=4.5 & Depth_m<5.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=5.0)

flora_3 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=8.5 & Depth_m<9.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=9.0)

flora_all <- rbind(flora_1,flora_2,flora_3)

flora_all <- flora_all %>% select(-c(Depth_m)) %>% mutate(Station = 50) %>% rename(Date = DateTime, Flora_ugL = TotalConc_ugL) %>% 
  arrange(Date,Station,Depth)

flora_all <- flora_all %>% relocate("Station",.after = "Date") %>% relocate("Depth",.after = "Station")

data_all <- full_join(data_all,flora_all,by=c("Date","Station","Depth"))
data_all <- data_all %>% arrange(Date,Station,Depth)

################################### Inflow data ###################################
inflow <- read_csv('./Data/Inflow.csv') %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d",tz='EST'))) %>% 
  filter(DateTime > as.POSIXct("2019-01-01") & DateTime < as.POSIXct("2020-01-01")) %>% 
  select(Site,DateTime,WVWA_Flow_cms) %>% 
  rename(Station = Site, Date = DateTime, Flow_cms = WVWA_Flow_cms)

inflow <- inflow %>%  group_by(Date) %>% summarise_all(funs(mean(.,na.rm=TRUE)))

# Load in RC inflow (for site 200)
rc_inflow <- read.csv("./Data/RC_Inflow.csv") %>% 
  filter(Reservoir=="FCR" & Site=="200") %>% select(Date,Site,Flow_cms) %>% 
  rename(Station = Site)

all_inflow <- rbind(inflow,rc_inflow)

all_inflow <- all_inflow %>% arrange(Date) %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%Y-%m-%d",tz='EST')))

data_all <- left_join(data_all,all_inflow,by=c("Date","Station"))

############################ Met data - Daily Precip; Max Solar radiation ########################
# Download met data from EDI
inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/389/4/c1db8816742823eba86696b29f106d0f"
infile1 <- paste0(getwd(),"/Data/MetData.csv")
download.file(inUrl1,infile1,method="curl")

met <- read_csv('./Data/MetData.csv') %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d",tz='EST'))) %>% 
  filter(DateTime > as.POSIXct("2019-01-01") & DateTime < as.POSIXct("2020-01-01"))

met <- met %>% select(Site,DateTime,Rain_Total_mm,ShortwaveRadiationUp_Average_W_m2,ShortwaveRadiationDown_Average_W_m2)

met_rain <- met %>% group_by(DateTime) %>% summarise_all(funs(sum(.,na.rm=TRUE))) %>% select(DateTime,Rain_Total_mm)

met_sw <- met %>% group_by(DateTime) %>% summarise_all(funs(max(.,na.rm=TRUE))) %>% select(DateTime,ShortwaveRadiationUp_Average_W_m2)

met_all <- left_join(met_rain,met_sw,by=c("DateTime"))

met_all <- met_all %>% mutate(Station = 50) %>% mutate(Depth = 0.1) %>% rename(Date = DateTime)

data_all <- left_join(data_all,met_all,by=c("Date","Station","Depth"))

######################## Export out data #############################
write_csv(data_all,"./Data/20201105_All_Data.csv")
