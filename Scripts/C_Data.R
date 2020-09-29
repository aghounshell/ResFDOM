### Script to look at 2019 data
### DOC, GHG, and other environmental parameters
### For use with EEMs/Abs data
### A Hounshell, 28 Sep 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load data
# DOC data from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/6/2b3dc84ae6b12d10bd5485f1c300af13" 
infile1 <- paste0(getwd(),"/Data/chem.csv")
download.file(inUrl1,infile1,method="curl")

fcrchem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

# Just curious...
ggplot(fcrchem,aes(DateTime,DOC_mgL,color=as.factor(Depth_m)))+
  geom_line()

# Select data for 2019
fcrchem_19 <- fcrchem %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Site == 50) %>% 
  filter(Depth_m == 0.1 | Depth_m == 5.0 | Depth_m == 9.0)

fcrdoc_19 <- fcrchem_19 %>% select(Reservoir:Depth_m,DOC_mgL)
fcrdoc_19 <- na.omit(fcrdoc_19)

# Select Inflow and Wetlands data
fcrchem_flow <- fcrchem %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Site == 100 | Site == 200)

fcrdoc_flow <- fcrchem_flow %>% 
  select(Reservoir:Depth_m,DOC_mgL)

fcrdoc_flow <- na.omit(fcrdoc_flow)

# Plot
s50 <- ggplot(fcrdoc_19,aes(DateTime,DOC_mgL,color=as.factor(Depth_m)))+
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
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  ylim(0,7)+
  scale_color_manual(breaks=c('0.1','5','9'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

inf <- ggplot(fcrdoc_flow,mapping=aes(DateTime,DOC_mgL,color=as.factor(Site)))+
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
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  ylim(0,7)+
  scale_color_manual(breaks=c('100','200'),values=c("#F0B670","#FE5F55"))+
  theme_classic(base_size=15)
  
ggarrange(s50,inf,ncol=1,nrow=2)

# GHG data from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/2/38d72673295864956cccd6bbba99a1a3" 
infile1 <- paste0(getwd(),"/Data/ghg.csv")
download.file(inUrl1,infile1,method="curl")

fcrghg <- read.csv("./Data/ghg.csv", header=T) %>%
  select(DateTime:co2_umolL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST")))

fcrghg_surf <- fcrghg %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1) %>% group_by(DateTime) %>% summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  arrange(DateTime) %>% mutate(grouping="FCR_epi")
  
fcrghg_5 <- fcrghg %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 5) %>% group_by(DateTime) %>% summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  arrange(DateTime) %>% mutate(grouping="FCR_meta")

fcrghg_9 <- fcrghg %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 9) %>% group_by(DateTime) %>% summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  arrange(DateTime) %>% mutate(grouping="FCR_hypo")

fcrghg_100 <- fcrghg %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 100) %>% group_by(DateTime) %>% summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  arrange(DateTime) %>% mutate(grouping="FCR_Inf")

fcrghg_200 <- fcrghg %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 200) %>% group_by(DateTime) %>% summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  arrange(DateTime) %>% mutate(grouping="FCR_wet")

# Plot
s50 <- ggplot()+
  geom_line(fcrghg_surf,mapping=aes(DateTime,ch4_umolL,color="epi"),size=1)+
  geom_point(fcrghg_surf,mapping=aes(DateTime,ch4_umolL,color="epi"),size=2)+
  geom_line(fcrghg_5,mapping=aes(DateTime,ch4_umolL,color="meta"),size=1)+
  geom_point(fcrghg_5,mapping=aes(DateTime,ch4_umolL,color="meta"),size=2)+
  geom_line(fcrghg_9,mapping=aes(DateTime,ch4_umolL,color="hypo"),size=1)+
  geom_point(fcrghg_9,mapping=aes(DateTime,ch4_umolL,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

inf <- ggplot()+
  geom_line(fcrghg_100,mapping=aes(DateTime,ch4_umolL,color="100"),size=1)+
  geom_point(fcrghg_100,mapping=aes(DateTime,ch4_umolL,color="100"),size=2)+
  geom_line(fcrghg_200,mapping=aes(DateTime,ch4_umolL,color="200"),size=1)+
  geom_point(fcrghg_200,mapping=aes(DateTime,ch4_umolL,color="200"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('100','200'),values=c("#F0B670","#FE5F55"))+
  theme_classic(base_size=15)

ggarrange(s50,inf,ncol=1,nrow=2)

# Plot CO2
s50 <- ggplot()+
  geom_line(fcrghg_surf,mapping=aes(DateTime,co2_umolL,color="epi"),size=1)+
  geom_point(fcrghg_surf,mapping=aes(DateTime,co2_umolL,color="epi"),size=2)+
  geom_line(fcrghg_5,mapping=aes(DateTime,co2_umolL,color="meta"),size=1)+
  geom_point(fcrghg_5,mapping=aes(DateTime,co2_umolL,color="meta"),size=2)+
  geom_line(fcrghg_9,mapping=aes(DateTime,co2_umolL,color="hypo"),size=1)+
  geom_point(fcrghg_9,mapping=aes(DateTime,co2_umolL,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

inf <- ggplot()+
  geom_line(fcrghg_100,mapping=aes(DateTime,co2_umolL,color="100"),size=1)+
  geom_point(fcrghg_100,mapping=aes(DateTime,co2_umolL,color="100"),size=2)+
  geom_line(fcrghg_200,mapping=aes(DateTime,co2_umolL,color="200"),size=1)+
  geom_point(fcrghg_200,mapping=aes(DateTime,co2_umolL,color="200"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('100','200'),values=c("#F0B670","#FE5F55"))+
  theme_classic(base_size=15)

ggarrange(s50,inf,ncol=1,nrow=2)

###########################################3

# Let's look at other environmental parameters as well!
# Start with Chla data from Flora
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/4/e7e3e6e513985a602d9a5f22687d4efc" 
infile1 <- paste0(getwd(),"/Data/flora.csv")
download.file(inUrl1,infile1,method="curl")

flora <- read.csv("./Data/flora.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Site == 50)

flora_1 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=0 & Depth_m<0.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(grouping="epi")

flora_2 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=4.5 & Depth_m<5.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(grouping="meta")

flora_3 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=8.5 & Depth_m<9.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(grouping="hypo")

# Plot
ggplot()+
  geom_line(flora_1,mapping=aes(DateTime,TotalConc_ugL,color="epi"),size=1)+
  geom_point(flora_1,mapping=aes(DateTime,TotalConc_ugL,color="epi"),size=2)+
  geom_line(flora_2,mapping=aes(DateTime,TotalConc_ugL,color="meta"),size=1)+
  geom_point(flora_2,mapping=aes(DateTime,TotalConc_ugL,color="meta"),size=2)+
  geom_line(flora_3,mapping=aes(DateTime,TotalConc_ugL,color="hypo"),size=1)+
  geom_point(flora_3,mapping=aes(DateTime,TotalConc_ugL,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

