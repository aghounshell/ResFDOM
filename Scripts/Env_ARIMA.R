### Script to look at environmental variables and conduct ARIMA modeling between 
### epi and hypo DOC concentrations and various limno. and met. parameters

### 18 Aug 2022, A. Hounshell

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,zoo,scales,plyr,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

###############################################################################

## Load in a format DOC concentration data from 2017-2021
chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DOC_mgL <= 15) %>% 
  drop_na(DOC_mgL)

###############################################################################

## Load in thermocline data - obtained from CTD and YSI casts
## Calculated using LakeAnalyzer
thermo <- read.csv("./Data/rev_FCR_results_LA.csv") %>% 
  select(datetime,thermo.depth) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = datetime)

###############################################################################

## Load in CTD + YSI data - temp, Sal, DO
## From merged spreadsheet in: LakeAnalyzer_thermo.R
casts <- read.csv("./Data/merged_YSI_CTD.csv") %>% 
  filter(depth >= 0) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = time)

## Calculate Epi and Hypo V.W. parameters using thermocline data
## Following Eco_DOC.R

# Create vector of different volumes for each depth: based on L&O-L 2020 paper
vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0), "Vol_m3" = c(138486.51,89053.28,59619.35,40197.90,13943.82,14038.52,1954.71))

thermo <- thermo %>% 
  mutate(thermo.depth = round(thermo.depth,digits=1))

# Use thermocline depth to calculate epi and hypo V.W. parameters (DO_mgL, Cond_uScm, Chla, Turb)
limno_data <- left_join(casts,thermo,by="DateTime")

limno_data <- limno_data %>% 
  mutate(depth = round(depth, digits=1)) %>% 
  mutate(layer = ifelse(is.na(thermo.depth), "Total",
                        ifelse(thermo.depth >= depth, "Epi",
                               ifelse(thermo.depth <= depth, "Hypo", NA))))

casts_depth <- limno_data %>% 
  group_by(DateTime,layer) %>% 
  dplyr::summarise(mean_Temp_C = mean(Temp_C, na.rm=TRUE),
                sd_Temp_C = sd(Temp_C, na.rm=TRUE),
                mean_DO_mgL = mean(DO_mgL, na.rm=TRUE),
                sd_DO_mgL = sd(DO_mgL, na.rm=TRUE),
                mean_Cond_uScm = mean(Cond_uScm, na.rm=TRUE),
                sd_Cond_uScm = sd(Cond_uScm, na.rm=TRUE),
                mean_Chla_ugL = mean(Chla_ugL, na.rm=TRUE),
                sd_Chla_ugL = sd(Chla_ugL, na.rm=TRUE),
                mean_Turb_NTU = mean(Turb_NTU, na.rm=TRUE),
                sd_Turb_NTU = sd(Turb_NTU, na.rm=TRUE)) %>% 
  filter(layer == "Epi" | layer == "Hypo")

## Plot data by epi and hypo
ggplot(casts_depth,mapping=aes(x=DateTime,y=mean_Temp_C))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(ymin=mean_Temp_C-sd_Temp_C,ymax=mean_Temp_C+sd_Temp_C,fill=layer),alpha=0.5)+
  geom_line(mapping=aes(color=layer),size=1)+
  geom_point(mapping=aes(color=layer),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Temp~(C^o)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

###############################################################################

## Load in Flora data - Chla and community analysis?

###############################################################################

## Load in inflow

###############################################################################

## Load in met data - shortwave radiation and rainfall

###############################################################################