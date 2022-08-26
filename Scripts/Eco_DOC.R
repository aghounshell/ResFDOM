### Updated script to include 'whole-ecosystem' DOC processing (epi and hypo)
### Following comments from PCH and CCC
### Script is drawing heavily from DOC_DO.R

### 11 Aug 2021, A Hounshell
### Revising and updating model and visualizations, 15 Aug 2022

###############################################################################

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,rMR,lme4,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics,truncnorm)

###############################################################################

## Load in various datastreams for box model
# Import Lake Analyzer thermocline results to determine median thermocline depth
la_results <- read.csv("./Data/rev_FCR_results_LA.csv") %>% 
  select(datetime,thermo.depth) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d", tz="EST"))) %>% 
  rename(DateTime = datetime)

### Load in Inflow data ----
# Weir discharge/temperature
# Downloaded on 12 Aug 22 for 2015-2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/8/cc045f9fe32501138d5f4e1e7f40d492"
#infile1 <- paste0(getwd(),"/Data/Inflow_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/Inflow_2013_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:VT_Temp_C)

# Plot WVWA vs. VT inflow
ggplot(inflow,aes(x=WVWA_Flow_cms,y=VT_Flow_cms))+
  geom_point()

# Find relationship between WVWA and VT flow
# Use to fill in any missing WVWA values with VT flow (where possible)
flow_lm <- lm(WVWA_Flow_cms ~ VT_Flow_cms, data = inflow)

inflow <- inflow %>% 
  mutate(WVWA_Flow_cms = ifelse(is.na(WVWA_Flow_cms), 0.8917528*VT_Flow_cms-0.0009302, WVWA_Flow_cms)) %>% 
  mutate(flow_diff = abs(VT_Flow_cms - WVWA_Flow_cms))

# Average inflow by day
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_at(vars("WVWA_Flow_cms"),funs(mean(.,na.rm=TRUE),sd)) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

inflow_daily_full <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
inflow_daily_full <- inflow_daily_full %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
inflow_daily_full <- left_join(inflow_daily_full, inflow_daily,by="DateTime")

# Calculate total variance - daily sd + difference between WVWA and VT inflow (0.002 cms)
inflow_daily_full <- inflow_daily_full %>% 
  mutate(total_sd = sqrt((sd^2)+((mean(inflow$flow_diff,na.rm=TRUE))^2)))

inflow_daily_full <- inflow_daily_full %>% 
  mutate(total_sd = ifelse(is.na(total_sd),mean(inflow_daily_full$total_sd,na.rm=TRUE),total_sd))

# Create 3D array of random variables from mean and total_sd
# Where each array contains the 2557 daily time-steps for inflow
inflow_model_input <- array(data = NA, dim=c(1000,2557,1))

for (j in 1:1000){
  for (i in 1:2557){
    inflow_model_input[j,i,1] <- rtruncnorm(1,a=0,mean=inflow_daily_full$mean[i],sd=inflow_daily_full$sd[i]/2)
  }
}

### Load DOC data ----
# Updated: 19 Apr 2021 with 2020 data!
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98" 
#infile1 <- paste0(getwd(),"/Data/chem.csv")
#download.file(inUrl1,infile1,method="curl")

chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

chem_100 <- chem %>% 
  filter(Site == 100) %>% 
  drop_na(DOC_mgL)

chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DOC_mgL <= 15) %>% 
  drop_na(DOC_mgL)

## Merge and plot measured inflow vs. measured DOC concentration
doc_100 <- left_join(chem_100,inflow_daily,by="DateTime")

docvinflow <- doc_100 %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  ggplot(mapping=aes(x=mean,y=DOC_mgL))+
  geom_point(size=1.5)+
  geom_smooth(method='lm')+ 
  xlab(expression(paste("Inflow (m"^3*"s"^-1*")")))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)

ggsave("./Fig_Output/SI_DOCvInflow.png",docvinflow,dpi=800,width=5,height=4)

doc_lm <- lm(DOC_mgL ~ mean, data = doc_100)
summary(doc_lm)

# Plot DOC concentrations with depth for the study period - separated by year
# For supplementary information
doc_2017 <- ggplot()+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=2)+
  geom_line(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=0.8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2017-12-31"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2017")+
  labs(color="Depth (m)")+
  theme_classic(base_size = 15)

doc_2018 <- ggplot()+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=2)+
  geom_line(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=0.8)+
  xlim(as.POSIXct("2018-01-01"),as.POSIXct("2018-12-31"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2018")+
  labs(color="Depth (m)")+
  theme_classic(base_size = 15)

doc_2019 <- ggplot()+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=2)+
  geom_line(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=0.8)+
  xlim(as.POSIXct("2019-01-01"),as.POSIXct("2019-12-31"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2019")+
  labs(color="Depth (m)")+
  theme_classic(base_size = 15)

doc_2020 <- ggplot()+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=2)+
  geom_line(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=0.8)+
  xlim(as.POSIXct("2020-01-01"),as.POSIXct("2020-12-31"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2020")+
  labs(color="Depth (m)")+
  theme_classic(base_size = 15)

doc_2021 <- ggplot()+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=2)+
  geom_line(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)),size=0.8)+
  geom_point(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=2)+
  geom_line(chem_100,mapping=aes(x=DateTime,y=DOC_mgL,color="Inflow"),size=0.8)+
  xlim(as.POSIXct("2021-01-01"),as.POSIXct("2021-12-31"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2021")+
  labs(color="Depth (m)")+
  theme_classic(base_size = 15)

doc_conc <- ggarrange(doc_2017,doc_2018,doc_2019,doc_2020,doc_2021,ncol=2,nrow=3,common.legend = TRUE)

ggsave("./Fig_Output/SI_DOC_Concentration.png",doc_conc,dpi=800,width=13,height=9)

## Plot thermocline depth through time and discharge at the primary inflow
thermo <- ggplot(la_results,mapping=aes(x=DateTime,y=-thermo.depth))+
  geom_hline(yintercept = -0.1,linetype="dashed",color="grey")+
  geom_hline(yintercept = -1.6, linetype="dashed",color="grey")+
  geom_hline(yintercept = -3.8, linetype="dashed",color="grey")+
  geom_hline(yintercept = -5, linetype="dashed",color="grey")+
  geom_hline(yintercept = -6.2, linetype="dashed",color="grey")+
  geom_hline(yintercept = -8, linetype="dashed",color="grey")+
  geom_hline(yintercept = -9, linetype="dashed",color="grey")+
  geom_point(size=2)+
  geom_line(size=0.8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("")+
  ylab("Thermocline depth (m)")+
  theme_classic(base_size = 15)

qcharge <- ggplot(inflow_daily_full,mapping=aes(x=DateTime,y=mean))+
  geom_line(inflow_daily,mapping=aes(x=DateTime,y=mean),size=0.8)+
  geom_ribbon(mapping=aes(ymin=mean-total_sd,ymax=mean+total_sd),alpha=0.6)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("")+
  ylab(expression(paste("Discharge (m"^3*" s"^-1*")")))+
  theme_classic(base_size = 15)

ggarrange(thermo,qcharge,ncol=1,nrow=2)

ggsave("./Fig_Output/SI_ThermoDepth_Q.png",dpi=800,width=9,height=8)

# Create vector of different volumes for each depth: based on L&O-L 2020 paper
vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0), "Vol_m3" = c(138486.51,89053.28,59619.35,40197.90,13943.82,14038.52,1954.71))

### Calculate DOC processing for Epi and Hypo -----
# First, separate by depth and convert to mass
doc_box <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_") %>% 
  mutate(DOC_0.1 = na.fill(na.approx(DOC_0.1,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_1.6 = na.fill(na.approx(DOC_1.6,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_3.8 = na.fill(na.approx(DOC_3.8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_5 = na.fill(na.approx(DOC_5,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_6.2 = na.fill(na.approx(DOC_6.2,na.rm=FALSE),"extend")) %>%   
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend"))

doc_box_full <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
doc_box_full <- doc_box_full %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
doc_box_full <- left_join(doc_box_full,doc_box,by="DateTime")

doc_box_full <- doc_box_full %>% 
  mutate(DOC_0.1 = na.fill(na.approx(DOC_0.1,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_1.6 = na.fill(na.approx(DOC_1.6,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_3.8 = na.fill(na.approx(DOC_3.8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_5 = na.fill(na.approx(DOC_5,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_6.2 = na.fill(na.approx(DOC_6.2,na.rm=FALSE),"extend")) %>%   
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend"))

## Create 'boot-strapped' values from [DOC] and MDL
# Create 3D array of random variables from the measured DOC and the MDL (0.11)
# For each timestep and each depth
doc_lake_model_input <- array(data = NA, dim=c(1000,2557,7))

for (j in 1:1000){
  for (i in 1:2557){
    for (k in 1:7){
      doc_lake_model_input[j,i,k] <- rtruncnorm(1,a=0,mean=doc_box_full[i,k+1],sd=0.11/2)
    }
  }
}

## Create boot-strapped parameters for volume - assuming volume +/-10%
vol_model_input <- array(data = NA, dim=c(1000,7,1))

for (j in 1:1000){
  for (i in 1:7){
    vol_model_input[j,i,1] <- rtruncnorm(1,a=0,mean=vol_depths$Vol_m3[i],sd=vol_depths$Vol_m3[i]*0.05)
  }
}

## Create boot-strapped parameters for DOC inflow - using LOQ for the SD
doc_inflow_full <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
doc_inflow_full <- doc_inflow_full %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
doc_inflow_full <- left_join(doc_inflow_full,chem_100,by="DateTime")

doc_inflow_full <- doc_inflow_full %>% 
  select(DateTime,DOC_mgL) %>% 
  mutate(DOC_mgL = na.fill(na.approx(DOC_mgL,na.rm=FALSE),"extend"))

doc_inflow_input <- array(data = NA, dim=c(1000,2558,1))

for (j in 1:1000){
  for (i in 1:2558){
    doc_inflow_input[j,i,1] <- rtruncnorm(1,a=0,mean=doc_inflow_full$DOC_mgL[i],sd=0.11/2)
  }
}

###############################################################################

## Determine location of the epi and hypo for each time point
# Add in thermocline depth information
thermocline_depth <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
thermocline_depth <- thermocline_depth %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
thermocline_depth <- left_join(thermocline_depth,la_results,by="DateTime")

thermocline_depth <- thermocline_depth %>% 
  select(DateTime,thermo.depth) %>% 
  mutate(thermo.depth = na.fill(na.approx(thermo.depth, na.rm=FALSE),"extend"))

# Use thermocline to find the location of the epi and hypo
thermocline_depth <- thermocline_depth %>% 
  mutate(thermo.depth = round(thermo.depth,digits=1)) %>% 
  mutate(epi_bottom_depth_m = ifelse(thermo.depth > 9.0, 9.0,
                                     ifelse(thermo.depth > 7.0, 8.0,
                                            ifelse(thermo.depth > 6.0, 6.2,
                                                   ifelse(thermo.depth > 4.4, 5.0,
                                                          ifelse(thermo.depth > 3.0, 3.8,
                                                                 ifelse(thermo.depth > 1.6, 1.6,
                                                                        ifelse(thermo.depth > 0.0, 0.1, NA)))))))) %>% 
  mutate(hypo_top_depth_m = ifelse(thermo.depth <= 0.0, 0.1,
                                   ifelse(thermo.depth <= 1.6, 1.6,
                                          ifelse(thermo.depth <= 3.0, 3.8,
                                                 ifelse(thermo.depth <= 4.4, 5.0,
                                                        ifelse(thermo.depth <= 6.0, 6.2,
                                                               ifelse(thermo.depth <= 7.0, 8.0,
                                                                      ifelse(thermo.depth <= 9.0, 9.0, NA))))))))

###############################################################################

## Have bootstrapped inputs for:
# In lake DOC concentrations for all depths (doc_lake_model_input)
# Inflow DOC concentrations (doc_inflow_input)
# Inflow rates (inflow_model_input)
# Volume of layers (vol_model_input)
# Then also have location of the thermocline (and therefore, of the epi and hypo)

## Then calculate processing

### Calculate total mass of DOC at each depth ###
doc_lake_mass <- array(data = NA, dim=c(1000,2557,7)) # Mass of DOC at each depth (DOC_mgL * Vol_m3 = DOC_g)

for (j in 1:1000){
  for (i in 1:2557){
    for (k in 1:7){
      doc_lake_mass[j,i,k] <- vol_model_input[j,k,1]*doc_lake_model_input[j,i,k]
    }
  }
}

## Calculate mass of inflow - for hypo and epi (will need to include parameter estimation at some point)

doc_inflow_mass <- array(data = NA, dim=c(1000,2557,1)) # Mass of DOC inflow (total) - DOC_mgL * m3_s = g/day

for (j in 1:1000){
  for (i in 1:2557){
    doc_inflow_mass[j,i,1] <- doc_inflow_input[j,i,1]*inflow_model_input[j,i,1]*60*60*24
  }
}

### Calculate outflow from epi and hypo ###

# Epi outflow = inflow * [DOC] at 0.1 m = g/d

doc_epi_mass_outflow <-  array(data = NA, dim=c(1000,2557,1)) # Mass of DOC outflow from the epi

for (j in 1:1000){
  for (i in 1:2557){
    doc_epi_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,1]*60*60*24
  }
}

# 'Outflow' to Epi from Hypo: concentration * outflow = mass/day
doc_hypo_mass_outflow <- array(data = NA, dim=c(1000, 2557, 1)) # Mass of DOC outflow from the hypo

for (j in 1:1000){
  for (i in 1:2557){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,2]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,3]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,4]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,5]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,6]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9) {
      doc_hypo_mass_outflow[j,i,1] <- inflow_model_input[j,i,1]*doc_lake_model_input[j,i,7]*60*60*24
    }
  }
}

### Calculate [DOC]/dt for each time point ###

## First calculate mass of epi and hypo at any time point
# Include changes due to the thermocline
doc_epi_mass <- array(data = NA, dim=c(1000,2557,1)) # Mass of epi at each time point

for (j in 1:1000){
  for (i in 1:2557){
    if(thermocline_depth$epi_bottom_depth_m[i] == 0.1){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 1.6){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]+doc_lake_mass[j,i,2]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 3.8){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]+doc_lake_mass[j,i,2]+doc_lake_mass[j,i,3]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 5){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]+doc_lake_mass[j,i,2]+doc_lake_mass[j,i,3]+doc_lake_mass[j,i,4]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 6.2){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]+doc_lake_mass[j,i,2]+doc_lake_mass[j,i,3]+doc_lake_mass[j,i,4]+doc_lake_mass[j,i,5]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 8){
      doc_epi_mass[j,i,1] <- doc_lake_mass[j,i,1]+doc_lake_mass[j,i,2]+doc_lake_mass[j,i,3]+doc_lake_mass[j,i,4]+doc_lake_mass[j,i,5]+doc_lake_mass[j,i,6]
    }
  }
}

doc_hypo_mass <- array(data=NA, dim=c(1000,2557,1)) # Mass of hypo at each time point

for (j in 1:1000){
  for (i in 1:2557){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,2]+doc_lake_mass[j,i,3]+doc_lake_mass[j,i,4]+doc_lake_mass[j,i,5]+doc_lake_mass[j,i,6]+doc_lake_mass[j,i,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,3]+doc_lake_mass[j,i,4]+doc_lake_mass[j,i,5]+doc_lake_mass[j,i,6]+doc_lake_mass[j,i,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,4]+doc_lake_mass[j,i,5]+doc_lake_mass[j,i,6]+doc_lake_mass[j,i,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,5]+doc_lake_mass[j,i,6]+doc_lake_mass[j,i,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,6]+doc_lake_mass[j,i,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9){
      doc_hypo_mass[j,i,1] <- doc_lake_mass[j,i,7]
    }
  }
}

## Then calculate changing volume: 
epi_vol <- array(data=NA, dim=c(1000,2557,1)) # Volume of epi for each time point

for (j in 1:1000){
  for (i in 1:2557){
    if(thermocline_depth$epi_bottom_depth_m[i] == 0.1){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 1.6){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]+vol_model_input[j,2,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 3.8){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]+vol_model_input[j,2,1]+vol_model_input[j,3,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 5){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]+vol_model_input[j,2,1]+vol_model_input[j,3,1]+vol_model_input[j,4,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 6.2){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]+vol_model_input[j,2,1]+vol_model_input[j,3,1]+vol_model_input[j,4,1]+vol_model_input[j,5,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 8){
      epi_vol[j,i,1] <- vol_model_input[j,1,1]+vol_model_input[j,2,1]+vol_model_input[j,3,1]+vol_model_input[j,4,1]+vol_model_input[j,5,1]+vol_model_input[j,6,1]
    }
  }
}

hypo_vol <- array(data=NA, dim=c(1000,2557,1)) # Volume of hypo for each time point

for (j in 1:1000){
  for (i in 1:2557){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6){
      hypo_vol[j,i,1] <- vol_model_input[j,2,1]+vol_model_input[j,3,1]+vol_model_input[j,4,1]+vol_model_input[j,5,1]+vol_model_input[j,6,1]+vol_model_input[j,7,1]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8){
      hypo_vol[j,i,1] <- vol_model_input[j,3,1]+vol_model_input[j,4,1]+vol_model_input[j,5,1]+vol_model_input[j,6,1]+vol_model_input[j,7,1]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5){
      hypo_vol[j,i,1] <- vol_model_input[j,4,1]+vol_model_input[j,5,1]+vol_model_input[j,6,1]+vol_model_input[j,7,1]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2){
      hypo_vol[j,i,1] <- vol_model_input[j,5,1]+vol_model_input[j,6,1]+vol_model_input[j,7,1]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8){
      hypo_vol[j,i,1] <- vol_model_input[j,6,1]+vol_model_input[j,7,1]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9){
      hypo_vol[j,i,1] <- vol_model_input[j,7,1]
    }
  }
}

## Now calculate DOC/dt for epi and hypo for each timepoint
doc_dt_epi <- array(data=NA, dim=c(1000,2557,1)) # Change in DOC/dt for the epi

for (j in 1:1000){
  for (i in 1:2556){
    doc_dt_epi[j,i+1,1] <- (doc_epi_mass[j,i+1,1]-doc_epi_mass[j,i,1])
  }
}

doc_dt_hypo <- array(data=NA, dim=c(1000,2557,1)) # Change in DOC/dt for the hypo

for (j in 1:1000){
  for (i in 1:2556){
    doc_dt_hypo[j,i+1,1] <- (doc_hypo_mass[j,i+1,1]-doc_hypo_mass[j,i,1])
  }
}

### Calculate entrainment from hypo to epi for each time point ###
# Loosely following FCR_DOCModel_edited_19May17 from Carey et al. 2018

# Entrainment
#if Entr is positive, then epi is getting bigger and hypo is getting smaller; 
#if Entr is negative, then hypo is getting bigger and epi is getting smaller

# Double check for calculating changes in thermocline depth
test <- as.data.frame(matrix(data=NA,nrow=2557,ncol=1))
for (i in 2:2557){
  test[i,1] <- thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[i-1]
}

doc_entr <- array(data=NA, dim=c(1000,2557,1)) # Entrainment for each time point

for (j in 1:1000){
  for (i in 2:2557){
    if(thermocline_depth$epi_bottom_depth_m[i] == thermocline_depth$epi_bottom_depth_m[(i-1)]){
      doc_entr[j,i,1] <- 0
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -7.9){
      doc_entr[j,i,1] <- (-sum(doc_lake_mass[j,i,2:6]))
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -4.9){
      doc_entr[j,i,1] <- (-sum(doc_lake_mass[j,i,2:4]))
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -3.4){
      doc_entr[j,i,1] <- (-sum(doc_lake_mass[j,i,3:4]))
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -3.0){
      doc_entr[j,i,1] <- (-sum(doc_lake_mass[j,i,5:6]))
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -2.4){
      doc_entr[j,i,1] <- (-sum(doc_lake_mass[j,i,4:5]))
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -2.2){
      doc_entr[j,i,1] <- -doc_lake_mass[j,i,3]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.8){
      doc_entr[j,i,1] <- -doc_lake_mass[j,i,6]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.5){
      doc_entr[j,i,1] <- -doc_lake_mass[j,i,2]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.2){
      doc_entr[j,i,1] <- -doc_lake_mass[j,i,5]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.2){
      doc_entr[j,i,1] <- doc_lake_mass[j,i,5]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.5){
      doc_entr[j,i,1] <- doc_lake_mass[j,i,2]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.8){
      doc_entr[j,i,1] <- doc_lake_mass[j,i,6]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 2.2){
      doc_entr[j,i,1] <- doc_lake_mass[j,i,3]
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 3.0){
      doc_entr[j,i,1] <- sum(doc_lake_mass[j,i,5:6])
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 3.4){
      doc_entr[j,i,1] <- sum(doc_lake_mass[j,i,3:4])
    } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 6.4){
      doc_entr[j,i,1] <- sum(doc_lake_mass[j,i,3:6])
    }
  }
}

###############################################################################

## Then calculate processing:
# NOTE: Does not include parameter tuning for epi vs. hypo fraction!!

# Calculate processing for epi
doc_epi_process_g <- array(data=NA, dim=c(1000,2557,1)) # DOC epi processing for each time point
doc_epi_process_mgL <- array(data=NA, dim=c(1000,2557,1))

p = 0.74 # Parameter tuning needed here!

for (j in 1:1000){
  for (i in 1:2557){
    doc_epi_process_g[j,i,1] = doc_dt_epi[j,i,1]-(doc_inflow_mass[j,i,1]*p)-(doc_hypo_mass_outflow[j,i,1]*(1-p))+doc_epi_mass_outflow[j,i,1]-doc_entr[j,i,1]
  }
}

for (j in 1:1000){
  for (i in 1:2557){
    doc_epi_process_mgL[j,i,1] = (doc_dt_epi[j,i,1]-(doc_inflow_mass[j,i,1]*p)-(doc_hypo_mass_outflow[j,i,1]*(1-p))+doc_epi_mass_outflow[j,i,1]-doc_entr[j,i,1])/epi_vol[j,i,1]
  }
}

doc_hypo_process_g <- array(data=NA, dim=c(1000,2557,1))
doc_hypo_process_mgL <- array(data=NA, dim=c(1000,2557,1))

for (j in 1:1000){
  for (i in 1:2557){
    doc_hypo_process_g[j,i,1] = doc_dt_hypo[j,i,1]-(doc_inflow_mass[j,i,1]*(1-p)+(doc_hypo_mass_outflow[j,i,1]*(1-p))+doc_entr[j,i,1])
  }
}

for (j in 1:1000){
  for (i in 1:2557){
    doc_hypo_process_mgL[j,i,1] = (doc_dt_hypo[j,i,1]-(doc_inflow_mass[j,i,1]*(1-p)+(doc_hypo_mass_outflow[j,i,1]*(1-p))+doc_entr[j,i,1]))/hypo_vol[j,i,1]
  }
}

### Average across model runs and calculate sd for reach of the various inputs and for processing
doc_inputs_g <- as.data.frame(matrix(data=NA,nrow=2557,ncol=20))

colnames(doc_inputs_g) <- c('mean_doc_inflow_g',
                            'sd_doc_inflow_g',
                            'mean_doc_hypo_outflow_g',
                            'sd_doc_hypo_outflow_g',
                            'mean_doc_dt_hypo_g',
                            'sd_doc_dt_hypo_g',
                            'mean_doc_entr_g',
                            'sd_doc_entr_g',
                            'mean_doc_dt_epi_g',
                            'sd_doc_dt_epi_g',
                            'mean_doc_epi_outflow_g',
                            'sd_doc_epi_outflow_g',
                            'mean_doc_epi_process_g',
                            'sd_doc_epi_process_g',
                            'mean_doc_epi_process_mgL',
                            'sd_doc_epi_process_mgL',
                            'mean_doc_hypo_process_g',
                            'sd_doc_hypo_process_g',
                            'mean_doc_hypo_process_mgL',
                            'sd_doc_hypo_process_mgL')

for (i in 1:2557){
  doc_inputs_g$mean_doc_inflow_g[i] = mean(doc_inflow_mass[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_inflow_g[i] = sd(doc_inflow_mass[,i,],na.rm=TRUE)
  
  doc_inputs_g$mean_doc_hypo_outflow_g[i] = mean(doc_hypo_mass_outflow[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_outflow_g[i] = sd(doc_hypo_mass_outflow[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_dt_hypo_g[i] = mean(doc_dt_hypo[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_dt_hypo_g[i] = sd(doc_dt_hypo[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_entr_g[i] = mean(doc_entr[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_entr_g[i] = sd(doc_entr[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_dt_epi_g[i] = mean(doc_dt_epi[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_dt_epi_g[i] = sd(doc_dt_epi[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_epi_outflow_g[i] = mean(doc_epi_mass_outflow[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_outflow_g[i] = sd(doc_epi_mass_outflow[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_epi_process_g[i] = mean(doc_epi_process_g[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_process_g[i] = sd(doc_epi_process_g[,i,],na.rm=TRUE)
  
  doc_inputs_g$mean_doc_epi_process_mgL[i] = mean(doc_epi_process_mgL[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_process_mgL[i] = sd(doc_epi_process_mgL[,i,],na.rm=TRUE)

  doc_inputs_g$mean_doc_hypo_process_g[i] = mean(doc_hypo_process_g[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_process_g[i] = sd(doc_hypo_process_g[,i,],na.rm=TRUE)
  
  doc_inputs_g$mean_doc_hypo_process_mgL[i] = mean(doc_hypo_process_mgL[,i,],na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_process_mgL[i] = sd(doc_hypo_process_mgL[,i,],na.rm=TRUE)
}

DateTime <- doc_box_full %>% 
  select(DateTime)

doc_inputs_g <- cbind(DateTime,doc_inputs_g)

doc_inputs_g <- na.omit(doc_inputs_g)

doc_inputs_g <- doc_inputs_g %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"))

## Plot epi and hypo processing
doc_proc_17 <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-06-01"),as.POSIXct("2017-11-15"))+
  ylim(-3,0.5)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2017")+
  theme_classic(base_size=15)+
  guides(fill=FALSE)+
  theme(legend.title=element_blank())

doc_proc_18 <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-11-15"))+
  ylim(-5,5)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2018")+
  theme_classic(base_size=15)+
  guides(fill=FALSE)+
  theme(legend.title=element_blank())

doc_proc_19 <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-11-15"))+
  ylim(-0.5,2)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2019")+
  theme_classic(base_size=15)+
  guides(fill=FALSE)+
  theme(legend.title=element_blank())

doc_proc_20 <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-11-15"))+
  #ylim(-0.5,2)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2020")+
  theme_classic(base_size=15)+
  guides(fill=FALSE)+
  theme(legend.title=element_blank())

doc_proc_21 <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1.5)+
  xlim(as.POSIXct("2021-06-01"),as.POSIXct("2021-11-15"))+
  ylim(-4,2.5)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("2021")+
  theme_classic(base_size=15)+
  guides(fill=FALSE)+
  theme(legend.title=element_blank())

ggarrange(doc_proc_17,doc_proc_18,doc_proc_19,doc_proc_20,doc_proc_21,ncol=1,nrow=5,common.legend=TRUE,
          labels = c("A.", "B.", "C.", "D.","E."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/EpiHypo_DOC_Processing.jpg",width=9,height=12,units="in",dpi=320)

## Plot hypo inputs
hypo_inflow <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=((mean_doc_inflow_g*0.26/1000)-(sd_doc_inflow_g*0.26/1000)),ymax=((mean_doc_inflow_g*0.26/1000)+(sd_doc_inflow_g*0.26/1000)),fill="Inflow"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.26/1000),color="Inflow"))+
  geom_point(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.26/1000),color="Inflow"))+
  ylab(expression(paste("Hypo Inflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Inflow"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("Inflow"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

hypo_outflow <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_hypo_outflow_g*0.26/1000)-(sd_doc_hypo_outflow_g*0.26/1000),ymax=(mean_doc_hypo_outflow_g*0.26/1000)+(sd_doc_hypo_outflow_g*0.26/1000),fill="Outflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Outflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Outflow"))+
  ylab(expression(paste("Hypo Outflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Outflow"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("Outflow"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

hypo_change <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_dt_hypo_g/1000-sd_doc_dt_hypo_g/1000,ymax=mean_doc_dt_hypo_g/1000+sd_doc_dt_hypo_g/1000,fill="DOC/dt"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_dt_hypo_g/1000,color="DOC/dt"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_dt_hypo_g/1000,color="DOC/dt"))+
  ylab(expression(paste("d[DOC]/dt (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("DOC/dt"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("DOC/dt"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

hypo_entr <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymax=(mean_doc_entr_g-sd_doc_entr_g)/1000,ymin=(mean_doc_entr_g+sd_doc_entr_g)/1000,fill="Entr"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_entr_g/1000,color="Entr"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_entr_g/1000,color="Entr"))+
  ylab(expression(paste("Entr. (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Entr"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("Entr"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

ggarrange(hypo_inflow,hypo_outflow,hypo_change,hypo_entr,nrow=4,ncol=1,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Hypo_Inputs.jpg",width=9,height=12,units="in",dpi=320)

## Calculate model inputs/outputs for Epi
epi_inflow <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=((mean_doc_inflow_g*0.74/1000)-(sd_doc_inflow_g*0.74/1000)),ymax=((mean_doc_inflow_g*0.74/1000)+(sd_doc_inflow_g*0.74/1000)),fill="Inflow"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.74/1000),color="Inflow"))+
  geom_point(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.74/1000),color="Inflow"))+
  ylab(expression(paste("Epi Inflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Inflow"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Inflow"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

epi_inflow_h <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_hypo_outflow_g*0.26/1000)-(sd_doc_hypo_outflow_g*0.26/1000),ymax=(mean_doc_hypo_outflow_g*0.26/1000)+(sd_doc_hypo_outflow_g*0.26/1000),fill="Outflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Outflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Outflow"))+
  ylab(expression(paste("Hypo Outflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Outflow"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Outflow"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

epi_outflow <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_epi_outflow_g/1000)-(sd_doc_epi_outflow_g/1000),ymax=(mean_doc_epi_outflow_g/1000)+(sd_doc_epi_outflow_g/1000),fill="Outflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_outflow_g/1000,color="Outflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_outflow_g/1000,color="Outflow"))+
  ylab(expression(paste("Epi Outflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Outflow"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Outflow"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

epi_change <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_dt_epi_g/1000-sd_doc_dt_epi_g/1000,ymax=mean_doc_dt_epi_g/1000+sd_doc_dt_epi_g/1000,fill="DOC/dt"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_dt_epi_g/1000,color="DOC/dt"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_dt_epi_g/1000,color="DOC/dt"))+
  ylab(expression(paste("d[DOC]/dt (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("DOC/dt"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("DOC/dt"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

epi_entr <- ggplot(doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymax=-(mean_doc_entr_g-sd_doc_entr_g)/1000,ymin=-(mean_doc_entr_g+sd_doc_entr_g)/1000,fill="Entr"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=-mean_doc_entr_g/1000,color="Entr"))+
  geom_point(mapping=aes(x=DateTime,y=-mean_doc_entr_g/1000,color="Entr"))+
  ylab(expression(paste("Entr. (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Entr"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Entr"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size=15)+
  theme(legend.position="none")

ggarrange(epi_inflow,epi_inflow_h,epi_outflow,epi_change,epi_entr,nrow=5,ncol=1,labels = c("A.", "B.", "C.", "D.","E."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Epi_Inputs.jpg",width=9,height=15,units="in",dpi=320)

###############################################################################

### Calculate volume weighted epi and hypo concentration using variable thermocline from model above

doc_wgt <- left_join(doc_box,thermocline_depth,by="DateTime")

doc_wgt <- doc_wgt %>% 
  mutate(doc_epi_mgL = ifelse(epi_bottom_depth_m == 0.1, DOC_0.1*vol_depths$Vol_m3[1]/sum(vol_depths$Vol_m3[1]), 
                              ifelse(epi_bottom_depth_m == 1.6, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2]))/sum(vol_depths$Vol_m3[1:2]),
                                     ifelse(epi_bottom_depth_m == 3.8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                            ifelse(epi_bottom_depth_m == 5.0, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                   ifelse(epi_bottom_depth_m == 6.2, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                          ifelse(epi_bottom_depth_m == 8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA))))))) %>% 
  mutate(doc_hypo_mgL = ifelse(hypo_top_depth_m == 1.6, ((DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]), 
                               ifelse(hypo_top_depth_m == 3.8, ((DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                      ifelse(hypo_top_depth_m == 5, ((DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                             ifelse(hypo_top_depth_m == 6.2, ((DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                    ifelse(hypo_top_depth_m == 8, ((DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                           ifelse(hypo_top_depth_m == 9, (DOC_9*vol_depths$Vol_m3[7])/sum(vol_depths$Vol_m3[7]),NA)))))))

doc_wgt <- full_join(doc_wgt,chem_100, by = "DateTime")

doc_wgt <- doc_wgt %>% 
  select(DateTime, doc_epi_mgL, doc_hypo_mgL, DOC_mgL) %>% 
  arrange(DateTime) %>% 
  pivot_longer(!DateTime, names_to = "Loc", values_to = "DOC_mgL")

doc_wgt <- doc_wgt %>% 
  mutate(Loc = ifelse(Loc == "doc_epi_mgL", "Epi",
                      ifelse(Loc == "doc_hypo_mgL", "Hypo", 
                             ifelse(Loc== "DOC_mgL", "Inflow", NA)))) %>% 
  drop_na()

## Plot vol weighted epi and hypo [DOC] and inflow [DOC]
vw_epi <- doc_wgt %>% 
  filter(Loc == "Epi") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  #geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed",color="darkgrey")+
  #geom_vline(xintercept = as.POSIXct("2016-10-09"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#7EBDC2")+
  geom_point(size=2,color="#7EBDC2")+
  xlab("")+
  ylab(expression(paste("Epi VW DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

vw_hypo <- doc_wgt %>% 
  filter(Loc == "Hypo") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  #geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed",color="darkgrey")+
  #geom_vline(xintercept = as.POSIXct("2016-10-09"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#393E41")+
  geom_point(size=2,color="#393E41")+
  xlab("")+
  ylab(expression(paste("Hypo VW DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

inflow_doc <- doc_wgt %>% 
  filter(Loc == "Inflow") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  #geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed",color="darkgrey")+
  #geom_vline(xintercept = as.POSIXct("2016-10-09"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#F0B670")+
  geom_point(size=2,color="#F0B670")+
  xlab("")+
  ylab(expression(paste("Inflow DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

ggarrange(vw_epi,vw_hypo,inflow_doc,nrow=3,ncol=1,labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig3_VW_DOC.jpg",width=7,height=10,units="in",dpi=320)

doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  ggplot(mapping=aes(x=Loc,y=DOC_mgL,fill=Loc))+
  geom_boxplot(size=0.8,alpha=0.5)+
  scale_fill_manual(breaks=c('Epi','Hypo','Inflow'),values=c("#7EBDC2","#393E41","#F0B670"))+
  xlab("")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

ggsave("./Fig_Output/SI_DOC_Boxplot.png",dpi=800,width=5,height=4)

doc_stats <- doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  group_by(Loc) %>% 
  summarise(min = min(DOC_mgL),
            max = max(DOC_mgL),
            median = median(DOC_mgL),
            sd = sd(DOC_mgL))

##### STOPPED HERE ######

#####################################################################################

# Plot stratification
ggplot(la_results)+
  geom_line(mapping=aes(x=DateTime,y=-SthermD))+
  geom_point(mapping=aes(x=DateTime,y=-SthermD))+
  xlim(as.POSIXct("2015-01-01"),as.POSIXct("2020-12-31"))+
  theme_classic(base_size = 15)

### Load in Flora data for Chla ----
# Load from EDI: downloaded on 10 June 2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/5/a24f4dbe9f0d856688f257547069d0a3" 
#infile1 <- paste0(getwd(),"/Data/FluoroProbe.csv")
#download.file(inUrl1,infile1,method="curl")
flora <- read.csv("./Data/FluoroProbe.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2015-01-01")) %>% 
  filter(Site == 50)

### Create a dataframe for Flora parameters at each sampling depth
depths<- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

#Initialize an empty matrix with the correct number of rows and columns 
temp<-matrix(data=NA, ncol=ncol(flora), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final<-matrix(data=NA, ncol=1, nrow=0)
dates<-unique(flora$DateTime)

#create a function to chose the matching depth closest to our focal depths
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j=dates[i]
  q <- subset(flora, flora$DateTime == j)
  
  layer1 <- q[q[, "Depth_m"] == closest(q$Depth_m,0.1),][1,]
  layer2<- q[q[, "Depth_m"] == closest(q$Depth_m,1.6),][1,]
  layer3<- q[q[, "Depth_m"] == closest(q$Depth_m,3.8),][1,]
  layer4<- q[q[, "Depth_m"] == closest(q$Depth_m,5.0),][1,]
  layer5<- q[q[, "Depth_m"] == closest(q$Depth_m,6.2),][1,]
  layer6<- q[q[, "Depth_m"] == closest(q$Depth_m,8.0),][1,]
  layer7<- q[q[, "Depth_m"] == closest(q$Depth_m,9.0),][1,]
  
  temp<-rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(flora))+1)] <- depths
  colnames(temp)[((ncol(flora))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
flora_depths <- as.data.frame(super_final) %>%
  select(-c(1,Depth_m)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

flora_depths$TotalConc_ugL <- as.numeric(flora_depths$TotalConc_ugL)
flora_depths$new_depth <- as.numeric(flora_depths$new_depth)

# Format long
flora_depths_2 <- flora_depths %>% 
  select(DateTime,TotalConc_ugL,new_depth) %>% 
  relocate(new_depth,.before=TotalConc_ugL) %>% 
  pivot_wider(names_from = new_depth, values_from = TotalConc_ugL, values_fil = NA, values_fn = mean, names_prefix = "Chla_") %>% 
  mutate(Chla_0.1_ug = Chla_0.1*vol_depths$Vol_m3[1]*1000) %>% 
  mutate(Chla_1.6_ug = Chla_1.6*vol_depths$Vol_m3[2]*1000) %>% 
  mutate(Chla_3.8_ug = Chla_3.8*vol_depths$Vol_m3[3]*1000) %>% 
  mutate(Chla_5_ug = Chla_5*vol_depths$Vol_m3[4]*1000) %>% 
  mutate(Chla_6.2_ug = Chla_6.2*vol_depths$Vol_m3[5]*1000) %>% 
  mutate(Chla_8_ug = Chla_8*vol_depths$Vol_m3[6]*1000) %>% 
  mutate(Chla_9_ug = Chla_9*vol_depths$Vol_m3[7]*1000)

# Add in thermocline depth information
flora_depths_3 <- left_join(flora_depths_2,la_results,by="DateTime") %>% 
  select(-St,-thermD,-N2,-SN2,-metaT,-metaB,-SmetaB,-SmetaT) %>% 
  mutate(SthermD = na.fill(na.approx(SthermD, na.rm=FALSE),"extend"))

### Thinking about how to designate Epi vs. Hypo what parameters depend on this:
# Mass of Epi vs. Mass of Hypo
# 'Outflow' to Epi from Hypo: concentration * outflow * scaled outflow = mass/day

# First, determine outflow concentration
flora_depths_3 <- flora_depths_3 %>% 
  mutate(SthermD = round(SthermD,digits=1)) %>% 
  mutate(epi_bottomg_depth_m = ifelse(SthermD > 9.0, 9.0,
                                      ifelse(SthermD > 7.0, 8.0,
                                             ifelse(SthermD > 6.0, 6.2,
                                                    ifelse(SthermD > 4.4, 5.0,
                                                           ifelse(SthermD > 3.0, 3.8,
                                                                  ifelse(SthermD > 1.6, 1.6,
                                                                         ifelse(SthermD > 0.0, 0.1, NA)))))))) %>% 
  mutate(hypo_top_depth_m = ifelse(SthermD <= 0.0, 0.1,
                                   ifelse(SthermD <= 1.6, 1.6,
                                          ifelse(SthermD <= 3.0, 3.8,
                                                 ifelse(SthermD <= 4.4, 5.0,
                                                        ifelse(SthermD <= 6.0, 6.2,
                                                               ifelse(SthermD <= 7.0, 8.0,
                                                                      ifelse(SthermD <= 9.0, 9.0, NA))))))))


### Calculate DOC/dt
# First need to figure out how the thermocline is changing
flora_depths_3 <- flora_depths_3 %>% 
  mutate(epi_ug = ifelse(epi_bottomg_depth_m==0.1,Chla_0.1_ug,
                        ifelse(epi_bottomg_depth_m==1.6,Chla_0.1_ug+Chla_1.6_ug,
                               ifelse(epi_bottomg_depth_m==3.8,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug,
                                      ifelse(epi_bottomg_depth_m==5.0,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug,
                                             ifelse(epi_bottomg_depth_m==6.2,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug,
                                                    ifelse(epi_bottomg_depth_m==8.0,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug+Chla_8_ug,
                                                           ifelse(epi_bottomg_depth_m==9.0,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug+Chla_8_ug+Chla_9_ug,NA)))))))) %>% 
  mutate(hypo_ug = ifelse(hypo_top_depth_m==0.1,Chla_0.1_ug+Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug+Chla_8_ug+Chla_9_ug,
                         ifelse(hypo_top_depth_m==1.6,Chla_1.6_ug+Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug+Chla_8_ug+Chla_9_ug,
                                ifelse(hypo_top_depth_m==3.8,Chla_3.8_ug+Chla_5_ug+Chla_6.2_ug+Chla_8_ug+Chla_9_ug,
                                       ifelse(hypo_top_depth_m==5.0,Chla_5_ug+Chla_6.2_ug+Chla_8_ug+Chla_9_ug,
                                              ifelse(hypo_top_depth_m==6.2,Chla_6.2_ug+Chla_8_ug+Chla_9_ug,
                                                     ifelse(hypo_top_depth_m==8.0,Chla_8_ug+Chla_9_ug,
                                                            ifelse(hypo_top_depth_m==9.0,Chla_9_ug,NA)))))))) %>% 
  mutate(epi_vol_m3 = ifelse(epi_bottomg_depth_m==0.1,sum(vol_depths$Vol_m3[1]),
                             ifelse(epi_bottomg_depth_m==1.6,sum(vol_depths$Vol_m3[1:2]),
                                    ifelse(epi_bottomg_depth_m==3.8,sum(vol_depths$Vol_m3[1:3]),
                                           ifelse(epi_bottomg_depth_m==5.0,sum(vol_depths$Vol_m3[1:4]),
                                                  ifelse(epi_bottomg_depth_m==6.2,sum(vol_depths$Vol_m3[1:5]),
                                                         ifelse(epi_bottomg_depth_m==8.0,sum(vol_depths$Vol_m3[1:6]),
                                                                ifelse(epi_bottomg_depth_m==9.0,sum(vol_depths$Vol_m3[1:7]),NA)))))))) %>% 
  mutate(hypo_vol_m3 = ifelse(hypo_top_depth_m==0.1,sum(vol_depths$Vol_m3[1:7]),
                              ifelse(hypo_top_depth_m==1.6,sum(vol_depths$Vol_m3[2:7]),
                                     ifelse(hypo_top_depth_m==3.8,sum(vol_depths$Vol_m3[3:7]),
                                            ifelse(hypo_top_depth_m==5.0,sum(vol_depths$Vol_m3[4:7]),
                                                   ifelse(hypo_top_depth_m==6.2,sum(vol_depths$Vol_m3[5:7]),
                                                          ifelse(hypo_top_depth_m==8.0,sum(vol_depths$Vol_m3[6:7]),
                                                                 ifelse(hypo_top_depth_m==9.0,sum(vol_depths$Vol_m3[7]),NA)))))))) %>% 
  mutate(epi_ugL = epi_ug/epi_vol_m3/1000,
         hypo_ugL = hypo_ug/hypo_vol_m3/1000)

### Load in CTD and YSI data ----
#need to import CTD observations from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/11/d771f5e9956304424c3bc0a39298a5ce" 
#infile1 <- paste0(getwd(),"/CTD_final_2013_2020.csv")
#download.file(inUrl1,infile1,method="curl")

ctd <- read.csv('CTD_final_2013_2020.csv') %>% #read in observed CTD data, which has multiple casts on the same day (problematic for comparison)
  filter(Reservoir=="FCR") %>%
  mutate(Date = as.POSIXct(strptime(Date, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:PAR_umolm2s)

ctd_50 <- ctd %>% 
  filter(Site==50) %>% 
  rename(time = Date)

# Import YSI observations from EDI
#inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/198/8/07ba1430528e01041435afc4c65fbeb6"
#infile1 <- paste0(getwd(),"/YSI_PAR_profiles_2013-2020.csv")
#download.file(inUrl1,infile1,method="curl")

ysi <- read_csv('YSI_PAR_profiles_2013-2020.csv') %>% 
  filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d",tz='EST'))) %>% 
  select(Reservoir:pH)

ysi_50 <- ysi %>% 
  filter(Site==50) %>% 
  rename(time = DateTime)

# Combine CTD and YSI data for site 50
# Select unique dates from both CTD and YSI casts
ysi_date_list <- as.data.frame(unique(as.Date(ysi_50$time)))
names(ysi_date_list)[1] <- "time"
ysi_date_list$ysi_fcr <- rep(-99,length(ysi_date_list$time))

ctd_date_list <- as.data.frame(unique(as.Date(ctd_50$time)))
names(ctd_date_list)[1] <- "time"
ctd_date_list$ctd_fcr <- rep(-99,length(ctd_date_list$time))

# Combine Unique dates list by date
fcr_dates <- merge(ysi_date_list, ctd_date_list, by="time", all.x=TRUE, all.y=TRUE)

### Merge data CTD and YSI datasets for FCR
fcr_merge <- merge(ctd_50, ysi_50, by="time", all.x=TRUE, all.y=TRUE)

# Find where there are Na values in the CTD data: need to do it for each column
ctd_fcr_na <- is.na(fcr_merge$Depth_m.x)
fcr_merge$Depth_m.x[ctd_fcr_na] <- fcr_merge$Depth_m.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$Temp_C.x)
fcr_merge$Temp_C.x[ctd_fcr_na] <- fcr_merge$Temp_C.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$DO_mgL.x)
fcr_merge$DO_mgL.x[ctd_fcr_na] <- fcr_merge$DO_mgL.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$DO_pSat)
fcr_merge$DO_pSat[ctd_fcr_na] <- fcr_merge$DOSat[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$Cond_uScm.x)
fcr_merge$Cond_uScm.x[ctd_fcr_na] <- fcr_merge$Cond_uScm.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$PAR_umolm2s.x)
fcr_merge$PAR_umolm2s.x[ctd_fcr_na] <- fcr_merge$PAR_umolm2s.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$ORP_mV.x)
fcr_merge$ORP_mV.x[ctd_fcr_na] <- fcr_merge$ORP_mV.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$pH.x)
fcr_merge$pH.x[ctd_fcr_na] <- fcr_merge$pH.y[ctd_fcr_na]

fcr_all <- fcr_merge %>% 
  select(time,Depth_m.x,Temp_C.x,DO_mgL.x,DO_pSat,Cond_uScm.x,Chla_ugL,Turb_NTU,pH.x,ORP_mV.x,PAR_umolm2s.x) %>% 
  rename(depth=Depth_m.x,Temp_C=Temp_C.x,DO_mgL=DO_mgL.x,Cond_uScm=Cond_uScm.x,pH=pH.x,ORP_mV=ORP_mV.x,PAR_umolm2s=PAR_umolm2s.x)

fcr_date_list <- as.data.frame(unique(as.Date(fcr_all$time)))

## Average across date and depth
fcr_all <- fcr_all %>% group_by(time,depth) %>% summarize_all(funs(mean),na.rm=TRUE)

### Create a dataframe for CTD/YSI parameters at each sampling depth
depths<- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

#Initialize an empty matrix with the correct number of rows and columns 
temp<-matrix(data=NA, ncol=ncol(fcr_all), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final<-matrix(data=NA, ncol=1, nrow=0)
dates<-unique(fcr_all$time)

#create a function to chose the matching depth closest to our focal depths
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j=dates[i]
  q <- subset(fcr_all, fcr_all$time == j)
  
  layer1 <- q[q[, "depth"] == closest(q$depth,0.1),][1,]
  layer2<- q[q[, "depth"] == closest(q$depth,1.6),][1,]
  layer3<- q[q[, "depth"] == closest(q$depth,3.8),][1,]
  layer4<- q[q[, "depth"] == closest(q$depth,5.0),][1,]
  layer5<- q[q[, "depth"] == closest(q$depth,6.2),][1,]
  layer6<- q[q[, "depth"] == closest(q$depth,8.0),][1,]
  layer7<- q[q[, "depth"] == closest(q$depth,9.0),][1,]
  
  temp<-rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(fcr_all))+1)] <- depths
  colnames(temp)[((ncol(fcr_all))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
ctd_50_depths <- as.data.frame(super_final) %>%
  select(-c(1,depth)) %>%
  rename(depth = new_depth) %>%
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST")))

ctd_50_depths$Temp_C <- as.numeric(ctd_50_depths$Temp_C)
ctd_50_depths$depth <- as.numeric(ctd_50_depths$depth)
ctd_50_depths$Cond_uScm <- as.numeric(ctd_50_depths$Cond_uScm)
ctd_50_depths$DO_mgL <- as.numeric(ctd_50_depths$DO_mgL)

# Calculate %DO Saturation using rMR
ctd_do_calc <- ctd_50_depths %>% 
  select(time,Temp_C,DO_mgL,Cond_uScm,depth)

ctd_do_calc <- na.omit(ctd_do_calc)

ctd_do_calc <- ctd_do_calc %>% 
  mutate(calc_DO_pSat = DO.saturation(DO.mgl = ctd_do_calc$DO_mgL,temp.C = ctd_do_calc$Temp_C,
                                      elevation.m = 509, bar.press = NULL, bar.units = NULL,
                                      ctd_do_calc$Cond_uScm,salinity.units = "uS")*100)

ctd_50_depths_2 <- left_join(ctd_50_depths,ctd_do_calc,by=c("time","depth","Temp_C","DO_mgL","Cond_uScm"))

ctd_50_depths_2 <- ctd_50_depths_2 %>% 
  mutate(calc_DO_pSat = ifelse(is.na(calc_DO_pSat),DO_pSat,calc_DO_pSat))

ctd_50_depths_2$calc_DO_pSat <- as.numeric(ctd_50_depths_2$calc_DO_pSat)

# Filter out 2019-05-30: Funky cast!
ctd_50_depths_2 <- ctd_50_depths_2 %>% 
  filter(time != as.POSIXct("2019-05-30"))

# Calculate VW Epi and VW Hypo
ctd_50_depths_3 <- ctd_50_depths_2 %>% 
  filter(time >= as.POSIXct("2015-01-01")) %>% 
  select(time,Temp_C,DO_mgL,calc_DO_pSat,depth) %>% 
  pivot_wider(names_from = depth, values_from = c(Temp_C,DO_mgL,calc_DO_pSat), values_fil = NA, values_fn = mean)

ctd_50_depths_3 <- na.omit(ctd_50_depths_3)

ctd_50_depths_3 <- ctd_50_depths_3 %>% 
  mutate(Temp_C_0.1 = Temp_C_0.1*vol_depths$Vol_m3[1]*1000) %>% 
  mutate(Temp_C_1.6 = Temp_C_1.6*vol_depths$Vol_m3[2]*1000) %>% 
  mutate(Temp_C_3.8 = Temp_C_3.8*vol_depths$Vol_m3[3]*1000) %>% 
  mutate(Temp_C_5 = Temp_C_5*vol_depths$Vol_m3[4]*1000) %>% 
  mutate(Temp_C_6.2 = Temp_C_6.2*vol_depths$Vol_m3[5]*1000) %>% 
  mutate(Temp_C_8 = Temp_C_8*vol_depths$Vol_m3[6]*1000) %>% 
  mutate(Temp_C_9 = Temp_C_9*vol_depths$Vol_m3[7]*1000) %>% 
  mutate(DO_0.1_mg = DO_mgL_0.1*vol_depths$Vol_m3[1]*1000) %>% 
  mutate(DO_1.6_mg = DO_mgL_1.6*vol_depths$Vol_m3[2]*1000) %>% 
  mutate(DO_3.8_mg = DO_mgL_3.8*vol_depths$Vol_m3[3]*1000) %>% 
  mutate(DO_5_mg = DO_mgL_5*vol_depths$Vol_m3[4]*1000) %>% 
  mutate(DO_6.2_mg = DO_mgL_6.2*vol_depths$Vol_m3[5]*1000) %>% 
  mutate(DO_8_mg = DO_mgL_8*vol_depths$Vol_m3[6]*1000) %>% 
  mutate(DO_9_mg = DO_mgL_9*vol_depths$Vol_m3[7]*1000) %>% 
  mutate(calc_DO_pSat_0.1 = calc_DO_pSat_0.1*vol_depths$Vol_m3[1]*1000) %>% 
  mutate(calc_DO_pSat_1.6 = calc_DO_pSat_1.6*vol_depths$Vol_m3[2]*1000) %>% 
  mutate(calc_DO_pSat_3.8 = calc_DO_pSat_3.8*vol_depths$Vol_m3[3]*1000) %>% 
  mutate(calc_DO_pSat_5 = calc_DO_pSat_5*vol_depths$Vol_m3[4]*1000) %>% 
  mutate(calc_DO_pSat_6.2 = calc_DO_pSat_6.2*vol_depths$Vol_m3[5]*1000) %>% 
  mutate(calc_DO_pSat_8 = calc_DO_pSat_8*vol_depths$Vol_m3[6]*1000) %>% 
  mutate(calc_DO_pSat_9 = calc_DO_pSat_9*vol_depths$Vol_m3[7]*1000) %>% 
  rename(DateTime = time)

# Add in thermocline depth information
ctd_depths <- left_join(ctd_50_depths_3,la_results,by="DateTime") %>% 
  select(-St,-thermD,-N2,-SN2,-metaT,-metaB,-SmetaB,-SmetaT) %>% 
  mutate(SthermD = na.fill(na.approx(SthermD, na.rm=FALSE),"extend"))

### Thinking about how to designate Epi vs. Hypo what parameters depend on this:
# Mass of Epi vs. Mass of Hypo
# 'Outflow' to Epi from Hypo: concentration * outflow * scaled outflow = mass/day

# First, determine outflow concentration
ctd_depths <- ctd_depths %>% 
  mutate(SthermD = round(SthermD,digits=1)) %>% 
  mutate(epi_bottomg_depth_m = ifelse(SthermD > 9.0, 9.0,
                                      ifelse(SthermD > 7.0, 8.0,
                                             ifelse(SthermD > 6.0, 6.2,
                                                    ifelse(SthermD > 4.4, 5.0,
                                                           ifelse(SthermD > 3.0, 3.8,
                                                                  ifelse(SthermD > 1.6, 1.6,
                                                                         ifelse(SthermD > 0.0, 0.1, NA)))))))) %>% 
  mutate(hypo_top_depth_m = ifelse(SthermD <= 0.0, 0.1,
                                   ifelse(SthermD <= 1.6, 1.6,
                                          ifelse(SthermD <= 3.0, 3.8,
                                                 ifelse(SthermD <= 4.4, 5.0,
                                                        ifelse(SthermD <= 6.0, 6.2,
                                                               ifelse(SthermD <= 7.0, 8.0,
                                                                      ifelse(SthermD <= 9.0, 9.0, NA))))))))


### Calculate DOC/dt
# First need to figure out how the thermocline is changing
ctd_depths <- ctd_depths %>% 
  mutate(epi_Temp_C = ifelse(epi_bottomg_depth_m==0.1,Temp_C_0.1,
                         ifelse(epi_bottomg_depth_m==1.6,Temp_C_0.1+Temp_C_1.6,
                                ifelse(epi_bottomg_depth_m==3.8,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8,
                                       ifelse(epi_bottomg_depth_m==5.0,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8+Temp_C_5,
                                              ifelse(epi_bottomg_depth_m==6.2,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8+Temp_C_5+Temp_C_6.2,
                                                     ifelse(epi_bottomg_depth_m==8.0,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8+Temp_C_5+Temp_C_6.2+Temp_C_8,
                                                            ifelse(epi_bottomg_depth_m==9.0,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8+Temp_C_5+Temp_C_6.2+Temp_C_8+Temp_C_9,NA)))))))) %>% 
  mutate(hypo_Temp_C = ifelse(hypo_top_depth_m==0.1,Temp_C_0.1+Temp_C_1.6+Temp_C_3.8+Temp_C_5+Temp_C_6.2+Temp_C_8+Temp_C_9,
                          ifelse(hypo_top_depth_m==1.6,Temp_C_1.6+Temp_C_3.8+Temp_C_5+Temp_C_6.2+Temp_C_8+Temp_C_9,
                                 ifelse(hypo_top_depth_m==3.8,Temp_C_3.8+Temp_C_5+Temp_C_6.2+Temp_C_8+Temp_C_9,
                                        ifelse(hypo_top_depth_m==5.0,Temp_C_5+Temp_C_6.2+Temp_C_8+Temp_C_9,
                                               ifelse(hypo_top_depth_m==6.2,Temp_C_6.2+Temp_C_8+Temp_C_9,
                                                      ifelse(hypo_top_depth_m==8.0,Temp_C_8+Temp_C_9,
                                                             ifelse(hypo_top_depth_m==9.0,Temp_C_9,NA)))))))) %>% 
  mutate(epi_DO_mg = ifelse(epi_bottomg_depth_m==0.1,DO_0.1_mg,
                             ifelse(epi_bottomg_depth_m==1.6,DO_0.1_mg+DO_1.6_mg,
                                    ifelse(epi_bottomg_depth_m==3.8,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg,
                                           ifelse(epi_bottomg_depth_m==5.0,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg+DO_5_mg,
                                                  ifelse(epi_bottomg_depth_m==6.2,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg+DO_5_mg+DO_6.2_mg,
                                                         ifelse(epi_bottomg_depth_m==8.0,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg+DO_5_mg+DO_6.2_mg+DO_8_mg,
                                                                ifelse(epi_bottomg_depth_m==9.0,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg+DO_5_mg+DO_6.2_mg+DO_8_mg+DO_9_mg,NA)))))))) %>% 
  mutate(hypo_DO_mg = ifelse(hypo_top_depth_m==0.1,DO_0.1_mg+DO_1.6_mg+DO_3.8_mg+DO_5_mg+DO_6.2_mg+DO_8_mg+DO_9_mg,
                              ifelse(hypo_top_depth_m==1.6,DO_1.6_mg+DO_3.8_mg+DO_5_mg+DO_6.2_mg+DO_8_mg+DO_9_mg,
                                     ifelse(hypo_top_depth_m==3.8,DO_3.8_mg+DO_5_mg+DO_6.2_mg+DO_8_mg+DO_9_mg,
                                            ifelse(hypo_top_depth_m==5.0,DO_5_mg+DO_6.2_mg+DO_8_mg+DO_9_mg,
                                                   ifelse(hypo_top_depth_m==6.2,DO_6.2_mg+DO_8_mg+DO_9_mg,
                                                          ifelse(hypo_top_depth_m==8.0,DO_8_mg+DO_9_mg,
                                                                 ifelse(hypo_top_depth_m==9.0,DO_9_mg,NA)))))))) %>% 
  mutate(epi_calc_DO_pSat = ifelse(epi_bottomg_depth_m==0.1,calc_DO_pSat_0.1,
                             ifelse(epi_bottomg_depth_m==1.6,calc_DO_pSat_0.1+calc_DO_pSat_1.6,
                                    ifelse(epi_bottomg_depth_m==3.8,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8,
                                           ifelse(epi_bottomg_depth_m==5.0,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5,
                                                  ifelse(epi_bottomg_depth_m==6.2,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2,
                                                         ifelse(epi_bottomg_depth_m==8.0,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8,
                                                                ifelse(epi_bottomg_depth_m==9.0,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,NA)))))))) %>% 
  mutate(hypo_calc_DO_pSat = ifelse(hypo_top_depth_m==0.1,calc_DO_pSat_0.1+calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,
                              ifelse(hypo_top_depth_m==1.6,calc_DO_pSat_1.6+calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,
                                     ifelse(hypo_top_depth_m==3.8,calc_DO_pSat_3.8+calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,
                                            ifelse(hypo_top_depth_m==5.0,calc_DO_pSat_5+calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,
                                                   ifelse(hypo_top_depth_m==6.2,calc_DO_pSat_6.2+calc_DO_pSat_8+calc_DO_pSat_9,
                                                          ifelse(hypo_top_depth_m==8.0,calc_DO_pSat_8+calc_DO_pSat_9,
                                                                 ifelse(hypo_top_depth_m==9.0,calc_DO_pSat_9,NA)))))))) %>% 
  mutate(epi_vol_m3 = ifelse(epi_bottomg_depth_m==0.1,sum(vol_depths$Vol_m3[1]),
                             ifelse(epi_bottomg_depth_m==1.6,sum(vol_depths$Vol_m3[1:2]),
                                    ifelse(epi_bottomg_depth_m==3.8,sum(vol_depths$Vol_m3[1:3]),
                                           ifelse(epi_bottomg_depth_m==5.0,sum(vol_depths$Vol_m3[1:4]),
                                                  ifelse(epi_bottomg_depth_m==6.2,sum(vol_depths$Vol_m3[1:5]),
                                                         ifelse(epi_bottomg_depth_m==8.0,sum(vol_depths$Vol_m3[1:6]),
                                                                ifelse(epi_bottomg_depth_m==9.0,sum(vol_depths$Vol_m3[1:7]),NA)))))))) %>% 
  mutate(hypo_vol_m3 = ifelse(hypo_top_depth_m==0.1,sum(vol_depths$Vol_m3[1:7]),
                              ifelse(hypo_top_depth_m==1.6,sum(vol_depths$Vol_m3[2:7]),
                                     ifelse(hypo_top_depth_m==3.8,sum(vol_depths$Vol_m3[3:7]),
                                            ifelse(hypo_top_depth_m==5.0,sum(vol_depths$Vol_m3[4:7]),
                                                   ifelse(hypo_top_depth_m==6.2,sum(vol_depths$Vol_m3[5:7]),
                                                          ifelse(hypo_top_depth_m==8.0,sum(vol_depths$Vol_m3[6:7]),
                                                                 ifelse(hypo_top_depth_m==9.0,sum(vol_depths$Vol_m3[7]),NA)))))))) %>% 
  mutate(epi_Temp_C_vw = epi_Temp_C/epi_vol_m3/1000,
         hypo_Temp_C_vw = hypo_Temp_C/hypo_vol_m3/1000,
         epi_DO_mgL = epi_DO_mg/epi_vol_m3/1000,
         hypo_DO_mgL = hypo_DO_mg/hypo_vol_m3/1000,
         epi_DO_pSat_vw = epi_calc_DO_pSat/epi_vol_m3/1000,
         hypo_DO_pSat_vw = hypo_calc_DO_pSat/hypo_vol_m3/1000)

# Remove 2019-05-30: weird cast!
ctd_depths <- ctd_depths %>% 
  filter(DateTime != as.POSIXct("2019-05-30 EST"))

# Plot!
temp <- ggplot(ctd_depths)+
  geom_line(mapping=aes(x=DateTime,y=epi_Temp_C_vw,color="Epi"))+
  geom_point(mapping=aes(x=DateTime,y=epi_Temp_C_vw,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_Temp_C_vw,color="Hypo"))+
  geom_point(mapping=aes(x=DateTime,y=hypo_Temp_C_vw,color="Hypo"))

do_mgL <- ggplot(ctd_depths)+
  geom_line(mapping=aes(x=DateTime,y=epi_DO_mgL,color="Epi"))+
  geom_point(mapping=aes(x=DateTime,y=epi_DO_mgL,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO_mgL,color="Hypo"))+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO_mgL,color="Hypo"))

do_sat <- ggplot(ctd_depths)+
  geom_line(mapping=aes(x=DateTime,y=epi_DO_pSat_vw,color="Epi"))+
  geom_point(mapping=aes(x=DateTime,y=epi_DO_pSat_vw,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO_pSat_vw,color="Hypo"))+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO_pSat_vw,color="Hypo"))

# Plus Flora plots
flora <- ggplot(flora_depths_3)+
  geom_line(mapping=aes(x=DateTime,y=epi_ugL,color="Epi"))+
  geom_point(mapping=aes(x=DateTime,y=epi_ugL,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_ugL,color="Hypo"))+
  geom_point(mapping=aes(x=DateTime,y=hypo_ugL,color="Hypo"))

ggarrange(temp,do_mgL,do_sat,flora,nrow=4,ncol=1)

ggsave("./Fig_Output/Env_Params.jpg",width=10,height=12,units="in",dpi=320)
