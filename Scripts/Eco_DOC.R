### Updated script to include 'whole-ecosystem' DOC processing (epi and hypo)
### Following comments from PCH and CCC
### Script is drawing heavily from DOC_DO.R
### 11 Aug 2021, A Hounshell

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,rMR,lme4,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics)

### Define boxes (epi vs. hypo) ----
# Important Lake Analyzer thermocline results to determine median thermocline depth
la_results <- read.csv("./Data/20210603_LA_FCR_results.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST")))

la_results_strat <- la_results %>% 
  mutate(month = month(DateTime)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(month %in% c(6,7,8,9))

thermo <- la_results_strat %>% 
  summarise_all(median,na.rm=TRUE)

thermo_year <- la_results_strat %>% 
  group_by(year) %>% 
  summarise_all(median,na.rm=TRUE)

thermo_month <- la_results_strat %>% 
  group_by(month) %>% 
  summarise_all(median,na.rm=TRUE)

### Load DOC data ----
# Updated: 19 Apr 2021 with 2020 data!
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98" 
#infile1 <- paste0(getwd(),"/Data/chem.csv")
#download.file(inUrl1,infile1,method="curl")

chem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

chem_100 <- chem %>% 
  filter(Site == 100)

chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DOC_mgL <= 15)

# Plot to see what type of data we have

chem_50 %>% 
  filter(year == 2020) %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)))+
  geom_point()+
  geom_line()

# Separate and calculate VW [DOC] for Epi and Hypo: from 2016 to 2020
doc_epi <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  drop_na() %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8)) %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_") %>% 
  mutate(DOC_0.1 = na.fill(na.approx(DOC_0.1,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_1.6 = na.fill(na.approx(DOC_1.6,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_3.8 = na.fill(na.approx(DOC_3.8,na.rm=FALSE),"extend")) %>% 
  mutate(VW_Epi_DOC_mgL = ((DOC_0.1*1.38*10^5)+(DOC_1.6*8.91*10^4)+(DOC_3.8*5.96*10^4))/((1.38*10^5)+(8.91*10^4)+(5.96*10^4))) %>% 
  mutate(month = month(DateTime))

doc_hypo <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  drop_na() %>% 
  filter(Depth_m %in% c(5.0,6.2,8.0,9.0)) %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_") %>% 
  mutate(DOC_5 = na.fill(na.approx(DOC_5,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_6.2 = na.fill(na.approx(DOC_6.2,na.rm=FALSE),"extend")) %>%   
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend")) %>% 
  mutate(VW_Hypo_DOC_mgL = ((DOC_5*4.02*10^4)+(DOC_6.2*1.40*10^4)+(DOC_8*1.40*10^4)+(DOC_9*1.95*10^3))/((4.02*10^4)+(1.40*10^4)+(1.40*10^4)+(1.95*10^3))) %>% 
  mutate(month = month(DateTime))

# Plot VW Epi and VW Hypo: WILL NEED TO UPDATE FOR MS
ggplot()+
  geom_point(doc_epi,mapping=aes(x=DateTime,y=VW_Epi_DOC_mgL,color="VW Epi"))+
  geom_line(doc_epi,mapping=aes(x=DateTime,y=VW_Epi_DOC_mgL,color="VW Epi"))+
  geom_point(doc_hypo,mapping=aes(x=DateTime,y=VW_Hypo_DOC_mgL,color="VW Hypo"))+
  geom_line(doc_hypo,mapping=aes(x=DateTime,y=VW_Hypo_DOC_mgL,color="VW Hypo"))+
  theme_classic(base_size=15)

### Load in Inflow data ----

