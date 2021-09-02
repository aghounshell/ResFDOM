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

### Load in Inflow data ----
# Weir discharge/temperature
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/7/f5fa5de4b49bae8373f6e7c1773b026e" 
#infile1 <- paste0(getwd(),"/Data/inflow_for_EDI_2013_10Jan2021.csv")
#download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/inflow_for_EDI_2013_10Jan2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:VT_Temp_C)

# Average inflow by day
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_all(funs(mean(.,na.rm=TRUE))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

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

doc_box_full <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2020-12-31",tz="EST"),by="days"))
doc_box_full <- doc_box_full %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2020-12-31", tz = "EST"), by = "days")`)
doc_box_full <- left_join(doc_box_full,doc_box,by="DateTime")

doc_box_full <- doc_box_full %>% 
  mutate(DOC_0.1 = na.fill(na.approx(DOC_0.1,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_1.6 = na.fill(na.approx(DOC_1.6,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_3.8 = na.fill(na.approx(DOC_3.8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_5 = na.fill(na.approx(DOC_5,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_6.2 = na.fill(na.approx(DOC_6.2,na.rm=FALSE),"extend")) %>%   
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend"))

doc_box_full <- doc_box_full %>% 
  mutate(DOC_0.1_g = DOC_0.1*vol_depths$Vol_m3[1],
         DOC_1.6_g = DOC_1.6*vol_depths$Vol_m3[2],
         DOC_3.8_g = DOC_3.8*vol_depths$Vol_m3[3],
         DOC_5_g = DOC_5*vol_depths$Vol_m3[4],
         DOC_6.2_g = DOC_6.2*vol_depths$Vol_m3[5],
         DOC_8_g = DOC_8*vol_depths$Vol_m3[6],
         DOC_9_g = DOC_9*vol_depths$Vol_m3[7])

# Second, extrapolate inflow concentrations and calculate daily mass
doc_inflow_full <- as.data.frame(seq(as.POSIXct("2015-01-01",tz="EST"),as.POSIXct("2020-12-31",tz="EST"),by="days"))
doc_inflow_full <- doc_inflow_full %>% 
  rename(DateTime = `seq(as.POSIXct("2015-01-01", tz = "EST"), as.POSIXct("2020-12-31", tz = "EST"), by = "days")`)

doc_inflow_full <- left_join(doc_inflow_full,inflow_daily,by="DateTime") %>% 
  select(DateTime,WVWA_Flow_cms,VT_Flow_cms) 

doc_inflow_full <- doc_inflow_full %>%  
  mutate(full_flow_cms = ifelse(WVWA_Flow_cms == "NaN", VT_Flow_cms, WVWA_Flow_cms))

doc_inflow_full <- left_join(doc_inflow_full,chem_100,by="DateTime") %>% 
  select(DateTime,full_flow_cms,DOC_mgL)

doc_inflow_full <- doc_inflow_full %>% 
  mutate(DOC_mgL = na.fill(na.approx(DOC_mgL,na.rm=FALSE),"extend"),
         full_flow_cms = na.fill(na.approx(full_flow_cms,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_100_g = DOC_mgL*full_flow_cms*60*60*24)

# DOC Box data
doc_box_data <- left_join(doc_box_full, doc_inflow_full,by="DateTime")

# Add in thermocline depth information
doc_box_data <- left_join(doc_box_data,la_results,by="DateTime") %>% 
  select(-St,-thermD,-N2,-SN2)

doc_box_data <- doc_box_data %>% 
  mutate(SthermD = na.fill(na.approx(SthermD, na.rm=FALSE),"extend"))

### Thinking about how to designate Epi vs. Hypo what parameters depend on this:
# Mass of Epi vs. Mass of Hypo
# 'Outflow' to Epi from Hypo: concentration * outflow * scaled outflow = mass/day

# First, determine outflow concentration
doc_box_data <- doc_box_data %>% 
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

# Calculate Epi and Hypo outflow mass
doc_box_data <- doc_box_data %>% 
  mutate(Epi_outflow_g = DOC_0.1*full_flow_cms*60*60*24,
         Hypo_outflow_g = ifelse(hypo_top_depth_m==0.1,DOC_0.1*full_flow_cms*60*60*24,
                                 ifelse(hypo_top_depth_m==1.6,DOC_1.6*full_flow_cms*60*60*24,
                                        ifelse(hypo_top_depth_m==3.8,DOC_3.8*full_flow_cms*60*60*24,
                                               ifelse(hypo_top_depth_m==5.0,DOC_5*full_flow_cms*60*60*24,
                                                      ifelse(hypo_top_depth_m==6.2,DOC_6.2*full_flow_cms*60*60*24,
                                                             ifelse(hypo_top_depth_m==8.0,DOC_8*full_flow_cms*60*60*24,
                                                                    ifelse(hypo_top_depth_m==9.0,DOC_9*full_flow_cms*60*60*24,0))))))))

### Calculate DOC/dt
# First need to figure out how the thermocline is changing
doc_box_data <- doc_box_data %>% 
  mutate(epi_g = ifelse(epi_bottomg_depth_m==0.1,DOC_0.1_g,
                        ifelse(epi_bottomg_depth_m==1.6,DOC_0.1_g+DOC_1.6_g,
                               ifelse(epi_bottomg_depth_m==3.8,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g,
                                      ifelse(epi_bottomg_depth_m==5.0,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g+DOC_5_g,
                                             ifelse(epi_bottomg_depth_m==6.2,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g+DOC_5_g+DOC_6.2_g,
                                                    ifelse(epi_bottomg_depth_m==8.0,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g+DOC_5_g+DOC_6.2_g+DOC_8_g,
                                                           ifelse(epi_bottomg_depth_m==9.0,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g+DOC_5_g+DOC_6.2_g+DOC_8_g+DOC_9_g,NA)))))))) %>% 
  mutate(hypo_g = ifelse(hypo_top_depth_m==0.1,DOC_0.1_g+DOC_1.6_g+DOC_3.8_g+DOC_5_g+DOC_6.2_g+DOC_8_g+DOC_9_g,
                         ifelse(hypo_top_depth_m==1.6,DOC_1.6_g+DOC_3.8_g+DOC_5_g+DOC_6.2_g+DOC_8_g+DOC_9_g,
                                ifelse(hypo_top_depth_m==3.8,DOC_3.8_g+DOC_5_g+DOC_6.2_g+DOC_8_g+DOC_9_g,
                                       ifelse(hypo_top_depth_m==5.0,DOC_5_g+DOC_6.2_g+DOC_8_g+DOC_9_g,
                                              ifelse(hypo_top_depth_m==6.2,DOC_6.2_g+DOC_8_g+DOC_9_g,
                                                     ifelse(hypo_top_depth_m==8.0,DOC_8_g+DOC_9_g,
                                                            ifelse(hypo_top_depth_m==9.0,DOC_9_g,NA)))))))) %>% 
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
  mutate(d_epi_g_dt = NA,
         d_hypo_g_dt = NA,
         Entr = NA,
         DOC_entr_g= NA)

# Calculate DOC/dt for epi and hypo
for(i in 2:length(doc_box_data$DateTime)){
  doc_box_data$d_epi_g_dt[i] = (doc_box_data$epi_g[i]-doc_box_data$epi_g[i-1])
}

for(i in 2:length(doc_box_data$DateTime)){
  doc_box_data$d_hypo_g_dt[i] = (doc_box_data$hypo_g[i]-doc_box_data$hypo_g[i-1])
}

### Calculate change due to entrainment (movement of thermocline)
# Loosely following FCR_DOCModel_edited_19May17 from Carey et al. 2018

# Entrainment
#if Entr is positive, then epi is getting bigger and hypo is getting smaller; 
#if Entr is negative, then hypo is getting bigger and epi is getting smaller

# First need to create data frame of mass by date and volume
doc_entr <- doc_box_data %>% 
  select(DateTime,DOC_0.1_g:DOC_9_g) %>% 
  rename(DOC_0.1 = DOC_0.1_g, DOC_1.6 = DOC_1.6_g, DOC_3.8 = DOC_3.8_g, DOC_5 = DOC_5_g, DOC_6.2 = DOC_6.2_g, DOC_8 = DOC_8_g, DOC_9 = DOC_9_g) %>% 
  rowwise() %>% 
  mutate(DOC_0.138 = sum(DOC_1.6,DOC_3.8)) %>% #Create a new 'depth' where bottom of epi goes from 0.1 to 3.8 m (0.138)
  mutate(DOC_0.15 = sum(DOC_1.6,DOC_3.8,DOC_5)) %>% # Create a new 'depth' where bottom of epi goes from 0.1 to 5 m (0.15)
  mutate(DOC_0.162 = sum(DOC_1.6,DOC_3.8,DOC_5,DOC_6.2)) %>% # Creating a new 'depth' where bottom of epi goes from 0.1 to 6.2 m (0.162)
  mutate(DOC_0.18 = sum(DOC_1.6,DOC_3.8,DOC_5,DOC_6.2,DOC_8)) %>% # Create a new depth where bottom of epi goes from 0.1 to 8 m (0.18)
  mutate(DOC_1.65 = sum(DOC_3.8,DOC_5)) %>% # Creating a new 'depth' where bottom of epi goes from 1.6 to 5 m (1.65)
  mutate(DOC_1.662 = sum(DOC_3.8,DOC_5,DOC_6.2)) %>%# Create a new 'depth' where bottom of epi goes from 1.6 to 6.2 m (1.662)
  mutate(DOC_1.68 = sum(DOC_3.8,DOC_5,DOC_6.2,DOC_8)) %>% # Create a new depth where bottom of epi goes from 1.6 to 8 m (1.68)
  mutate(DOC_3.862 = sum(DOC_5,DOC_6.2)) %>%# Creating a new 'depth' where bottom of epi goes from 3.8 to 6.2 m (3.862)
  mutate(DOC_5.8 = sum(DOC_6.2,DOC_8)) %>% # Create a new 'depth' where bottom of epi goes from 5 to 8 m (5.8)
  pivot_longer(!DateTime,names_to = "depth", values_to = "DOC_g",names_prefix ="DOC_") %>% 
  mutate(vol_m3 = ifelse(depth == 0.1, vol_depths$Vol_m3[1],
                         ifelse(depth == 0.138, sum(vol_depths$Vol_m3[2:3]),
                                ifelse(depth == 0.15, sum(vol_depths$Vol_m3[2:4]),
                                       ifelse(depth == 0.162, sum(vol_depths$Vol_m3[2:5]),
                                              ifelse(depth == 0.18, sum(vol_depths$Vol_m3[2:6]),
                                                     ifelse(depth == 1.6, vol_depths$Vol_m3[2],
                                                            ifelse(depth == 1.65, sum(vol_depths$Vol_m3[3:4]),
                                                                   ifelse(depth == 1.662, sum(vol_depths$Vol_m3[3:5]),
                                                                          ifelse(depth == 1.68, sum(vol_depths$Vol_m3[3:6]),
                                                                                 ifelse(depth == 3.8, vol_depths$Vol_m3[3],
                                                                                        ifelse(depth == 3.862, sum(vol_depths$Vol_m3[4:5]),
                                                                                               ifelse(depth == 5, vol_depths$Vol_m3[4],
                                                                                                      ifelse(depth == 5.8, sum(vol_depths$Vol_m3[5:6]),
                                                                                                             ifelse(depth == 6.2, vol_depths$Vol_m3[5],
                                                                                                                    ifelse(depth == 8, vol_depths$Vol_m3[6],
                                                                                                                           ifelse(depth == 9, vol_depths$Vol_m3[7], NA)))))))))))))))))

# Need to figure out: 1. When the thermocline is moving; 2. How far the thermocline moves; 3. Then calculate
# how much mass is moved b/c of the change in thermocline
doc_entr <- as.data.frame(doc_entr)
doc_box_data <- as.data.frame(doc_box_data)
for(i in 2:length(doc_box_data$DateTime)){
  # 1. Figure out if the thermocline has moved
  if(doc_box_data$epi_vol_m3[i] > doc_box_data$epi_vol_m3[i-1]){
    doc_box_data$Entr[i] = doc_box_data$epi_vol_m3[i] - doc_box_data$epi_vol_m3[i-1]
    x = which(doc_box_data$DateTime[i]==doc_entr$DateTime & abs(round(doc_box_data$Entr[i],2))==round(doc_entr$vol_m3,2))
    doc_box_data$DOC_entr_g[i] = doc_entr$DOC_g[x]
  }
  else if(doc_box_data$epi_vol_m3[i] < doc_box_data$epi_vol_m3[i-1]){
    doc_box_data$Entr[i] = doc_box_data$epi_vol_m3[i] - doc_box_data$epi_vol_m3[i-1]
    x = which(doc_box_data$DateTime[i]==doc_entr$DateTime & abs(round(doc_box_data$Entr[i],2))==round(doc_entr$vol_m3,2))
    doc_box_data$DOC_entr_g[i] = -doc_entr$DOC_g[x]
  }
  else{
    doc_box_data$Entr[i] = 0
    doc_box_data$DOC_entr_g[i] = 0
  }
}

### Then calculate 'processing' for epi and hypo ----
doc_box_data <- doc_box_data %>% 
  mutate(epi_processing_g = d_epi_g_dt-DOC_100_g*0.74-Hypo_outflow_g*0.26+Epi_outflow_g-DOC_entr_g,
         epi_processing_mgL = epi_processing_g/epi_vol_m3,
         epi_DOC_mgL = epi_g/epi_vol_m3,
         hypo_processing_g = d_hypo_g_dt-DOC_100_g*0.26+Hypo_outflow_g*0.26+DOC_entr_g,
         hypo_processing_mgL = hypo_processing_g/hypo_vol_m3,
         hypo_DOC_mgL = hypo_g/hypo_vol_m3,
         DOC_entr_mgL = DOC_entr_g/abs(Entr))

# Plot Epi and Hypo processing for days w/ data
doc_box_sel <- chem_50 %>% 
  distinct(DateTime)

doc_box_sel <- left_join(doc_box_sel,doc_box_data,by="DateTime")

# Plot Epi and Hypo DOC VW concentrations (mg/L)
doc_conc <- ggplot(chem_50,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)))+
  geom_line()+
  geom_point()

doc_box_conc <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=epi_DOC_mgL,color="Epi"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=epi_DOC_mgL,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_DOC_mgL,color="Hypo"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=hypo_DOC_mgL,color="Hypo"))

ggarrange(doc_conc,doc_box_conc,nrow=2,ncol=1)

ggsave("./Fig_Output/DOC_Concentrations.jpg",width=10,height=12,units="in",dpi=320)

# Plot all processes for the epi
epi_flows <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=DOC_100_g*0.74,color="Epi_Inflow"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=DOC_100_g*0.74,color="Epi_Inflow"))+
  geom_line(mapping=aes(x=DateTime,y=Epi_outflow_g,color="Epi_Outflow"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=Epi_outflow_g,color="Epi_Outflow"))+
  geom_line(mapping=aes(x=DateTime,y=Hypo_outflow_g*0.26,color="Hypo_Inflow"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=Hypo_outflow_g*0.26,color="Hypo_Inflow"))
  
epi_change <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=d_epi_g_dt,color="DOC/dt"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=d_epi_g_dt,color="DOC/dt"))+
  geom_line(mapping=aes(x=DateTime,y=DOC_entr_g,color="Entrain"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=DOC_entr_g,color="Entrain"))

epi_process <- ggplot(doc_box_data,mapping=aes(x=DateTime,y=epi_processing_mgL,color="Epi_Processing"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=epi_processing_mgL,color="Epi_Processing"))

ggarrange(epi_flows,epi_change,epi_process,nrow=3,ncol=1)

ggsave("./Fig_Output/Epi_DOC_Processes.jpg",width=10,height=12,units="in",dpi=320)

# Plot all processes for the hypo
hypo_flows <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=DOC_100_g*0.26,color="Hypo_Inflow"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=DOC_100_g*0.26,color="Hypo_Inflow"))+
  geom_line(mapping=aes(x=DateTime,y=Hypo_outflow_g*0.26,color="Hypo_Outflow"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=Hypo_outflow_g*0.26,color="Hypo_Outflow"))

hypo_change <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=d_hypo_g_dt,color="DOC/dt"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=d_hypo_g_dt,color="DOC/dt"))+
  geom_line(mapping=aes(x=DateTime,y=-1*DOC_entr_g,color="Entrain"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=-1*DOC_entr_g,color="Entrain"))

hypo_process <- ggplot(doc_box_data)+
  geom_line(mapping=aes(x=DateTime,y=hypo_processing_mgL,color="Hypo_Processing"))+
  geom_point(doc_box_sel,mapping=aes(x=DateTime,y=hypo_processing_mgL,color="Hypo_Processing"))

ggarrange(hypo_flows,hypo_change,hypo_process,nrow=3,ncol=1)

ggsave("./Fig_Output/Hypo_DOC_Processes.jpg",width=10,height=12,units="in",dpi=320)

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

# Plot!
ggplot(flora_depths_3)+
  geom_line(mapping=aes(x=DateTime,y=epi_ugL,color="Epi"))+
  geom_point(mapping=aes(x=DateTime,y=epi_ugL,color="Epi"))+
  geom_line(mapping=aes(x=DateTime,y=hypo_ugL,color="Hypo"))+
  geom_point(mapping=aes(x=DateTime,y=hypo_ugL,color="Hypo"))

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

# Calculate VW Epi and VW Hypo
ctd_50_epi <- ctd_50_depths_2 %>% 
  filter(depth %in% c(0.1,1.6,3.8)) %>% 
  filter(time >= as.POSIXct("2015-01-01")) %>% 
  select(time,Temp_C,DO_mgL,calc_DO_pSat,depth) %>% 
  pivot_wider(names_from = depth, values_from = c(Temp_C,DO_mgL,calc_DO_pSat), values_fil = NA, values_fn = mean)

ctd_50_epi <- na.omit(ctd_50_epi)

ctd_50_epi <- ctd_50_epi %>% 
  mutate(VW_temp_c = ((Temp_C_0.1*1.38*10^5)+(Temp_C_1.6*8.91*10^4)+(Temp_C_3.8*5.96*10^4))/((1.38*10^5)+(8.91*10^4)+(5.96*10^4))) %>% 
  mutate(VW_DO_mgL = ((DO_mgL_0.1*1.38*10^5)+(DO_mgL_1.6*8.91*10^4)+(DO_mgL_3.8*5.96*10^4))/((1.38*10^5)+(8.91*10^4)+(5.96*10^4))) %>% 
  mutate(VW_pSat_DO = ((calc_DO_pSat_0.1*1.38*10^5)+(calc_DO_pSat_1.6*8.91*10^4)+(calc_DO_pSat_3.8*5.96*10^4))/((1.38*10^5)+(8.91*10^4)+(5.96*10^4))) %>% 
  mutate(month = month(time))

ctd_50_hypo <- ctd_50_depths_2 %>% 
  filter(depth %in% c(5.0,6.2,8.0,9.0)) %>% 
  filter(time >= as.POSIXct("2015-01-01")) %>% 
  select(time,Temp_C,DO_mgL,calc_DO_pSat,depth) %>% 
  pivot_wider(names_from = depth, values_from = c(Temp_C,DO_mgL,calc_DO_pSat), values_fil = NA, values_fn = mean)

ctd_50_hypo <- na.omit(ctd_50_hypo)

ctd_50_hypo <- ctd_50_hypo %>% 
  mutate(VW_temp_c = ((Temp_C_5*4.02*10^4)+(Temp_C_6.2*1.40*10^4)+(Temp_C_8*1.40*10^4)+(Temp_C_9*1.95*10^3))/((4.02*10^4)+(1.40*10^4)+(1.40*10^4)+(1.95*10^3))) %>% 
  mutate(VW_DO_mgL = ((DO_mgL_5*4.02*10^4)+(DO_mgL_6.2*1.40*10^4)+(DO_mgL_8*1.40*10^4)+(DO_mgL_9*1.95*10^3))/((4.02*10^4)+(1.40*10^4)+(1.40*10^4)+(1.95*10^3))) %>% 
  mutate(VW_pSat_DO = ((calc_DO_pSat_5*4.02*10^4)+(calc_DO_pSat_6.2*1.40*10^4)+(calc_DO_pSat_8*1.40*10^4)+(calc_DO_pSat_9*1.95*10^3))/((4.02*10^4)+(1.40*10^4)+(1.40*10^4)+(1.95*10^3))) %>% 
  mutate(month = month(time))

# Plot for now: VW Temp
ggplot()+
  geom_point(ctd_50_epi,mapping=aes(x=time,y=VW_temp_c,color="VW Epi"))+
  geom_line(ctd_50_epi,mapping=aes(x=time,y=VW_temp_c,color="VW Epi"))+
  geom_point(ctd_50_hypo,mapping=aes(x=time,y=VW_temp_c,color="VW Hypo"))+
  geom_line(ctd_50_hypo,mapping=aes(x=time,y=VW_temp_c,color="VW Hypo"))+
  theme_classic(base_size=15)

# Plot VW DO_Sat
ggplot()+
  geom_point(ctd_50_epi,mapping=aes(x=time,y=VW_pSat_DO,color="VW Epi"))+
  geom_line(ctd_50_epi,mapping=aes(x=time,y=VW_pSat_DO,color="VW Epi"))+
  geom_point(ctd_50_hypo,mapping=aes(x=time,y=VW_pSat_DO,color="VW Hypo"))+
  geom_line(ctd_50_hypo,mapping=aes(x=time,y=VW_pSat_DO,color="VW Hypo"))+
  theme_classic(base_size=15)

# DO_mgL
ggplot()+
  geom_point(ctd_50_epi,mapping=aes(x=time,y=VW_DO_mgL,color="VW Epi"))+
  geom_line(ctd_50_epi,mapping=aes(x=time,y=VW_DO_mgL,color="VW Epi"))+
  geom_point(ctd_50_hypo,mapping=aes(x=time,y=VW_DO_mgL,color="VW Hypo"))+
  geom_line(ctd_50_hypo,mapping=aes(x=time,y=VW_DO_mgL,color="VW Hypo"))+
  theme_classic(base_size=15)