### Script to calculate WRT for FCR and BVR on RC days
### 13 Oct 2020, A Hounshell

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in data from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/454/4/a18421fd2e95c15d6f97009d5fef3e59" 
infile1 <- paste0(getwd(),"/Data/RC_Inflow.csv")
download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/RC_Inflow.csv")

# Separate by site
fcr_99 <- inflow %>% filter(Reservoir=="FCR" & Site=="99") %>% select(Date,Site,Flow_cms)
fcr_200 <- inflow %>% filter(Reservoir=="FCR" & Site=="200") %>% select(Date,Site,Flow_cms)
bvr_100 <- inflow %>% filter(Reservoir=="BVR"& Site=="100") %>% select(Date,Site,Flow_cms)
bvr_200 <- inflow %>% filter(Reservoir=="BVR" & Site=="200") %>% select(Date,Site,Flow_cms)

# Combine by reservoir
fcr <- rbind.data.frame(fcr_99,fcr_200)
fcr <- fcr %>% group_by(Date) %>% summarise_all(funs(sum)) %>% filter(Site=="299")

bvr <- rbind.data.frame(bvr_100,bvr_200)
bvr <- bvr %>% group_by(Date) %>% summarise_all(funs(sum)) %>% filter(Site=="300")

# Calculate WRT
# FCR: Assume full pond volume (3.1E5 m3)
fcr <- fcr %>% mutate(wrt_d = (3.1E5/Flow_cms)*(1/60)*(1/60)*(1/24))

# BVR: Averaged water level for 2019 (10.91m) which corresponds to 8.9E5 m3
bvr <- bvr %>% mutate(wrt_d = (8.9E5/Flow_cms)*(1/60)*(1/60)*(1/24))

# Plot
ggplot()+
  geom_line(bvr,mapping=aes(x=Date,y=wrt_d,color='BVR',group=1),size=1)+
  geom_point(bvr,mapping=aes(x=Date,y=wrt_d,color='BVR'),size=2)+
  geom_line(fcr,mapping=aes(x=Date,y=wrt_d,color='FCR',group=1),size=1)+
  geom_point(fcr,mapping=aes(x=Date,y=wrt_d,color='FCR'),size=2)+
  theme_classic(base_size=15)

# Plot water level for BVR
bvr_wl <- read.csv("C:/Users/ahoun/OneDrive/Desktop/BVR-GLM/BVR-GLM/Data_Output/09Apr20_BVR_WaterLevelDailyVol_2.csv")
bvr_wl$Date <-as.POSIXct(strptime(bvr_wl$Date, "%m/%d/%Y", tz = "EST"))

# Plot
ggplot(bvr_wl,mapping=aes(x=Date,y=BVR_WaterLevel_m))+
  geom_line()+
  xlim(as.POSIXct("2019-01-01"),as.POSIXct("2019-12-31"))+
  geom_vline(xintercept = as.POSIXct("2019-04-29"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-05-30"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-06-20"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-07-18"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-07-24"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-08-22"),linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2019-10-04"),linetype="dashed")+
  theme_classic(base_size=15)
