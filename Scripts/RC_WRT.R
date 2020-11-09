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
fcr_1 <- inflow %>% filter(Reservoir=="FCR" & Site=="1") %>% select(Date,Site,Flow_cms)
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
bvr_wl <- read.csv("C:/Users/ahoun/OneDrive/Desktop/BVR-GLM/BVR-GLM/Data_Output/09Nov20_BVR_WaterLevelDaily.csv")
bvr_wl$Date <-as.POSIXct(strptime(bvr_wl$Date, "%Y-%m-%d", tz = "EST"))

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

# Load in all continuum data (collated by WMW)
data <- read.csv("./Data/continuum_ww.csv")
data$DateTime <-as.POSIXct(strptime(data$DateTime, "%Y-%m-%d", tz = "EST"))

# Separate out by inflows/outflows
fcr_1_data <- data %>% 
  filter(Reservoir=="FCR" & Site == "1") %>% 
  select(DateTime,TN_ugL:Flow_cms)
fcr_99_data <- data %>% filter(Reservoir=="FCR" & Site == "99")%>% 
  select(DateTime,TN_ugL:Flow_cms)
fcr_200_data <- data %>% filter(Reservoir=="FCR" & Site == "200")%>% 
  select(DateTime,TN_ugL:Flow_cms)
fcr_50_data <- data %>% filter(Reservoir=="FCR" & Site == "50")%>% 
  select(DateTime,TN_ugL:Flow_cms)
fcr_102_data <- data %>% filter(Reservoir=="FCR" & Site == "102")%>% 
  select(DateTime,TN_ugL:Flow_cms)

bvr_100_data <- data %>% filter(Reservoir=="BVR" & Site == "100")%>% 
  select(DateTime,TN_ugL:Flow_cms)
bvr_200_data <- data %>% filter(Reservoir=="BVR" & Site == "200")%>% 
  select(DateTime,TN_ugL:Flow_cms)
bvr_50_data <- data %>% filter(Reservoir=="BVR" & Site == "50")%>% 
  select(DateTime,TN_ugL:Flow_cms)

# Calculate processing along the stream sites (from FCR 102 to FCR 99)
fcr_102_data <- fcr_102_data %>% 
  mutate(TN_load_mgs=TN_ugL*Flow_cms) %>% 
  mutate(TP_load_mgs=TP_ugL*Flow_cms) %>% 
  mutate(NH4_load_mgs=NH4_ugL*Flow_cms) %>% 
  mutate(NO3NO2_load_mgs=NO3NO2_ugL*Flow_cms) %>% 
  mutate(SRP_load_mgs=SRP_ugL*Flow_cms) %>% 
  mutate(DOC_load_gs=DOC_mgL*Flow_cms)

fcr_99_data <- fcr_99_data %>% 
  mutate(TN_load_mgs=TN_ugL*Flow_cms) %>% 
  mutate(TP_load_mgs=TP_ugL*Flow_cms) %>% 
  mutate(NH4_load_mgs=NH4_ugL*Flow_cms) %>% 
  mutate(NO3NO2_load_mgs=NO3NO2_ugL*Flow_cms) %>% 
  mutate(SRP_load_mgs=SRP_ugL*Flow_cms) %>% 
  mutate(DOC_load_gs=DOC_mgL*Flow_cms)

stream_p <- fcr_102_data %>% select(DateTime)
# (+) = production; (-) = consumption
stream_p$TN_mgs <- fcr_99_data$TN_load_mgs-fcr_102_data$TN_load_mgs
stream_p$TP_mgs <- fcr_99_data$TP_load_mgs-fcr_102_data$TP_load_mgs
stream_p$NH4_mgs <- fcr_99_data$NH4_load_mgs-fcr_102_data$NH4_load_mgs
stream_p$NO3NO2_mgs <- fcr_99_data$NO3NO2_load_mgs-fcr_102_data$NO3NO2_load_mgs
stream_p$SRP_mgs <- fcr_99_data$SRP_load_mgs-fcr_102_data$SRP_load_mgs
stream_p$DOC_gs <- fcr_99_data$DOC_load_gs-fcr_102_data$DOC_load_gs

fcr_200_data <- fcr_200_data %>% 
  mutate(TN_load_mgs=TN_ugL*Flow_cms) %>% 
  mutate(TP_load_mgs=TP_ugL*Flow_cms) %>% 
  mutate(NH4_load_mgs=NH4_ugL*Flow_cms) %>% 
  mutate(NO3NO2_load_mgs=NO3NO2_ugL*Flow_cms) %>% 
  mutate(SRP_load_mgs=SRP_ugL*Flow_cms) %>% 
  mutate(DOC_load_gs=DOC_mgL*Flow_cms)

fcr_1_data <- fcr_1_data %>% 
  mutate(TN_load_mgs=TN_ugL*Flow_cms) %>% 
  mutate(TP_load_mgs=TP_ugL*Flow_cms) %>% 
  mutate(NH4_load_mgs=NH4_ugL*Flow_cms) %>% 
  mutate(NO3NO2_load_mgs=NO3NO2_ugL*Flow_cms) %>% 
  mutate(SRP_load_mgs=SRP_ugL*Flow_cms) %>% 
  mutate(DOC_load_gs=DOC_mgL*Flow_cms)

fcr_res_p <- fcr_102_data %>% select(DateTime)
fcr_res_p$TN_mgs <- fcr_1_data$TN_load_mgs-fcr_99_data$TN_load_mgs-fcr_200_data$TN_load_mgs
fcr_res_p$TP_mgs <- fcr_1_data$TP_load_mgs-fcr_99_data$TP_load_mgs-fcr_200_data$TP_load_mgs
fcr_res_p$NH4_mgs <- fcr_1_data$NH4_load_mgs-fcr_99_data$NH4_load_mgs-fcr_200_data$NH4_load_mgs
fcr_res_p$NO3NO2_mgs <- fcr_1_data$NO3NO2_load_mgs-fcr_99_data$NO3NO2_load_mgs-fcr_200_data$NO3NO2_load_mgs
fcr_res_p$SRP_mgs <- fcr_1_data$SRP_load_mgs-fcr_99_data$SRP_load_mgs-fcr_200_data$SRP_load_mgs
fcr_res_p$DOC_gs <- fcr_1_data$DOC_load_gs-fcr_99_data$DOC_load_gs-fcr_200_data$DOC_load_gs
fcr_res_p$DOC_in <- fcr_99_data$DOC_load_gs+fcr_200_data$DOC_load_gs

ggplot()+
  geom_point(fcr_102_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 102"))+
  geom_line(fcr_102_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 102"))+
  geom_point(fcr_99_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 99"))+
  geom_line(fcr_99_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 99"))+
  geom_point(stream_p,mapping=aes(x=DateTime,y=TN_mgs,color="Stream Processing"))+
  geom_line(stream_p,mapping=aes(x=DateTime,y=TN_mgs,color="Stream Processing"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_200_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 200"))+
  geom_line(fcr_200_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 200"))+
  geom_point(fcr_99_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 99"))+
  geom_line(fcr_99_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 99"))+
  geom_point(fcr_1_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 1"))+
  geom_line(fcr_1_data,mapping=aes(x=DateTime,y=TN_load_mgs,color="FCR 1"))+
  geom_point(fcr_res_p,mapping=aes(x=DateTime,y=TN_mgs,color="FCR Processing"))+
  geom_line(fcr_res_p,mapping=aes(x=DateTime,y=TN_mgs,color="FCR Processing"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_102_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 102"))+
  geom_line(fcr_102_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 102"))+
  geom_point(fcr_99_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 99"))+
  geom_line(fcr_99_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 99"))+
  geom_point(stream_p,mapping=aes(x=DateTime,y=TP_mgs,color="Stream Processing"))+
  geom_line(stream_p,mapping=aes(x=DateTime,y=TP_mgs,color="Stream Processing"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_200_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 200"))+
  geom_line(fcr_200_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 200"))+
  geom_point(fcr_99_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 99"))+
  geom_line(fcr_99_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 99"))+
  geom_point(fcr_1_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 1"))+
  geom_line(fcr_1_data,mapping=aes(x=DateTime,y=TP_load_mgs,color="FCR 1"))+
  geom_point(fcr_res_p,mapping=aes(x=DateTime,y=TP_mgs,color="FCR Processing"))+
  geom_line(fcr_res_p,mapping=aes(x=DateTime,y=TP_mgs,color="FCR Processing"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_102_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 102"))+
  geom_line(fcr_102_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 102"))+
  geom_point(fcr_99_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 99"))+
  geom_line(fcr_99_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 99"))+
  geom_point(stream_p,mapping=aes(x=DateTime,y=DOC_gs,color="Stream Processing"))+
  geom_line(stream_p,mapping=aes(x=DateTime,y=DOC_gs,color="Stream Processing"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_res_p,mapping=aes(x=DateTime,y=DOC_in,color="FCR In"))+
  geom_line(fcr_res_p,mapping=aes(x=DateTime,y=DOC_in,color="FCR In"))+
  geom_point(fcr_1_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 1"))+
  geom_line(fcr_1_data,mapping=aes(x=DateTime,y=DOC_load_gs,color="FCR 1"))+
  geom_point(fcr_res_p,mapping=aes(x=DateTime,y=DOC_gs,color="FCR Processing"))+
  geom_line(fcr_res_p,mapping=aes(x=DateTime,y=DOC_gs,color="FCR Processing"))+
  theme_classic(base_size=15)

