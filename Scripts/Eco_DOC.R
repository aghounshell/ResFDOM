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

### !!! NEED TO CALCULATE VW EPI AND HYPO FLORA !!! ###


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
