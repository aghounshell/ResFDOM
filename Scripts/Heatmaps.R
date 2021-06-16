### Script to plot heatmaps for DO and Temp from 2019-2020
# 13 May 2021, A Hounshell
# Updated: 16 Jue 2021

wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(akima,dplyr,ggplot2,tidyverse,reshape2,gridExtra,grid,colorRamps,RColorBrewer,lubridate,ggpubr)

# CTD and YSI casts - combine together for most complete time-period
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

## Select dates from 2018-2020
fcr_all_2 <- fcr_all %>% 
  filter(time>as.POSIXct("2018-01-01")&time<as.POSIXct("2020-12-31"))

# Filter out depths in the CTD cast that are closest to these specified values.
df.final<-data.frame()
ctd1<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 0.1)))
ctd2<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 0.4)))
ctd3<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 0.7)))
ctd4<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 1)))
ctd5<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 1.3)))
ctd6<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 1.6)))
ctd7<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 1.9)))
ctd8<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 2.3)))
ctd9<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 2.6)))
ctd10<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 2.9)))
ctd11<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 3.2)))
ctd12<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 3.5)))
ctd13<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 3.8)))
ctd14<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 4.1)))
ctd15<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 4.4)))
ctd16<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 4.7)))
ctd17<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 5)))
ctd18<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 5.3)))
ctd19<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 5.6)))
ctd20<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 5.9)))
ctd21<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 6.2)))
ctd22<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 6.5)))
ctd23<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 6.8)))
ctd24<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 7.1)))
ctd25<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 7.4)))
ctd26<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 7.7)))
ctd27<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 8)))
ctd28<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 8.3)))
ctd29<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 8.7)))
ctd30<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 9)))
ctd31<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 9.3)))
ctd32<-fcr_all_2 %>% group_by(time) %>% slice(which.min(abs(as.numeric(depth) - 9.6)))

# Bind each of the data layers together.
df.final = rbind(ctd1,ctd2,ctd3,ctd4,ctd5,ctd6,ctd7,ctd8,ctd9,ctd10,ctd11,ctd12,ctd13,ctd14,ctd15,ctd16,ctd17,ctd18,ctd19,
                 ctd20,ctd21,ctd22,ctd23,ctd24,ctd25,ctd26,ctd27,ctd28,ctd29,ctd30,ctd31,ctd32)

# Re-arrange the data frame by date
ctd_all <- arrange(df.final, time)

# Round each extracted depth to the nearest 10th. 
ctd_all$depth <- round(as.numeric(ctd_all$depth), digits = 0.5)
ctd_all <- ctd_all %>% group_by(time,depth) %>% summarise_each(funs(mean))
ctd_all$DOY <- yday(ctd_all$time)


# Select and make each CTD variable a separate dataframe
temp_18 <- ctd_all %>% 
  select(time, depth, Temp_C,DOY) %>% 
  filter(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31")) %>% 
  arrange(Temp_C,depth)
temp_18 <- na.omit(temp_18)

temp_19 <- ctd_all %>% 
  select(time, depth, Temp_C,DOY) %>% 
  filter(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31")) %>% 
  arrange(Temp_C,depth)
temp_19 <- na.omit(temp_19)

temp_20 <- ctd_all %>% 
  select(time, depth, Temp_C,DOY) %>% 
  filter(time >= as.POSIXct("2020-01-01") & time <= as.POSIXct("2020-12-31")) %>% 
  arrange(Temp_C,depth)
temp_20 <- na.omit(temp_20)

# Interpolate Temp
interp_temp_18 <- interp(x=temp_18$DOY, y = temp_18$depth, z = temp_18$Temp_C,
                         xo = seq(min(temp_18$DOY), max(temp_18$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                            extrap = F, linear = T, duplicate = "strip")
interp_temp_18 <- interp2xyz(interp_temp_18, data.frame=T)

interp_temp_19 <- interp(x=temp_19$DOY, y = temp_19$depth, z = temp_19$Temp_C,
                         xo = seq(min(temp_19$DOY), max(temp_19$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                         extrap = F, linear = T, duplicate = "strip")
interp_temp_19 <- interp2xyz(interp_temp_19, data.frame=T)

interp_temp_20 <- interp(x=temp_20$DOY, y = temp_20$depth, z = temp_20$Temp_C,
                         xo = seq(min(temp_20$DOY), max(temp_20$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                         extrap = F, linear = T, duplicate = "strip")
interp_temp_20 <- interp2xyz(interp_temp_20, data.frame=T)

# Plot
interp_temp_18$date <- as.Date(interp_temp_18$x,origin="2018-01-01")
interp_temp_19$date <- as.Date(interp_temp_19$x,origin="2019-01-01")
interp_temp_20$date <- as.Date(interp_temp_20$x,origin="2020-01-01")

## DO
DO_18 <- ctd_all %>% 
  select(time, depth, DO_mgL,DOY) %>% 
  filter(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31")) %>% 
  arrange(DO_mgL,depth)
DO_18 <- na.omit(DO_18)

DO_19 <- ctd_all %>% 
  select(time, depth, DO_mgL,DOY) %>% 
  filter(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31")) %>% 
  arrange(DO_mgL,depth)
DO_19 <- na.omit(DO_19)

DO_20 <- ctd_all %>% 
  select(time, depth, DO_mgL,DOY) %>% 
  filter(time >= as.POSIXct("2020-01-01") & time <= as.POSIXct("2020-12-31")) %>% 
  arrange(DO_mgL,depth)
DO_20 <- na.omit(DO_20)

# Interpolate DO
interp_DO_18 <- interp(x=DO_18$DOY, y = DO_18$depth, z = DO_18$DO_mgL,
                         xo = seq(min(DO_18$DOY), max(DO_18$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                         extrap = F, linear = T, duplicate = "strip")
interp_DO_18 <- interp2xyz(interp_DO_18, data.frame=T)

interp_DO_19 <- interp(x=DO_19$DOY, y = DO_19$depth, z = DO_19$DO_mgL,
                         xo = seq(min(DO_19$DOY), max(DO_19$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                         extrap = F, linear = T, duplicate = "strip")
interp_DO_19 <- interp2xyz(interp_DO_19, data.frame=T)

interp_DO_20 <- interp(x=DO_20$DOY, y = DO_20$depth, z = DO_20$DO_mgL,
                         xo = seq(min(DO_20$DOY), max(DO_20$DOY), by = .1), 
                         yo = seq(0.1, 9.6, by = 0.01),
                         extrap = F, linear = T, duplicate = "strip")
interp_DO_20 <- interp2xyz(interp_DO_20, data.frame=T)

# Plot
interp_DO_18$date <- as.Date(interp_DO_18$x,origin="2018-01-01")
interp_DO_19$date <- as.Date(interp_DO_19$x,origin="2019-01-01")
interp_DO_20$date <- as.Date(interp_DO_20$x,origin="2020-01-01")

# All the plots!
temp_18 <- ggplot(interp_temp_18, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2018-10-21"),linetype="dotted",colour="black",size=.7)+ #Turnover FCR
  geom_vline(xintercept = as.Date("2018-04-23"),colour="black",size=.7)+
  annotate("text",x=as.Date("2018-04-16"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2018-07-30"),colour="black",size=0.7)+
  annotate("text",x=as.Date("2018-07-23"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2019", y = "Depth (m)", fill=expression(''*~degree*C*''))+
  scale_fill_gradientn(colours = blue2green2red(60), na.value="gray")+
  xlim(as.Date("2018-03-01"),("2018-12-07"))+
  theme_classic(base_size = 15)

temp_19 <- ggplot(interp_temp_19, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2019-10-23"),linetype="dotted",size=0.7)+ #Turnover FCR
  geom_vline(xintercept = as.Date("2019-06-03"),size=0.7)+
  annotate("text",x=as.Date("2019-05-27"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-06-17"),size=0.7)+
  annotate("text",x=as.Date("2019-06-10"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-07-08"),size=0.7)+
  annotate("text",x=as.Date("2019-07-01"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-07-19"),size=0.7)+
  annotate("text",x=as.Date("2019-07-13"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-08-05"),size=0.7)+
  annotate("text",x=as.Date("2019-07-30"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-08-19"),size=0.7)+
  annotate("text",x=as.Date("2019-08-12"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-09-02"),size=0.7)+
  annotate("text",x=as.Date("2019-08-26"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-12-01"),size=0.7)+
  annotate("text",x=as.Date("2019-11-23"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2019", y = "Depth (m)", fill=expression(''*~degree*C*''))+
  scale_fill_gradientn(colours = blue2green2red(60), na.value="gray")+
  xlim(as.Date("2019-03-01"),("2019-12-07"))+
  theme_classic(base_size = 15)
  
temp_20 <- ggplot(interp_temp_20, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2020-11-01"),linetype="dotted",size=0.7)+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.Date("2020-06-29"),size=0.7)+
  annotate("text",x=as.Date("2020-06-22"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2020-09-11"),size=0.7)+
  annotate("text",x=as.Date("2020-09-02"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2020-09-25"),size=0.7)+
  annotate("text",x=as.Date("2020-09-19"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2020-12-02"),size=0.7)+
  annotate("text",x=as.Date("2020-11-25"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2020", y = "Depth (m)", fill=expression(''*~degree*C*''))+
  scale_fill_gradientn(colours = blue2green2red(60), na.value="gray")+
  xlim(as.Date("2020-03-01"),("2020-12-07"))+
  theme_classic(base_size=15)

ggarrange(temp_18,temp_19,temp_20,ncol=1,nrow=3,common.legend=TRUE,legend="right",labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Temp_heatmaps.png",width=10,height=12,units="in",dpi=320)

DO_18 <- ggplot(interp_DO_18, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2018-10-21"),linetype="dotted",colour="black",size=.7)+ #Turnover FCR
  geom_vline(xintercept = as.Date("2018-04-23"),colour="black",size=.7)+
  annotate("text",x=as.Date("2018-04-16"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2018-07-30"),colour="black",size=0.7)+
  annotate("text",x=as.Date("2018-07-23"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2018", y = "Depth (m)", fill=expression("DO (mg L"^-1*")"))+
  scale_fill_gradientn(colours = rev(blue2green2red(60)), na.value="gray")+
  xlim(as.Date("2018-03-01"),("2018-12-07"))+
  theme_classic(base_size = 15)

DO_19 <- ggplot(interp_DO_19, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2019-10-23"),linetype="dotted",size=0.7)+ #Turnover FCR
  geom_vline(xintercept = as.Date("2019-06-03"),size=0.7)+
  annotate("text",x=as.Date("2019-05-27"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-06-17"),size=0.7)+
  annotate("text",x=as.Date("2019-06-10"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-07-08"),size=0.7)+
  annotate("text",x=as.Date("2019-07-01"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-07-19"),size=0.7)+
  annotate("text",x=as.Date("2019-07-13"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-08-05"),size=0.7)+
  annotate("text",x=as.Date("2019-07-30"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-08-19"),size=0.7)+
  annotate("text",x=as.Date("2019-08-12"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2019-09-02"),size=0.7)+
  annotate("text",x=as.Date("2019-08-26"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2019-12-01"),size=0.7)+
  annotate("text",x=as.Date("2019-11-23"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2019", y = "Depth (m)", fill=expression("DO (mg L"^-1*")"))+
  scale_fill_gradientn(colours = rev(blue2green2red(60)), na.value="gray")+
  xlim(as.Date("2019-03-01"),("2019-12-07"))+
  theme_classic(base_size = 15)

DO_20 <- ggplot(interp_DO_20, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_vline(xintercept = as.Date("2020-11-01"),linetype="dotted",size=0.7)+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.Date("2020-06-29"),size=0.7)+
  annotate("text",x=as.Date("2020-06-22"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2020-09-11"),size=0.7)+
  annotate("text",x=as.Date("2020-09-02"),y=-0.3,label="Off",size=5)+
  geom_vline(xintercept = as.Date("2020-09-25"),size=0.7)+
  annotate("text",x=as.Date("2020-09-19"),y=-0.3,label="On",size=5)+
  geom_vline(xintercept = as.Date("2020-12-02"),size=0.7)+
  annotate("text",x=as.Date("2020-11-25"),y=-0.3,label="Off",size=5)+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 8, linetype="dashed", colour="white",size=1)+
  geom_hline(yintercept = 9, linetype="dashed", colour="white",size=1)+
  labs(x = "2020", y = "Depth (m)", fill=expression("DO (mg L"^-1*")"))+
  scale_fill_gradientn(colours = rev(blue2green2red(60)), na.value="gray")+
  xlim(as.Date("2020-03-01"),("2020-12-07"))+
  theme_classic(base_size=15)

ggarrange(DO_18,DO_19,DO_20,ncol=1,nrow=3,common.legend=TRUE,legend="right",labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/DO_heatmaps.png",width=10,height=12,units="in",dpi=320)



### OLD CODE ----

# Heatmaps in R ----



# When ready for 'fine' interpolation
#interp_temp <- interp(x=temp$DOY, y = temp$Depth_m, z = temp$Temp_C,
#                      xo = seq(min(temp$DOY), max(temp$DOY), by = .1), 
#                      yo = seq(0.1, 9.6, by = 0.01),
#                      extrap = F, linear = T, duplicate = "strip")
#interp_temp <- interp2xyz(interp_temp, data.frame=T)

# Formatting for Matlab ----
# Go from long to wide format
temp_wide <- temp %>% 
  pivot_wider(names_from = depth,values_from = Temp_C,names_prefix = "Depth_",values_fn = mean) %>% 
  arrange(time) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d")))

do_wide <- do %>% 
  pivot_wider(names_from = depth,values_from = DO_mgL,names_prefix = "Depth_",values_fn = mean) %>% 
  arrange(time)%>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d")))

write_csv(temp_wide,"./Data/Heatmap_Temp.csv")
write_csv(do_wide,"./Data/Heatmap_DO.csv")
interp_temp_18$date <- as.Date(interp_temp_18$x,origin="2018-01-01")