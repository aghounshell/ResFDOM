### Script to plot heatmaps for DO and Temp from 2015-2020
# 13 May 2021, A Hounshell

wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(akima,dplyr,ggplot2,tidyverse,reshape2,gridExtra,grid,colorRamps,RColorBrewer,lubridate)

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

## Select dates from 2015-2020
fcr_all_2 <- fcr_all %>% 
  filter(time>as.POSIXct("2015-01-01")&time<as.POSIXct("2020-12-31"))

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
#ctd_all$DOY <- yday(ctd_all$time)


# Select and make each CTD variable a separate dataframe
temp <- select(ctd_all, time, depth, Temp_C)
temp <- na.omit(temp)
temp <- arrange(temp,depth)
do <- select(ctd_all, time, depth, DO_mgL)
do <- na.omit(do)
do <- arrange(do,depth)

# Formatting for Matlab ----
# Go from long to wide format
temp_wide <- temp %>% 
  pivot_wider(names_from = depth,values_from = Temp_C,names_prefix = "Depth_",values_fn = mean) %>% 
  arrange(time)

do_wide <- do %>% 
  pivot_wider(names_from = depth,values_from = DO_mgL,names_prefix = "Depth_",values_fn = mean) %>% 
  arrange(time)

write_csv(temp_wide,"./Data/Heatmap_Temp.csv")
write_csv(do_wide,"./Data/Heatmap_DO.csv")

# Ignore for now: heatmaps in R ----

# Interpolate Temp
interp_temp_rough <- interp(x=temp$time, y = temp$depth, z = temp$Temp_C,
                      xo = seq(min(temp$time), max(temp$time), by = 100), 
                      yo = seq(0.1, 9.6, by = 0.3),
                      extrap = F, linear = T, duplicate = "strip")
interp_temp_rough <- interp2xyz(interp_temp_rough, data.frame=T)

# Plot
interp_temp_rough$date <- as.POSIXct(interp_temp_rough$x)

# Temp
ggplot(interp_temp_rough, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 1.6, linetype="dashed", colour="white")+
  geom_hline(yintercept = 3.8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 6.2, linetype="dashed", colour="white")+
  geom_hline(yintercept = 8, linetype="dashed", colour="white")+
  geom_hline(yintercept = 9, linetype="dashed", colour="white")+
  labs(x = "Date", y = "Depth (m)", fill=expression(''*~degree*C*''))+
  scale_fill_gradientn(colours = blue2green2red(60), na.value="gray")+
  theme_classic()

# When ready for 'fine' interpolation
#interp_temp <- interp(x=temp$DOY, y = temp$Depth_m, z = temp$Temp_C,
#                      xo = seq(min(temp$DOY), max(temp$DOY), by = .1), 
#                      yo = seq(0.1, 9.6, by = 0.01),
#                      extrap = F, linear = T, duplicate = "strip")
#interp_temp <- interp2xyz(interp_temp, data.frame=T)
