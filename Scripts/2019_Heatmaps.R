### Script to create heatmaps for FCR; Summer 2019
### Start with DO and Temp

### A Hounshell, 20 Nov 2019; Following: FCR_heatmaps

# Rfile: FCR_Heatmap_2019

# Load libraries
pacman::p_load(akima,dplyr,ggplot2,tidyverse,reshape2,gridExtra,grid,colorRamps,RColorBrewer,lubridate)

# Load in merged CTD and YSI casts for the 2016-2017 period
# Will need to divide into Summer 2016 and Summer 2017
casts <- read_csv("C:/Users/ahoun/Dropbox/ResFDOM/ResFDOM/Data/20Nov19_CTD2019.csv")
casts$Date <- as.POSIXct(strptime(casts$Date, "%Y-%m-%d", tz = "EST"))

# Filter out depths in the CTD cast that are closest to these specified values.
df.final<-data.frame()
ctd1<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.1)))
ctd2<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.4)))
ctd3<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.7)))
ctd4<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1)))
ctd5<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.3)))
ctd6<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.6)))
ctd7<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.9)))
ctd8<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.3)))
ctd9<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.6)))
ctd10<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.9)))
ctd11<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.2)))
ctd12<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.5)))
ctd13<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.8)))
ctd14<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.1)))
ctd15<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.4)))
ctd16<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.7)))
ctd17<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5)))
ctd18<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.3)))
ctd19<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.6)))
ctd20<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.9)))
ctd21<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.2)))
ctd22<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.5)))
ctd23<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.8)))
ctd24<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.1)))
ctd25<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.4)))
ctd26<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.7)))
ctd27<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8)))
ctd28<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.3)))
ctd29<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.7)))
ctd30<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9)))
ctd31<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9.3)))
ctd32<-casts %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9.6)))

# Bind each of the data layers together.
df.final = rbind(ctd1,ctd2,ctd3,ctd4,ctd5,ctd6,ctd7,ctd8,ctd9,ctd10,ctd11,ctd12,ctd13,ctd14,ctd15,ctd16,ctd17,ctd18,ctd19,
                 ctd20,ctd21,ctd22,ctd23,ctd24,ctd25,ctd26,ctd27,ctd28,ctd29,ctd30,ctd31,ctd32)

# Re-arrange the data frame by date
ctd_19 <- arrange(df.final, Date)

# Round each extracted depth to the nearest 10th. 
ctd_19$Depth_m <- round(as.numeric(ctd_19$Depth_m), digits = 0.5)
ctd_19 <- ctd_19 %>% group_by(Date,Depth_m) %>% summarise_each(funs(mean))
ctd_19$DOY <- yday(ctd_19$Date)


# Select and make each CTD variable a separate dataframe
temp_19 <- select(ctd_19, Date, DOY, Depth_m, Temp_C)
temp_19 <- na.omit(temp_19)
do_19 <- select(ctd_19, Date, DOY, Depth_m, DO_mgL)
do_19 <- na.omit(do_19)

# Complete data interpolation for the heatmaps
# interative processes here

#temperature
interp_temp <- interp(x=temp_19$DOY, y = temp_19$Depth_m, z = temp_19$Temp_C,
                      xo = seq(min(temp_19$DOY), max(temp_19$DOY), by = .1), 
                      yo = seq(0.1, 9.6, by = 0.01),
                      extrap = F, linear = T, duplicate = "strip")
interp_temp <- interp2xyz(interp_temp, data.frame=T)

#dissolved oxygen
interp_do <- interp(x=do_19$DOY, y = do_19$Depth_m, z = do_19$DO_mgL,
                    xo = seq(min(do_19$DOY), max(do_19$DOY), by = .1), 
                    yo = seq(0.1, 10.2, by = 0.01),
                    extrap = F, linear = T, duplicate = "strip")
interp_do <- interp2xyz(interp_do, data.frame=T)

interp_temp$date <- as.Date(interp_temp$x,origin="2019-01-01")
interp_do$date <- as.Date(interp_do$x,origin="2019-01-01")

# Plot DO
ggplot(interp_do, aes(x=date, y=y))+
  geom_raster(aes(fill=z))+
  scale_y_reverse()+
  geom_hline(yintercept = 0.1, linetype="dashed", colour="white")+
  geom_hline(yintercept = 5, linetype="dashed", colour="white")+
  geom_hline(yintercept = 9, linetype="dashed", colour="white")+
  geom_vline(xintercept = as.Date("2019-05-27"), color="black")+
  geom_vline(xintercept = as.Date("2019-07-03"), color="black")+
  geom_vline(xintercept = as.Date("2019-07-15"), color="black")+
  geom_vline(xintercept = as.Date("2019-08-14"), color="black")+
  geom_point(aes(x = as.Date("2019-06-03"),y = 0),shape=25,color="black",fill="black",size=2)+ #Oxygen on
  geom_point(aes(x = as.Date("2019-06-17"),y = 0),shape=0,size=2)+ #Oxygen off
  geom_point(aes(x = as.Date("2019-07-08"),y = 0),shape=25,color="black",fill="black",size=2)+ # Oxygen on
  geom_point(aes(x = as.Date("2019-07-19"),y = 0),shape=0,size=2)+ #Oxygen off
  geom_point(aes(x = as.Date("2019-08-05"),y = 0),shape=25,color="black",fill="black",size=2)+ # Oxygen on
  geom_point(aes(x = as.Date("2019-08-19"),y = 0),shape=0,size=2)+ #Oxygen off
  geom_point(aes(x = as.Date("2019-09-02"),y = 0),shape=25,color="black",fill="black",size=2)+ # Oxygen on
  geom_point(aes(x = as.Date("2019-11-02"),y = 0),shape=15, size=2)+ #Turnover
  labs(x = "2019", y = "Depth (m)", fill=expression("DO (mg L"^-1*")"))+
  scale_fill_gradientn(colours = rev(blue2green2red(60)), na.value="gray")+
  theme_classic(base_size=15)
