### Script to merge and plot CTD and YSI data from 2019
### A Hounshell, 29 Sep 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in CTD and YSI data
#need to import CTD observations from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/10/2461524a7da8f1906bfc3806d594f94c" 
infile1 <- paste0(getwd(),"/Data/CTD_final_2013_2019.csv")
download.file(inUrl1,infile1,method="curl")

#read in CTD temp file from EDI to create field file, but first need to subset CTD data per each day to depths
ctd<-read_csv('./Data/CTD_final_2013_2019.csv',col_types = cols(Reservoir = col_character(),Date = col_character(),pH=col_double())) %>% #read in observed CTD data, which has multiple casts on the same day (problematic for comparison)
  filter(Reservoir=="FCR") %>%
  filter(Site==50) %>%
  rename(time=Date, depth=Depth_m, temp=Temp_C, DO=DO_mgL, chla = Chla_ugL) %>%
  mutate(time = as.POSIXct(strptime(time, "%m/%d/%Y", tz="EST"))) %>%
  filter(time > as.POSIXct("2018-12-31") & time < as.POSIXct("2020-01-01")) %>% 
  select(time, depth, temp, DO, chla, pH) %>%
  na.omit() 

# Check pH variability
ctd_1 <- ctd %>%  
  filter(depth>=0 & depth<0.2) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="epi")

ctd_2 <- ctd %>% 
  filter(depth>=4.9 & depth<5.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="meta")

ctd_3 <- ctd %>% 
  filter(depth>=8.9 & depth<9.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="hypo")

# Plot pH
ggplot()+
  geom_line(ctd_1,mapping=aes(time,pH,color="epi"),size=1)+
  geom_point(ctd_1,mapping=aes(time,pH,color="epi"),size=2)+
  geom_line(ctd_2,mapping=aes(time,pH,color="meta"),size=1)+
  geom_point(ctd_2,mapping=aes(time,pH,color="meta"),size=2)+
  geom_line(ctd_3,mapping=aes(time,pH,color="hypo"),size=1)+
  geom_point(ctd_3,mapping=aes(time,pH,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

# Import YSI observations from EDI
inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/198/7/25b5e8b7f4291614d5c6d959a08148d8"
infile1 <- paste0(getwd(),"/Data/YSI_PAR_profiles_2013-2019.csv")
download.file(inUrl1,infile1,method="curl")

ysi <- read_csv('./Data/YSI_PAR_profiles_2013-2019.csv') %>% 
  filter(Reservoir=="FCR") %>% 
  filter(Site==50) %>% 
  rename(time=DateTime,depth=Depth_m,temp=Temp_C,DO=DO_mgL) %>% 
  mutate(time = as.POSIXct(strptime(time, "%m/%d/%y",tz='EST'))) %>% 
  filter(time > as.POSIXct("2018-12-31") & time < as.POSIXct("2020-01-01")) %>% 
  select(time,depth,temp,DO) %>% 
  na.omit()

# Select unique dates from both CTD and YSI casts
ysi_date_list <- as.data.frame(unique(as.Date(ysi$time)))
names(ysi_date_list)[1] <- "time"
ysi_date_list$ysi_bvr <- rep(-99,length(ysi_date_list$time))

ctd_date_list <- as.data.frame(unique(as.Date(ctd$time)))
names(ctd_date_list)[1] <- "time"
ctd_date_list$ctd_bvr <- rep(-99,length(ctd_date_list$time))

# Combine Unique dates list by date
fcr_dates <- merge(ysi_date_list, ctd_date_list, by="time", all.x=TRUE, all.y=TRUE)

### Merge data CTD and YSI datasets for BVR
fcr_merge <- merge(ctd, ysi, by="time", all.x=TRUE, all.y=TRUE)

# Find where there are Na values in the CTD data: need to do it for each column
ctd_fcr_na <- is.na(fcr_merge$depth.x)
fcr_merge$depth.x[ctd_fcr_na] <- fcr_merge$depth.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$temp.x)
fcr_merge$temp.x[ctd_fcr_na] <- fcr_merge$temp.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$DO.x)
fcr_merge$DO.x[ctd_fcr_na] <- fcr_merge$DO.y[ctd_fcr_na]

fcr_all <- fcr_merge %>% select(time,depth.x,temp.x,DO.x,chla) %>% 
  rename(depth=depth.x,temp=temp.x,DO=DO.x)

fcr_date_list <- as.data.frame(unique(as.Date(fcr_all$time)))

## Average across date and depth
fcr_all <- fcr_all %>% group_by(time,depth) %>% summarize_all(funs(mean))

# Save as .csv
write.csv(fcr_all, "./Data/YSICTD_Merge.csv", row.names = F)

########################################################
# Then average around 0.1 m, 5.0 m, and 9.0 m
ctd_1 <- fcr_all %>%  
  filter(depth>=0 & depth<0.2) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="epi")

ctd_2 <- fcr_all %>% 
  filter(depth>=4.9 & depth<5.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="meta")

ctd_3 <- fcr_all %>% 
  filter(depth>=8.9 & depth<9.1) %>% 
  group_by(time) %>% summarize_all(funs(mean)) %>% arrange(time) %>% 
  mutate(grouping="hypo")

# Plot DO
ggplot()+
  geom_line(ctd_1,mapping=aes(time,DO,color="epi"),size=1)+
  geom_point(ctd_1,mapping=aes(time,DO,color="epi"),size=2)+
  geom_line(ctd_2,mapping=aes(time,DO,color="meta"),size=1)+
  geom_point(ctd_2,mapping=aes(time,DO,color="meta"),size=2)+
  geom_line(ctd_3,mapping=aes(time,DO,color="hypo"),size=1)+
  geom_point(ctd_3,mapping=aes(time,DO,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

# Plot Temp
ggplot()+
  geom_line(ctd_1,mapping=aes(time,temp,color="epi"),size=1)+
  geom_point(ctd_1,mapping=aes(time,temp,color="epi"),size=2)+
  geom_line(ctd_2,mapping=aes(time,temp,color="meta"),size=1)+
  geom_point(ctd_2,mapping=aes(time,temp,color="meta"),size=2)+
  geom_line(ctd_3,mapping=aes(time,temp,color="hypo"),size=1)+
  geom_point(ctd_3,mapping=aes(time,temp,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

# Plot Chla (for the heck of it!)
ggplot()+
  geom_line(ctd_1,mapping=aes(time,chla,color="epi"),size=1)+
  geom_point(ctd_1,mapping=aes(time,chla,color="epi"),size=2)+
  geom_line(ctd_2,mapping=aes(time,chla,color="meta"),size=1)+
  geom_point(ctd_2,mapping=aes(time,chla,color="meta"),size=2)+
  geom_line(ctd_3,mapping=aes(time,chla,color="hypo"),size=1)+
  geom_point(ctd_3,mapping=aes(time,chla,color="hypo"),size=2)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('epi','meta','hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)


##########################################################
# Time for inflow data!!!
inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/202/6/96bdffa73741ec6b43a98f2c5d15daeb"
infile1 <- paste0(getwd(),"/Data/Inflow.csv")
download.file(inUrl1,infile1,method="curl")

inflow <- read_csv('./Data/Inflow.csv') %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M",tz='EST'))) %>% 
  filter(DateTime > as.POSIXct("2019-01-01") & DateTime < as.POSIXct("2020-01-01")) 

# Plot WVWA flow
ggplot(inflow,mapping=aes(DateTime,WVWA_Flow_cms))+
  geom_line(size=1)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  theme_classic(base_size=15)
