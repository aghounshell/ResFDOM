# Script to calculate thermocline depth with rLakeAnalyzer
# A Hounshell, 30 May 2021

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate,zoo)

#devtools::install_github("GLEON/rLakeAnalyzer")

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

layer <- fcr_all %>% 
  rename(Date = time, Depth_m = depth)

### Formatting for LakeAnalyzer ----
df.final<-data.frame()

layer1<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.1)))
layer1$Depth_m <- 0.1
layer2<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.5)))
layer2$Depth_m <- 0.5
layer3<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1)))
layer3$Depth_m <- 1.0
layer4<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.5)))
layer4$Depth_m <- 1.5
layer5<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2)))
layer5$Depth_m <- 2.0
layer6<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.5)))
layer6$Depth_m <- 2.5
layer7<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3)))
layer7$Depth_m <- 3.0
layer8<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.5)))
layer8$Depth_m <- 3.5
layer9<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4)))
layer9$Depth_m <- 4.0
layer10<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.5)))
layer10$Depth_m <- 4.5
layer11<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5)))
layer11$Depth_m <- 5.0
layer12<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.5)))
layer12$Depth_m <- 5.5
layer13<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6)))
layer13$Depth_m <- 6.0
layer14<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.5)))
layer14$Depth_m <- 6.5
layer15<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7)))
layer15$Depth_m <- 7.0
layer16<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.5)))
layer16$Depth_m <- 7.5
layer17<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8)))
layer17$Depth_m <- 8.0
layer18<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.5)))
layer18$Depth_m <- 8.5
layer19<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9)))
layer19$Depth_m <- 9.0

df.final = rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7,layer8,layer9,layer10,layer11,layer12,layer13,
                 layer14,layer15,layer16,layer17,layer18,layer19)

fcr_layers <- arrange(df.final, Date)
fcr_layers$Depth_m <- round(as.numeric(fcr_layers$Depth_m), digits = .5)

fcr_layers_temp <- fcr_layers %>% select(Date,Depth_m,Temp_C) %>% group_by(Date,Depth_m) %>% summarise_each(funs(mean))

fcr_new <- fcr_layers_temp %>% spread(Depth_m,Temp_C)

fcr_new <- fcr_new[complete.cases(fcr_new),]

names(fcr_new)[1] <- "dateTime"
names(fcr_new)[2] <- "temp0.1"
names(fcr_new)[3] <- "temp0.5"
names(fcr_new)[4] <- "temp1.0"
names(fcr_new)[5] <- "temp1.5"
names(fcr_new)[6] <- "temp2.0"
names(fcr_new)[7] <- "temp2.5"
names(fcr_new)[8] <- "temp3.0"
names(fcr_new)[9] <- "temp3.5"
names(fcr_new)[10] <- "temp4.0"
names(fcr_new)[11] <- "temp4.5"
names(fcr_new)[12] <- "temp5.0"
names(fcr_new)[13] <- "temp5.5"
names(fcr_new)[14] <- "temp6.0"
names(fcr_new)[15] <- "temp6.5"
names(fcr_new)[16] <- "temp7.0"
names(fcr_new)[17] <- "temp7.5"
names(fcr_new)[18] <- "temp8.0"
names(fcr_new)[19] <- "temp8.5"
names(fcr_new)[20] <- "temp9.0"

# Export out for LakeAnalyzer in Matlab
write.csv(fcr_new,"./Data/LA_FCR.wrt")

# Re-import results from Lake Analyzer in Matlab
# Looking specifically at thermocline depth for the study period
# Load in results from Lake Analyzer in Matlab
la_results <- read.csv("./Data/20210603_LA_FCR_results.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2020-01-01"))

la_results_strat <- la_results %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(4,5,6,7,8,9,10))

# Plot to check?
ggplot(la_results)+
  geom_line(la_results_strat,mapping=aes(x=DateTime,y=-thermD,color="thermD"))+
  geom_point(la_results_strat,mapping=aes(x=DateTime,y=-thermD,color="thermD"))+
  geom_line(la_results_strat,mapping=aes(x=DateTime,y=-SthermD,color="SthermD"))+
  geom_point(la_results_strat,mapping=aes(x=DateTime,y=-SthermD,color="SthermD"))+
  theme_classic(base_size=15)

# Average across the stratified period?
thermo <- la_results_strat %>% 
  summarise_all(mean,na.rm=TRUE)

### DID NOT USE BELOW ----
# Could not figure out how to use LA in R... : (

wtr <- fcr_new %>% 
  ungroup() %>% 
  select(temp0.1:temp9.0)

### Calculate Thermocline depth ----
pacman::p_load(rLakeAnalyzer)

depths <- c(0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9,0)
depths=depths[1:length(depths)-1]

wtr <- as.numeric(unlist(wtr))

thermo <- thermo.depth(wtr,depths,seasonal=FALSE)


### OLD CODE ----

### Create a dataframe for CTD/YSI parameters at each sampling depth
depths <- seq(0.1, 9.0, by = 0.1)

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
  layer2<- q[q[, "depth"] == closest(q$depth,0.2),][1,]
  layer3<- q[q[, "depth"] == closest(q$depth,0.3),][1,]
  layer4<- q[q[, "depth"] == closest(q$depth,0.4),][1,]
  layer5<- q[q[, "depth"] == closest(q$depth,0.5),][1,]
  layer6<- q[q[, "depth"] == closest(q$depth,0.6),][1,]
  layer7<- q[q[, "depth"] == closest(q$depth,0.7),][1,]
  layer8 <- q[q[, "depth"] == closest(q$depth,0.8),][1,]
  layer9<- q[q[, "depth"] == closest(q$depth,0.9),][1,]
  layer10<- q[q[, "depth"] == closest(q$depth,1.0),][1,]
  layer11<- q[q[, "depth"] == closest(q$depth,1.1),][1,]
  layer12<- q[q[, "depth"] == closest(q$depth,1.2),][1,]
  layer13<- q[q[, "depth"] == closest(q$depth,1.3),][1,]
  layer14<- q[q[, "depth"] == closest(q$depth,1.4),][1,]
  layer15 <- q[q[, "depth"] == closest(q$depth,1.5),][1,]
  layer16<- q[q[, "depth"] == closest(q$depth,1.6),][1,]
  layer17<- q[q[, "depth"] == closest(q$depth,1.7),][1,]
  layer18<- q[q[, "depth"] == closest(q$depth,1.8),][1,]
  layer19<- q[q[, "depth"] == closest(q$depth,1.9),][1,]
  layer20<- q[q[, "depth"] == closest(q$depth,2.0),][1,]
  layer21<- q[q[, "depth"] == closest(q$depth,2.1),][1,]
  layer22 <- q[q[, "depth"] == closest(q$depth,2.2),][1,]
  layer23<- q[q[, "depth"] == closest(q$depth,2.3),][1,]
  layer24<- q[q[, "depth"] == closest(q$depth,2.4),][1,]
  layer25<- q[q[, "depth"] == closest(q$depth,2.5),][1,]
  layer26<- q[q[, "depth"] == closest(q$depth,2.6),][1,]
  layer27<- q[q[, "depth"] == closest(q$depth,2.7),][1,]
  layer28<- q[q[, "depth"] == closest(q$depth,2.8),][1,]
  layer29 <- q[q[, "depth"] == closest(q$depth,2.9),][1,]
  layer30<- q[q[, "depth"] == closest(q$depth,3.0),][1,]
  layer31<- q[q[, "depth"] == closest(q$depth,3.1),][1,]
  layer32<- q[q[, "depth"] == closest(q$depth,3.2),][1,]
  layer33<- q[q[, "depth"] == closest(q$depth,3.3),][1,]
  layer34<- q[q[, "depth"] == closest(q$depth,3.4),][1,]
  layer35<- q[q[, "depth"] == closest(q$depth,3.5),][1,]
  layer36 <- q[q[, "depth"] == closest(q$depth,3.6),][1,]
  layer37<- q[q[, "depth"] == closest(q$depth,3.7),][1,]
  layer38<- q[q[, "depth"] == closest(q$depth,3.8),][1,]
  layer39<- q[q[, "depth"] == closest(q$depth,3.9),][1,]
  layer40<- q[q[, "depth"] == closest(q$depth,4.0),][1,]
  layer41<- q[q[, "depth"] == closest(q$depth,4.1),][1,]
  layer42<- q[q[, "depth"] == closest(q$depth,4.2),][1,]
  layer43 <- q[q[, "depth"] == closest(q$depth,4.3),][1,]
  layer44<- q[q[, "depth"] == closest(q$depth,4.4),][1,]
  layer45<- q[q[, "depth"] == closest(q$depth,4.5),][1,]
  layer46<- q[q[, "depth"] == closest(q$depth,4.6),][1,]
  layer47<- q[q[, "depth"] == closest(q$depth,4.7),][1,]
  layer48<- q[q[, "depth"] == closest(q$depth,4.8),][1,]
  layer49<- q[q[, "depth"] == closest(q$depth,4.9),][1,]
  layer50 <- q[q[, "depth"] == closest(q$depth,5.0),][1,]
  layer51<- q[q[, "depth"] == closest(q$depth,5.1),][1,]
  layer52<- q[q[, "depth"] == closest(q$depth,5.2),][1,]
  layer53<- q[q[, "depth"] == closest(q$depth,5.3),][1,]
  layer54<- q[q[, "depth"] == closest(q$depth,5.4),][1,]
  layer55<- q[q[, "depth"] == closest(q$depth,5.5),][1,]
  layer56<- q[q[, "depth"] == closest(q$depth,5.6),][1,]
  layer57 <- q[q[, "depth"] == closest(q$depth,5.7),][1,]
  layer58<- q[q[, "depth"] == closest(q$depth,5.8),][1,]
  layer59<- q[q[, "depth"] == closest(q$depth,5.9),][1,]
  layer60<- q[q[, "depth"] == closest(q$depth,6.0),][1,]
  layer61<- q[q[, "depth"] == closest(q$depth,6.1),][1,]
  layer62<- q[q[, "depth"] == closest(q$depth,6.2),][1,]
  layer63<- q[q[, "depth"] == closest(q$depth,6.3),][1,]
  layer64 <- q[q[, "depth"] == closest(q$depth,6.4),][1,]
  layer65<- q[q[, "depth"] == closest(q$depth,6.5),][1,]
  layer66<- q[q[, "depth"] == closest(q$depth,6.6),][1,]
  layer67<- q[q[, "depth"] == closest(q$depth,6.7),][1,]
  layer68<- q[q[, "depth"] == closest(q$depth,6.8),][1,]
  layer69<- q[q[, "depth"] == closest(q$depth,6.9),][1,]
  layer70<- q[q[, "depth"] == closest(q$depth,7.0),][1,]
  layer71 <- q[q[, "depth"] == closest(q$depth,7.1),][1,]
  layer72<- q[q[, "depth"] == closest(q$depth,7.2),][1,]
  layer73<- q[q[, "depth"] == closest(q$depth,7.3),][1,]
  layer74<- q[q[, "depth"] == closest(q$depth,7.4),][1,]
  layer75<- q[q[, "depth"] == closest(q$depth,7.5),][1,]
  layer76<- q[q[, "depth"] == closest(q$depth,7.6),][1,]
  layer77<- q[q[, "depth"] == closest(q$depth,7.7),][1,]
  layer78 <- q[q[, "depth"] == closest(q$depth,7.8),][1,]
  layer79<- q[q[, "depth"] == closest(q$depth,7.9),][1,]
  layer80<- q[q[, "depth"] == closest(q$depth,8.0),][1,]
  layer81<- q[q[, "depth"] == closest(q$depth,8.1),][1,]
  layer82<- q[q[, "depth"] == closest(q$depth,8.2),][1,]
  layer83<- q[q[, "depth"] == closest(q$depth,8.3),][1,]
  layer84<- q[q[, "depth"] == closest(q$depth,8.4),][1,]
  layer85 <- q[q[, "depth"] == closest(q$depth,8.5),][1,]
  layer86<- q[q[, "depth"] == closest(q$depth,8.6),][1,]
  layer87<- q[q[, "depth"] == closest(q$depth,8.7),][1,]
  layer88<- q[q[, "depth"] == closest(q$depth,8.8),][1,]
  layer89<- q[q[, "depth"] == closest(q$depth,8.9),][1,]
  layer90<- q[q[, "depth"] == closest(q$depth,9.0),][1,]
  
  temp<-rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7,layer8,layer9,layer10,
              layer11,layer12,layer13,layer14,layer15,layer16,layer17,layer18,layer19,
              layer20,layer21,layer22,layer23,layer24,layer25,layer26,layer27,layer28,
              layer29,layer30,layer31,layer32,layer33,layer34,layer35,layer36,layer37,
              layer38,layer39,layer40,layer41,layer42,layer43,layer44,layer45,layer46,
              layer47,layer48,layer49,layer50,layer51,layer52,layer53,layer54,layer55,
              layer56,layer57,layer58,layer59,layer60,layer61,layer62,layer63,layer64,
              layer65,layer66,layer67,layer68,layer69,layer70,layer71,layer72,layer73,
              layer74,layer75,layer76,layer77,layer78,layer79,layer80,layer81,layer82,
              layer83,layer84,layer85,layer86,layer87,layer88,layer89,layer90)
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