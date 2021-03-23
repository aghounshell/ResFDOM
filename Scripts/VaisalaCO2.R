### Script to explore dissolved CO2 sensor deployed by Mark Johnson ###
# Following script from ABP: https://github.com/CareyLabVT/ManualDownloadsSCCData/blob/master/Scripts/Co2%20sensor.R
# 09 Mar 2021, A Hounshell
# Updated: 23 Mar 2021 to add FDOM data and plot both CO2 and FDOM data w/ ice on/off data for winter 2020/2021

pacman::p_load(tidyverse,ggplot2)

# NOTE: Will need to pull SCCManual Repo for most updated data!
mydir = "C:/Users/ahoun/Desktop/ManualDownloadsSCCData/FCR_CO2sensors/Vaisala Sensor"
myfiles = list.files(path=mydir, pattern="CR1000_CO2_Sensor_FCR_CO2*", full.names=TRUE)

#combine the files
#create an out.file for the combined data
out.file<-""
#for loop to combine the files
for(i in 1:length(myfiles)){
  file <- read.csv(myfiles[i], header=F, skip=4)
  out.file <- rbind(out.file, file)
}

#Naming the header because they have to be eliminated above to combine the files
colnames(out.file)=c("TIMESTAMP","RECORD","batt_volt_Min","PTemp","CO2_1_Avg","CO2_2_Avg")
#make the dates cooperate
out.file$TIMESTAMP <- as.POSIXct(out.file$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+4")
#add a month column to make sorting easier
out.file$Month=months(out.file$TIMESTAMP)
#change from character to Numeric
out.file$CO2_2_Avg=as.numeric(out.file$CO2_2_Avg)

#create the graph (just to see!)
ggplot(out.file, aes(x=TIMESTAMP, y=CO2_2_Avg, col=CO2_2_Avg)) +
  geom_point()

## Load in ice on/off data from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/456/3/ddb604975b844125b3bd542862ca0a3b" 
infile1 <- paste0(getwd(),"/Data/Ice_Data.csv")
download.file(inUrl1,infile1,method="curl")

ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))
ice_on <- ice %>% 
  filter(Date>"2020-01-01" & IceOn == 1)
ice_off <- ice %>% 
  filter(Date>"2020-01-01", IceOff == 1)

# Create graph to look at ice on/off
ggplot()+
  geom_point(out.file,mapping=aes(x=TIMESTAMP,y=CO2_2_Avg),size=0.1)+
  geom_line(out.file,mapping=aes(x=TIMESTAMP,y=CO2_2_Avg))+
  geom_vline(data = ice_on,mapping=aes(xintercept = Date), linetype = "dashed", color="cadetblue")+
  geom_vline(data = ice_off,mapping=aes(xintercept = Date), linetype = "dashed", color="blue")+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylim(0,2500)+
  theme_classic(base_size = 10)

# Sum CO2 to daily?
daily <- out.file %>% 
  mutate(Date = format(as.POSIXct(out.file$TIMESTAMP,"%Y-%m-%d"),format="%Y-%m-%d")) %>% 
  group_by(Date) %>% 
  summarize_all(funs(mean))
daily$Date <- as.POSIXct(strptime(daily$Date,"%Y-%m-%d"))

# Plot to see?
ggplot(daily,mapping=aes(x=Date,y=CO2_2_Avg))+
  geom_line()+
  geom_point()+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylim(0,1500)+
  theme_classic(base_size=10)

# Plot all together
ggplot()+
  geom_line(out.file,mapping=aes(x=TIMESTAMP,y=CO2_2_Avg,color="All"),color="grey")+
  geom_vline(data = ice_on,mapping=aes(xintercept = Date), linetype = "dashed", color="red")+
  geom_vline(data = ice_off,mapping=aes(xintercept = Date), linetype = "dashed", color="blue")+
  geom_line(daily,mapping=aes(x=Date,y=CO2_2_Avg,color="Daily mean"),size=1)+
  geom_point(daily,mapping=aes(x=Date,y=CO2_2_Avg,color="Daily mean"),size=2)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylim(0,2500)+
  theme_classic(base_size = 10)

## Load and plot FDOM data as well (from 1.6 m)
# Download RAW data from GitHub
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data/Catwalk.csv",paste0(getwd(), "/Data/Catwalk_2020.csv"))

fdom <- read.csv("./Data/Catwalk_2020.csv",skip=1) #%>% 
fdom <- fdom %>% 
  mutate(TIMESTAMP = as.POSIXct(strptime(fdom$TIMESTAMP,"%Y-%m-%d %H:%M:%S"))) %>% 
  filter(TIMESTAMP>"2020-12-20 00:00:00")

fdom_2 <- fdom %>% 
  mutate(fDOM_RFU_1 = ifelse(fDOM_RFU_1 == "NAN", NA, fDOM_RFU_1))

fdom_2 <- fdom_2[!is.na(fdom_2$fDOM_RFU_1),]
fdom_2$fDOM_RFU_1 <- as.numeric(fdom_2$fDOM_RFU_1)

ggplot()+
  geom_point(fdom_2,mapping=aes(x=TIMESTAMP,y=fDOM_RFU_1))+
  geom_vline(data = ice_on,mapping=aes(xintercept = Date), linetype = "dashed", color="red")+
  geom_vline(data = ice_off,mapping=aes(xintercept = Date), linetype = "dashed", color="blue")+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylim(5,7)+
  theme_classic(base_size = 10)

