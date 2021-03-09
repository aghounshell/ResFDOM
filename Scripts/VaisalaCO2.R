### Script to explore dissolved CO2 sensor deployed by Mark Johnson ###
# Following script from ABP: https://github.com/CareyLabVT/ManualDownloadsSCCData/blob/master/Scripts/Co2%20sensor.R
# 09 Mar 2021, A Hounshell

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

#create the graph
ggplot(out.file, aes(x=TIMESTAMP, y=CO2_2_Avg, col=CO2_2_Avg)) +
  geom_point()

ggplot(out.file,aes(x=TIMESTAMP,y=CO2_2_Avg))+
  geom_point(size=0.1)+
  geom_line()+
  xlim(as.POSIXct("2020-11-01"),as.POSIXct("2021-03-01"))+
  theme_classic(base_size = 10)

# Create graph to look at ice on/off
ggplot(out.file,aes(x=TIMESTAMP,y=CO2_2_Avg))+
  geom_point(size=0.1)+
  geom_line()+
  xlim(as.POSIXct("2020-12-27"),as.POSIXct("2021-03-01"))+
  ylim(0,2500)+
  theme_classic(base_size = 10)

ggplot(out.file,aes(x=TIMESTAMP,y=CO2_2_Avg))+
  geom_point(size=0.1)+
  geom_line()+
  xlim(as.POSIXct("2021-01-27"),as.POSIXct("2021-02-10"))+
  ylim(0,2500)+
  theme_classic(base_size = 10)
