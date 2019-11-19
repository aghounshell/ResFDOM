### Script to start exploring Reservoir FDOM data from Summer 2019
### A Hounshell, 19 Nov 19

## Rfile saved as: FI_Depth

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in data: EEMs 'results' files
data <- read_csv("C:/Users/ahoun/Dropbox/ResFDOM/ResFDOM/Data/ResultsFiles_ResEEMs2019.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%m/%d/%Y", tz = "EST"))

# Just start plotting?
ggplot(data,aes(x=Date,y=HIX,group=Depth,color=Depth))+
  geom_point()+
  geom_line()+
  theme_classic()

ggplot(data,aes(x=Date,y=BIX,group=Depth,color=Depth))+
  geom_point()+
  geom_line()+
  theme_classic()

date_char <- data
date_char$Date <- as.character(date_char$Date)

ggplot(date_char,aes(x=HIX,y=Depth,group=Date,color=Date))+
  geom_point()+
  geom_line()+
  scale_y_reverse()+
  theme_classic()

data$Depth <- -1 * data$Depth

# Filter for days with Epi, Meta, and Hypo: 27May19, 3Jul19, 15Jul19, 14Aug19
may <- data %>% filter(Date==as.Date("2019-05-27"))
jul_an <- data %>% filter(Date==as.Date("2019-07-03"))
jul_an <- jul_an[-3,]
jul_ox <- data %>% filter(Date==as.Date("2019-07-15"))
jul_ox <- jul_ox[-3,]
aug <- data %>% filter(Date==as.Date("2019-08-14"))
aug <- aug[-4,]

# Plot by day
date <- rbind.data.frame(may,jul_an,jul_ox,aug)
date$Date <- as.character(date$Date)

HIX <- ggplot(date,aes(x=Depth,y=HIX,group=Date,color=Date))+
  geom_line(size=1)+
  geom_point(size=2)+
  coord_flip()+
  scale_color_manual(labels=c("Anoxic May","Anoxic Jul","Oxic Jul","Oxic Aug"),
                     values=c('#F0750F','#F0A150','#85C0F9','#0F2080'))+
  labs(color="")+
  theme_classic(base_size=15)

BIX <- ggplot(date,aes(x=Depth,y=BIX,group=Date,color=Date))+
  geom_line(size=1)+
  geom_point(size=2)+
  coord_flip()+
  scale_color_manual(labels=c("Anoxic May","Anoxic Jul","Oxic Jul","Oxic Aug"),
                     values=c('#F0750F','#F0A150','#85C0F9','#0F2080'))+
  labs(color="")+
  theme_classic(base_size=15)

ggarrange(HIX,BIX,common.legend=TRUE,legend="right")
