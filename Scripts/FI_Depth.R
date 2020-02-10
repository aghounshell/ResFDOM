### Script to start exploring Reservoir FDOM data from Summer 2019
### A Hounshell, 19 Nov 19

## Rfile saved as: FI_Depth

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in data: EEMs 'results' files
data <- read_csv("C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20200210_ResultsFiles_ResEEMs2019.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%m/%d/%Y", tz = "EST"))
fcr <- data %>% filter(Reservoir == "FCR")
fcr_epi <- fcr %>% filter(Depth == 0.1)
fcr_meta <- fcr %>% filter(Depth == 5.0)
fcr_hypo <- fcr %>% filter(Depth == 9.0)

# Just start plotting?
ggplot(fcr,aes(x=Date,y=HIX,group=Depth,color=Depth))+
  geom_point()+
  geom_line()+
  theme_classic()

ggplot(fcr,aes(x=Date,y=BIX,group=Depth,color=Depth))+
  geom_point()+
  geom_line()+
  theme_classic()

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=HIX),color="#F5793A",size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=HIX),color="#F5793A",size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=HIX),color="#A95AA1",size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=HIX),color="#A95AA1",size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=HIX),color="#85C0F9",size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=HIX),color="#85C0F9",size=1)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  theme_classic(base_size=15)

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=BIX),color="blue")+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=BIX),color="blue")+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=BIX),color="orange")+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=BIX),color="orange")+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=BIX),color="green")+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=BIX),color="green")+
  theme_classic()

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=T),color="blue")+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=T),color="blue")+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=T),color="orange")+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=T),color="orange")+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=T),color="green")+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=T),color="green")+
  theme_classic()

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=A),color="blue")+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=A),color="blue")+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=A),color="orange")+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=A),color="orange")+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=A),color="green")+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=A),color="green")+
  theme_classic()

date_char <- fcr
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
