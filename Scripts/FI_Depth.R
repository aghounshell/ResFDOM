### Script to start exploring Reservoir FDOM data from Summer 2019
### A Hounshell, 19 Nov 19

## Rfile saved as: FI_Depth

## Updated with more data
## A Hounshell, 08 Aug 2020

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in data: EEMs 'results' files
data <- read_csv("C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/20200928_ResultsFiles_ResEEMs2019.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%m/%d/%Y", tz = "EST"))
data$HIX <- as.numeric(data$HIX)
data$BIX <- as.numeric(data$BIX)
fcr <- data %>% filter(Reservoir == "FCR")
fcr_epi <- fcr %>% filter(Station == 50 & Depth == 0.1) %>% group_by(Date) %>% 
  summarise_all(funs(mean))
fcr_epi_sd <- fcr %>% filter(Station == 50 & Depth == 0.1) %>% group_by(Date) %>% 
  summarise_all(funs(sd))

# Plot
ggplot(fcr_epi,mapping=aes(x=Date,y=HIX))+
  geom_point()+
  geom_line()+
  geom_errorbar(fcr_epi_sd,mapping=aes(ymin=fcr_epi$HIX-HIX,ymax=fcr_epi$HIX+HIX))+
  theme_classic()

fcr_meta <- fcr %>% filter(Depth == 5.0) %>% group_by(Date) %>% summarise_all(funs(mean))
fcr_meta_sd <- fcr %>% filter(Depth == 5.0) %>% group_by(Date) %>% summarise_all(funs(sd))

ggplot(fcr_meta,mapping=aes(x=Date,y=HIX))+
  geom_point()+
  geom_line()+
  geom_errorbar(fcr_meta_sd,mapping=aes(ymin=fcr_meta$HIX-HIX,ymax=fcr_meta$HIX+HIX))+
  theme_classic()

fcr_hypo <- fcr %>% filter(Depth == 9.0) %>% group_by(Date) %>% summarise_all(funs(mean))
fcr_hypo_sd <- fcr %>% filter(Depth == 9.0) %>% group_by(Date) %>% summarise_all(funs(sd))

ggplot(fcr_hypo,mapping=aes(x=Date,y=HIX))+
  geom_point()+
  geom_line()+
  geom_errorbar(fcr_hypo_sd,mapping=aes(ymin=fcr_hypo$HIX-HIX,ymax=fcr_hypo$HIX+HIX))+
  theme_classic()

# Plot all
ggplot()+
  geom_point(fcr_epi,mapping=aes(x=Date,y=HIX,color="Epi"))+
  geom_line(fcr_epi,mapping=aes(x=Date,y=HIX,color="Epi"))+
  geom_errorbar(fcr_epi_sd,mapping=aes(x=Date,y=fcr_epi$HIX,ymin=fcr_epi$HIX-HIX,ymax=fcr_epi$HIX+HIX,color="Epi"))+
  geom_point(fcr_meta,mapping=aes(x=Date,y=HIX,color="Meta"))+
  geom_line(fcr_meta,mapping=aes(x=Date,y=HIX,color="Meta"))+
  geom_errorbar(fcr_meta_sd,mapping=aes(x=Date,y=fcr_meta$HIX,ymin=fcr_meta$HIX-HIX,ymax=fcr_meta$HIX+HIX,color="Meta"))+
  geom_point(fcr_hypo,mapping=aes(x=Date,y=HIX,color="Hypo"))+
  geom_line(fcr_hypo,mapping=aes(x=Date,y=HIX,color="Hypo"))+
  geom_errorbar(fcr_hypo_sd,mapping=aes(x=Date,y=fcr_hypo$HIX,ymin=fcr_hypo$HIX-HIX,ymax=fcr_hypo$HIX+HIX,color="Hypo"))+
  theme_classic()

ggplot()+
  geom_point(fcr_epi,mapping=aes(x=Date,y=BIX,color="Epi"))+
  geom_line(fcr_epi,mapping=aes(x=Date,y=BIX,color="Epi"))+
  geom_errorbar(fcr_epi_sd,mapping=aes(x=Date,y=fcr_epi$BIX,ymin=fcr_epi$BIX-BIX,ymax=fcr_epi$BIX+BIX,color="Epi"))+
  geom_point(fcr_meta,mapping=aes(x=Date,y=BIX,color="Meta"))+
  geom_line(fcr_meta,mapping=aes(x=Date,y=BIX,color="Meta"))+
  geom_errorbar(fcr_meta_sd,mapping=aes(x=Date,y=fcr_meta$BIX,ymin=fcr_meta$BIX-BIX,ymax=fcr_meta$BIX+BIX,color="Meta"))+
  geom_point(fcr_hypo,mapping=aes(x=Date,y=BIX,color="Hypo"))+
  geom_line(fcr_hypo,mapping=aes(x=Date,y=BIX,color="Hypo"))+
  geom_errorbar(fcr_hypo_sd,mapping=aes(x=Date,y=fcr_hypo$BIX,ymin=fcr_hypo$BIX-BIX,ymax=fcr_hypo$BIX+BIX,color="Hypo"))+
  theme_classic()

fcr_inf <- fcr %>% filter(Station == 100)
fcr_wet <- fcr %>% filter(Station == 200)

# Plot
ggplot()+
  geom_point(fcr_inf,mapping=aes(x=Date,y=HIX,color="Inf"))+
  geom_line(fcr_inf,mapping=aes(x=Date,y=HIX,color="Inf"))+
  geom_point(fcr_wet,mapping=aes(x=Date,y=HIX,color="Wet"))+
  geom_line(fcr_wet,mapping=aes(x=Date,y=HIX,color="Wet"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(fcr_inf,mapping=aes(x=Date,y=BIX,color="Inf"))+
  geom_line(fcr_inf,mapping=aes(x=Date,y=BIX,color="Inf"))+
  geom_point(fcr_wet,mapping=aes(x=Date,y=BIX,color="Wet"))+
  geom_line(fcr_wet,mapping=aes(x=Date,y=BIX,color="Wet"))+
  theme_classic(base_size=15)

s50 <- ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=HIX,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=HIX,color='Epi'),size=1)+
  geom_errorbar(fcr_epi_sd,mapping=aes(x=Date,y=fcr_epi$HIX,ymin=fcr_epi$HIX-HIX,ymax=fcr_epi$HIX+HIX,color="Epi"))+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=HIX,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=HIX,color='Meta'),size=1)+
  geom_errorbar(fcr_meta_sd,mapping=aes(x=Date,y=fcr_meta$HIX,ymin=fcr_meta$HIX-HIX,ymax=fcr_meta$HIX+HIX,color="Meta"))+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=HIX,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=HIX,color='Hypo'),size=1)+
  geom_errorbar(fcr_hypo_sd,mapping=aes(x=Date,y=fcr_hypo$HIX,ymin=fcr_hypo$HIX-HIX,ymax=fcr_hypo$HIX+HIX,color="Hypo"))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  geom_hline(yintercept = 6, color="grey")+
  ylim(0,10)+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('Epi','Meta','Hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+
  theme_classic(base_size=15)

inf <- ggplot()+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=HIX,color='Inf'),size=3)+
  geom_line(data=fcr_inf,mapping=aes(x=Date,y=HIX,color='Inf'),size=1)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=HIX,color='Wet'),size=3)+
  geom_line(data=fcr_wet,mapping=aes(x=Date,y=HIX,color='Wet'),size=1)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  geom_hline(yintercept = 6, color="grey")+
  ylim(0,10)+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('Inf','Wet'),values=c("#F0B670","#FE5F55"))+
  theme_classic(base_size=15)

ggarrange(s50,inf,ncol=1,nrow=2)

s50 <- ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=BIX,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=BIX,color='Epi'),size=1)+
  geom_errorbar(fcr_epi_sd,mapping=aes(x=Date,y=fcr_epi$BIX,ymin=fcr_epi$BIX-BIX,ymax=fcr_epi$BIX+BIX,color="Epi"))+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=BIX,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=BIX,color='Meta'),size=1)+
  geom_errorbar(fcr_meta_sd,mapping=aes(x=Date,y=fcr_meta$BIX,ymin=fcr_meta$BIX-BIX,ymax=fcr_meta$BIX+BIX,color="Meta"))+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=BIX,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=BIX,color='Hypo'),size=1)+
  geom_errorbar(fcr_hypo_sd,mapping=aes(x=Date,y=fcr_hypo$BIX,ymin=fcr_hypo$BIX-BIX,ymax=fcr_hypo$BIX+BIX,color="Hypo"))+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,0.9)+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('Epi','Meta','Hypo'),values=c("#7FC6A4","#7EBDC2","#393E41"))+  
  theme_classic(base_size=15)

inf <- ggplot()+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=BIX,color='Inf'),size=3)+
  geom_line(data=fcr_inf,mapping=aes(x=Date,y=BIX,color='Inf'),size=1)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=BIX,color='Wet'),size=3)+
  geom_line(data=fcr_wet,mapping=aes(x=Date,y=BIX,color='Wet'),size=1)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  ylim(0,0.9)+
  xlim(as.POSIXct("2019-04-28"),as.POSIXct("2019-11-09"))+
  scale_color_manual(breaks=c('Inf','Wet'),values=c("#F0B670","#FE5F55"))+  
  theme_classic(base_size=15)

ggarrange(s50,inf,ncol=1,nrow=2)

# Let's look at RC day data...
fcr_surf <- fcr %>% filter(Depth==0.1)
bvr_surf <- data %>% filter(Reservoir == "BVR")

# Sepearate FCR by inflows and reservoir
fcr_inflows <- fcr_surf %>% filter(Station == 99 | Station == 100 | Station == 101 | Station == 102 | Station == 200)
fcr_res <- fcr_surf %>% filter(Station == 20 | Station == 30 | Station == 45 | Station == 50 | Station == 1)

# Plot FCR
ggplot(fcr_surf,mapping=aes(x=Date,y=HIX,color=as.factor(Station)))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(breaks=c('1','20','30','45','50','99','100','101','102','200'),
                     values=c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#E7804B","#DA2C38","#7DAF4B"))+ 
  ylim(0,9)+
  theme_classic(base_size=15)

ggplot(fcr_surf,mapping=aes(x=Date,y=BIX,color=as.factor(Station)))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(breaks=c('1','20','30','45','50','99','100','101','102','200'),
                     values=c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#E7804B","#DA2C38","#7DAF4B"))+ 
  ylim(0,1)+
  theme_classic(base_size=15)

ggplot(fcr_inflows,mapping=aes(x=Date,y=HIX,color=as.factor(Station)))+
  geom_point()+
  geom_line()+
  theme_classic(base_size=15)

ggplot(fcr_inflows,mapping=aes(x=Date,y=BIX,color=as.factor(Station)))+
  geom_point()+
  geom_line()+
  theme_classic(base_size=15)

ggplot(fcr_res,mapping=aes(x=Date,y=HIX,color=as.factor(Station)))+
  geom_point()+
  geom_line()+
  theme_classic(base_size=15)

ggplot(fcr_res,mapping=aes(x=Date,y=BIX,color=as.factor(Station)))+
  geom_point()+
  geom_line()+
  theme_classic(base_size=15)

# Plot BVR
ggplot(bvr_surf,mapping=aes(x=Date,y=HIX,color=as.factor(Station)))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(breaks=c('1','20','30','45','50','99','100','200'),
                     values=c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#7DAF4B"))+ 
  ylim(0,9)+
  theme_classic(base_size=15)

ggplot(bvr_surf,mapping=aes(x=Date,y=BIX,color=as.factor(Station)))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(breaks=c('1','20','30','45','50','99','100','200'),
                     values=c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#7DAF4B"))+ 
  ylim(0,1)+
  theme_classic(base_size=15)


###############################################


ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=`A/T`,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=`A/T`,color='Epi'),size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=`A/T`,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=`A/T`,color='Meta'),size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=`A/T`,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=`A/T`,color='Hypo'),size=1)+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=`A/T`,color='Inf'),size=3)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=`A/T`,color='Wet'),size=3)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  #ylim(0,0.9)+
  scale_color_manual(breaks=c('Epi','Meta','Hypo','Inf','Wet'),values=c("#F5793A","#A95AA1","#85C0F9","#003366","#006400"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=M,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=M,color='Epi'),size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=M,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=M,color='Meta'),size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=M,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=M,color='Hypo'),size=1)+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=M,color='Inf'),size=3)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=M,color='Wet'),size=3)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  #ylim(0,2)+
  scale_color_manual(breaks=c('Epi','Meta','Hypo','Inf','Wet'),values=c("#F5793A","#A95AA1","#85C0F9","#003366","#006400"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=`T/M`,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=`T/M`,color='Epi'),size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=`T/M`,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=`T/M`,color='Meta'),size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=`T/M`,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=`T/M`,color='Hypo'),size=1)+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=`T/M`,color='Inf'),size=3)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=`T/M`,color='Wet'),size=3)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  #ylim(0,0.9)+
  scale_color_manual(breaks=c('Epi','Meta','Hypo','Inf','Wet'),values=c("#F5793A","#A95AA1","#85C0F9","#003366","#006400"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=M/T,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=M/T,color='Epi'),size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=M/T,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=M/T,color='Meta'),size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=M/T,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=M/T,color='Hypo'),size=1)+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=M/T,color='Inf'),size=3)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=M/T,color='Wet'),size=3)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  #ylim(0,0.9)+
  scale_color_manual(breaks=c('Epi','Meta','Hypo','Inf','Wet'),values=c("#F5793A","#A95AA1","#85C0F9","#003366","#006400"))+
  theme_classic(base_size=15)

ggplot()+
  geom_point(data=fcr_epi,mapping=aes(x=Date,y=T,color='Epi'),size=2)+
  geom_line(data=fcr_epi,mapping=aes(x=Date,y=T,color='Epi'),size=1)+
  geom_point(data=fcr_meta,mapping=aes(x=Date,y=T,color='Meta'),size=2)+
  geom_line(data=fcr_meta,mapping=aes(x=Date,y=T,color='Meta'),size=1)+
  geom_point(data=fcr_hypo,mapping=aes(x=Date,y=T,color='Hypo'),size=2)+
  geom_line(data=fcr_hypo,mapping=aes(x=Date,y=T,color='Hypo'),size=1)+
  geom_point(data=fcr_inf,mapping=aes(x=Date,y=T,color='Inf'),size=3)+
  geom_point(data=fcr_wet,mapping=aes(x=Date,y=T,color='Wet'),size=3)+
  geom_vline(xintercept = as.POSIXct("2019-06-03"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-06-17"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-07-08"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-07-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-08-05"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-08-19"), color="black",linetype="dashed")+ # Oxygen off
  geom_vline(xintercept = as.POSIXct("2019-09-02"), color="black")+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2019-11-02"), color="black",linetype="dashed")+ # Turnover
  #ylim(0,0.9)+
  scale_color_manual(breaks=c('Epi','Meta','Hypo','Inf','Wet'),values=c("#F5793A","#A95AA1","#85C0F9","#003366","#006400"))+
  theme_classic(base_size=15)

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
