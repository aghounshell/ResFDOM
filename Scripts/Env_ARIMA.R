### Script to look at environmental variables and conduct ARIMA modeling between 
### epi and hypo DOC concentrations and various limno. and met. parameters

### 18 Aug 2022, A. Hounshell

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,zoo,scales,plyr,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

###############################################################################

## Load in a format DOC concentration data from 2017-2021
chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DOC_mgL <= 15) %>% 
  drop_na(DOC_mgL)

###############################################################################

## Load in thermocline data - obtained from CTD and YSI casts
## Calculated using LakeAnalyzer
thermo <- read.csv("./Data/rev_FCR_results_LA.csv") %>% 
  select(datetime,thermo.depth) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = datetime)

## Estimate the location of the epi and hypo using thermocline depth
thermo <- thermo %>% 
  mutate(thermo.depth = round(thermo.depth,digits=1)) %>% 
  mutate(epi_bottomg_depth_m = ifelse(thermo.depth > 9.0, 9.0,
                                      ifelse(thermo.depth > 7.0, 8.0,
                                             ifelse(thermo.depth > 6.0, 6.2,
                                                    ifelse(thermo.depth > 4.4, 5.0,
                                                           ifelse(thermo.depth > 3.0, 3.8,
                                                                  ifelse(thermo.depth > 1.6, 1.6,
                                                                         ifelse(thermo.depth > 0.0, 0.1, NA)))))))) %>% 
  mutate(hypo_top_depth_m = ifelse(thermo.depth <= 0.0, 0.1,
                                   ifelse(thermo.depth <= 1.6, 1.6,
                                          ifelse(thermo.depth <= 3.0, 3.8,
                                                 ifelse(thermo.depth <= 4.4, 5.0,
                                                        ifelse(thermo.depth <= 6.0, 6.2,
                                                               ifelse(thermo.depth <= 7.0, 8.0,
                                                                      ifelse(thermo.depth <= 9.0, 9.0, NA))))))))

###############################################################################

## Load in CTD + YSI data - temp, Sal, DO
## From merged spreadsheet in: LakeAnalyzer_thermo.R
casts <- read.csv("./Data/merged_YSI_CTD.csv") %>% 
  filter(depth >= 0) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = time)

## Calculate Epi and Hypo V.W. parameters using thermocline data
## Following Eco_DOC.R

# Create vector of different volumes for each depth: based on L&O-L 2020 paper
vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0), "Vol_m3" = c(138486.51,89053.28,59619.35,40197.90,13943.82,14038.52,1954.71))

### Create a dataframe for cast parameters at each sampling depth
depths <- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

#Initialize an empty matrix with the correct number of rows and columns 
temp <- matrix(data=NA, ncol=ncol(casts), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final <- matrix(data=NA, ncol=1, nrow=0)
dates <- unique(casts$DateTime)

#create a function to chose the matching depth closest to our focal depths
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j = dates[i]
  q <- subset(casts, casts$DateTime == j)
  
  layer1 <- q[q[, "depth"] == closest(q$depth,0.1),][1,]
  layer2 <- q[q[, "depth"] == closest(q$depth,1.6),][1,]
  layer3 <- q[q[, "depth"] == closest(q$depth,3.8),][1,]
  layer4 <- q[q[, "depth"] == closest(q$depth,5.0),][1,]
  layer5 <- q[q[, "depth"] == closest(q$depth,6.2),][1,]
  layer6 <- q[q[, "depth"] == closest(q$depth,8.0),][1,]
  layer7 <- q[q[, "depth"] == closest(q$depth,9.0),][1,]
  
  temp <- rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(casts))+1)] <- depths
  colnames(temp)[((ncol(casts))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
casts_depths <- as.data.frame(super_final) %>%
  select(-c(1,depth)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

# Separate by Env parameter of interest and pivot longer
# Temperature
temp_c <- casts_depths %>% 
  select(DateTime,new_depth,Temp_C) %>% 
  mutate(Temp_C = as.numeric(Temp_C)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Temp_C, values_fil = NA, values_fn = mean, names_prefix = "Temp_")

temp_c <- left_join(temp_c, thermo, by="DateTime")

temp_c <- temp_c %>% 
  mutate(epi_temp = ifelse(is.na(epi_bottomg_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Temp_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Temp_0.1*vol_depths$Vol_m3[1]+Temp_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_temp = ifelse(is.na(hypo_top_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Temp_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# temp_c <- temp_c[!(temp_c$DateTime = as.POSIXct("2021-08-20") | temp_c$DateTime == as.POSIXct("2020-07-08") | temp_c$DateTime == as.POSIXct("2019-05-30") | temp_c$DateTime == as.POSIXct("2019-04-29") | temp_c$DateTime == as.POSIXct("2017-09-17")),]

temp_c <- temp_c[-c(267,348,353,418,484),]

## Plot data by epi and hypo
temp_c %>%  
  drop_na(epi_temp,hypo_temp) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Temp~(C^o)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

## Dissolved oxygen
do_mgL <- casts_depths %>% 
  select(DateTime,new_depth,DO_mgL) %>% 
  mutate(DO_mgL = as.numeric(DO_mgL)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = DO_mgL, values_fil = NA, values_fn = mean, names_prefix = "DO_")

do_mgL <- left_join(do_mgL, thermo, by="DateTime")

do_mgL <- do_mgL %>% 
  mutate(epi_DO = ifelse(is.na(epi_bottomg_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, DO_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (DO_0.1*vol_depths$Vol_m3[1]+DO_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DO = ifelse(is.na(hypo_top_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, DO_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2021-08-20;479
# 2020-07-08; 413

do_mgL <- do_mgL[-c(413,479),]

do_mgL %>%  
  drop_na(epi_DO,hypo_DO) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~D.O.~(mg~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

## Conductivity
cond_uScm <- casts_depths %>% 
  select(DateTime,new_depth,Cond_uScm) %>% 
  mutate(Cond_uScm = as.numeric(Cond_uScm)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Cond_uScm, values_fil = NA, values_fn = mean, names_prefix = "Cond_")

cond_uScm <- left_join(cond_uScm, thermo, by="DateTime")

cond_uScm <- cond_uScm %>% 
  mutate(epi_Cond = ifelse(is.na(epi_bottomg_depth_m), ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                         ifelse(epi_bottomg_depth_m == 0.1, Cond_0.1,
                                ifelse(epi_bottomg_depth_m == 1.6, (Cond_0.1*vol_depths$Vol_m3[1]+Cond_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                       ifelse(epi_bottomg_depth_m == 3.8, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                              ifelse(epi_bottomg_depth_m == 5, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                     ifelse(epi_bottomg_depth_m == 6.2, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                            ifelse(epi_bottomg_depth_m == 8, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Cond = ifelse(is.na(hypo_top_depth_m), ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(hypo_top_depth_m == 1.6, ((Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                 ifelse(hypo_top_depth_m == 3.8, ((Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                        ifelse(hypo_top_depth_m == 5, ((Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                               ifelse(hypo_top_depth_m == 6.2, ((Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                      ifelse(hypo_top_depth_m == 8, ((Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                             ifelse(hypo_top_depth_m == 9, Cond_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2017-03-12; 208
cond_uScm <- cond_uScm[-c(208),]

cond_uScm %>%  
  drop_na(epi_Cond,hypo_Cond) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Cond,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Cond,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Cond,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Cond,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Cond.~(~mu*s~cm^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Turbidity
turb_ntu <- casts_depths %>% 
  select(DateTime,new_depth,Turb_NTU) %>% 
  mutate(Turb_NTU = as.numeric(Turb_NTU)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Turb_NTU, values_fil = NA, values_fn = mean, names_prefix = "Turb_")

turb_ntu <- left_join(turb_ntu, thermo, by="DateTime")

turb_ntu <- turb_ntu %>% 
  mutate(epi_Turb = ifelse(is.na(epi_bottomg_depth_m), ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Turb_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Turb_0.1*vol_depths$Vol_m3[1]+Turb_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Turb = ifelse(is.na(hypo_top_depth_m), ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Turb_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2019-07-03; 208
turb_ntu <- turb_ntu[-c(270),]

turb_ntu %>%  
  drop_na(epi_Turb,hypo_Turb) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Turb,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Turb,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Turb,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Turb,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylim(0,50)+
  ylab(expression(V.W.~Turb.~(~NTU)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

###############################################################################

## Load in Flora data - Chla and community analysis

#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/5/a24f4dbe9f0d856688f257547069d0a3" 
#infile1 <- paste0(getwd(),"/Data/FluoroProbe.csv")
#download.file(inUrl1,infile1,method="curl")
flora <- read.csv("./Data/FluoroProbe_2014_2021.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2015-01-01")) %>% 
  filter(Site == 50)

###############################################################################

## Load in inflow

###############################################################################

## Load in met data - shortwave radiation and rainfall

###############################################################################