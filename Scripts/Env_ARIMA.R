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

## Check average thermo depth from April-Oct from 2017-2021
thermo %>% 
  mutate(month = month(DateTime)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & month %in% c(4,5,6,7,8,9)) %>% 
  summarise_all(median,na.rm=TRUE)

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

## Calculate VW epi and hypo DOC concentrations

# Create vector of different volumes for each depth: based on L&O-L 2020 paper
vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0), "Vol_m3" = c(138486.51,89053.28,59619.35,40197.90,13943.82,14038.52,1954.71))

# Select DOC concentrations and pivot longer
doc_mgL <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_")

doc_mgL <- left_join(doc_mgL, thermo, by="DateTime")

doc_mgL <- doc_mgL %>% 
  mutate(epi_DOC = ifelse(is.na(epi_bottomg_depth_m), ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(epi_bottomg_depth_m == 0.1, DOC_0.1,
                                 ifelse(epi_bottomg_depth_m == 1.6, (DOC_0.1*vol_depths$Vol_m3[1]+DOC_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                        ifelse(epi_bottomg_depth_m == 3.8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                               ifelse(epi_bottomg_depth_m == 5, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                      ifelse(epi_bottomg_depth_m == 6.2, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                             ifelse(epi_bottomg_depth_m == 8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DOC = ifelse(is.na(hypo_top_depth_m), ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(hypo_top_depth_m == 1.6, ((DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                  ifelse(hypo_top_depth_m == 3.8, ((DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                         ifelse(hypo_top_depth_m == 5, ((DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                ifelse(hypo_top_depth_m == 6.2, ((DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                       ifelse(hypo_top_depth_m == 8, ((DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                              ifelse(hypo_top_depth_m == 9, DOC_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))
  
# Format for ARIMA modeling
final_doc_mgL <- doc_mgL %>% 
  select(DateTime,epi_DOC,hypo_DOC) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_DOC_mgL") %>% 
  mutate(Depth = ifelse(Depth == "epi_DOC", "Epi",
                        ifelse(Depth == "hypo_DOC", "Hypo", NA)))

###############################################################################

## Load in CTD + YSI data - temp, Sal, DO
## From merged spreadsheet in: LakeAnalyzer_thermo.R
casts <- read.csv("./Data/merged_YSI_CTD.csv") %>% 
  filter(depth >= 0) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = time)

## Calculate Epi and Hypo V.W. parameters using thermocline data
## Following Eco_DOC.R

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

# Format for ARIMA modeling
final_temp_c <- temp_c %>% 
  select(DateTime,epi_temp,hypo_temp) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Temp_C") %>% 
  mutate(Depth = ifelse(Depth == "epi_temp", "Epi",
                        ifelse(Depth == "hypo_temp", "Hypo", NA)))
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

# Format for ARIMA modeling
final_do_mgL <- do_mgL %>% 
  select(DateTime,epi_DO,hypo_DO) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_DO_mgL") %>% 
  mutate(Depth = ifelse(Depth == "epi_DO", "Epi",
                        ifelse(Depth == "hypo_DO", "Hypo", NA)))

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

# Format for ARIMA modeling
final_cond_uScm <- cond_uScm %>% 
  select(DateTime,epi_Cond,hypo_Cond) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Cond_uScm") %>% 
  mutate(Depth = ifelse(Depth == "epi_Cond", "Epi",
                        ifelse(Depth == "hypo_Cond", "Hypo", NA)))

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

# Format for ARIMA modeling
final_turb_ntu <- turb_ntu %>% 
  select(DateTime,epi_Turb,hypo_Turb) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Turb_NTU") %>% 
  mutate(Depth = ifelse(Depth == "epi_Turb", "Epi",
                        ifelse(Depth == "hypo_Turb", "Hypo", NA)))

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

## Create dataframe of Flora data that is closest to pre-defined sampling depths; then calculate epi and hypo VW Chla concentrations
#Initialize an empty matrix with the correct number of rows and columns 
temp <- matrix(data=NA, ncol=ncol(flora), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final <- matrix(data=NA, ncol=1, nrow=0)
dates <- unique(flora$DateTime)

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j = dates[i]
  q <- subset(flora, flora$DateTime == j)
  
  layer1 <- q[q[, "Depth_m"] == closest(q$Depth_m,0.1),][1,]
  layer2 <- q[q[, "Depth_m"] == closest(q$Depth_m,1.6),][1,]
  layer3 <- q[q[, "Depth_m"] == closest(q$Depth_m,3.8),][1,]
  layer4 <- q[q[, "Depth_m"] == closest(q$Depth_m,5.0),][1,]
  layer5 <- q[q[, "Depth_m"] == closest(q$Depth_m,6.2),][1,]
  layer6 <- q[q[, "Depth_m"] == closest(q$Depth_m,8.0),][1,]
  layer7 <- q[q[, "Depth_m"] == closest(q$Depth_m,9.0),][1,]
  
  temp <- rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(flora))+1)] <- depths
  colnames(temp)[((ncol(flora))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
flora_depths <- as.data.frame(super_final) %>%
  select(-c(1,Depth_m)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Calculate VW Epi and Hypo concentrations
chla_ugL <- flora_depths %>% 
  select(DateTime,new_depth,TotalConc_ugL) %>% 
  mutate(TotalConc_ugL = as.numeric(TotalConc_ugL)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = TotalConc_ugL, values_fil = NA, values_fn = mean, names_prefix = "Chla_")

chla_ugL <- left_join(chla_ugL, thermo, by="DateTime")

chla_ugL <- chla_ugL %>% 
  mutate(epi_Chla = ifelse(is.na(epi_bottomg_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Chla_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Chla_0.1*vol_depths$Vol_m3[1]+Chla_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Chla = ifelse(is.na(hypo_top_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Chla_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Plot
chla_ugL %>%  
  drop_na(epi_Chla,hypo_Chla) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Chla~(~mu*g~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_chla_ugL <- chla_ugL %>% 
  select(DateTime,epi_Chla,hypo_Chla) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Chla_ugL") %>% 
  mutate(Depth = ifelse(Depth == "epi_Chla", "Epi",
                        ifelse(Depth == "hypo_Chla", "Hypo", NA)))

###############################################################################

## Load in inflow
inflow <- read.csv("./Data/Inflow_2013_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:VT_Temp_C)

# Find relationship between WVWA and VT flow
# Use to fill in any missing WVWA values with VT flow (where possible)
flow_lm <- lm(WVWA_Flow_cms ~ VT_Flow_cms, data = inflow)

inflow <- inflow %>% 
  mutate(WVWA_Flow_cms = ifelse(is.na(WVWA_Flow_cms), 0.8917528*VT_Flow_cms-0.0009302, WVWA_Flow_cms)) %>% 
  mutate(flow_diff = abs(VT_Flow_cms - WVWA_Flow_cms))

# Average inflow by day
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_at(vars("WVWA_Flow_cms"),funs(mean(.,na.rm=TRUE),sd)) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

# Plot daily inflow for the study period
inflow_daily %>% 
  na.omit(mean) %>% 
  ggplot(mapping=aes(x=DateTime,y=mean))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  geom_point(size=2)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Inflow~(~m^3~s^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_inflow_m3s <- inflow_daily %>% 
  select(DateTime,mean) %>% 
  dplyr::rename(Inflow_m3s = mean)
  
###############################################################################

## Load in met data - shortwave radiation and rainfall
# Downloaded on 23 Aug 2022
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb"
#infile1 <- paste0(getwd(),"/Data/FCR_Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

## Load in Met data
met <- read.csv("./Data/FCR_Met_final_2015_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01 00:00:00"))

## Average to daily
met_daily <- met %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  dplyr::summarise(rain_tot_mm = sum(Rain_Total_mm, na.rm = TRUE),
            ShortwaveRadiationUp_Average_W_m2 = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE))

## Then plot total rainfall and shorwave radiation
met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=rain_tot_mm))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Total~Rainfall~(mm)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=ShortwaveRadiationUp_Average_W_m2))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(S.W.~Radiation~(W~m^2)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

###############################################################################

## Load in dissolved GHG data - as a potential explanatory variable for DOC concentrations...
# Downloaded on 23 Aug 2022
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/6/38d72673295864956cccd6bbba99a1a3"
#infile1 <- paste0(getwd(),"/Data/final_GHG_2015-2021.csv")
#download.file(inUrl1,infile1,method="curl")

## Load in Dissolved GHG data
ghg <- read.csv("./Data/final_GHG_2015-2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01 00:00:00") & Reservoir == "FCR" & Site == 50) %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST"))

ghg <- ghg %>% 
  group_by(DateTime,Depth_m) %>% 
  dplyr::summarise(mean_ch4_umolL = mean(ch4_umolL, na.rm=TRUE),
            mean_co2_umolL = mean(co2_umolL, na.rm=TRUE))

# Pivot_longer for each GHG and each sampling time point
ch4_umolL <- ghg %>% 
  select(DateTime,Depth_m,mean_ch4_umolL) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = mean_ch4_umolL, values_fil = NA, values_fn = mean, names_prefix = "ch4_") %>% 
  select(ch4_0.1,ch4_1.6,ch4_3.8,ch4_5,ch4_6.2,ch4_8,ch4_9)

ch4_umolL <- left_join(ch4_umolL, thermo, by="DateTime")

ch4_umolL <- ch4_umolL %>% 
  mutate(epi_ch4 = ifelse(is.na(epi_bottomg_depth_m), ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, ch4_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (ch4_0.1*vol_depths$Vol_m3[1]+ch4_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_ch4 = ifelse(is.na(hypo_top_depth_m), ((ch4_0.1*vol_depths$Vol_m3[1])+(ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((ch4_1.6*vol_depths$Vol_m3[2])+(ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((ch4_3.8*vol_depths$Vol_m3[3])+(ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((ch4_5*vol_depths$Vol_m3[4])+(ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((ch4_6.2*vol_depths$Vol_m3[5])+(ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((ch4_8*vol_depths$Vol_m3[6])+(ch4_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, ch4_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## CO2
co2_umolL <- ghg %>% 
  select(DateTime,Depth_m,mean_co2_umolL) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = mean_co2_umolL, values_fil = NA, values_fn = mean, names_prefix = "co2_") %>% 
  select(co2_0.1,co2_1.6,co2_3.8,co2_5,co2_6.2,co2_8,co2_9)

co2_umolL <- left_join(co2_umolL, thermo, by="DateTime")

co2_umolL <- co2_umolL %>% 
  mutate(epi_co2 = ifelse(is.na(epi_bottomg_depth_m), ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(epi_bottomg_depth_m == 0.1, co2_0.1,
                                 ifelse(epi_bottomg_depth_m == 1.6, (co2_0.1*vol_depths$Vol_m3[1]+co2_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                        ifelse(epi_bottomg_depth_m == 3.8, ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                               ifelse(epi_bottomg_depth_m == 5, ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                      ifelse(epi_bottomg_depth_m == 6.2, ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                             ifelse(epi_bottomg_depth_m == 8, ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_co2 = ifelse(is.na(hypo_top_depth_m), ((co2_0.1*vol_depths$Vol_m3[1])+(co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(hypo_top_depth_m == 1.6, ((co2_1.6*vol_depths$Vol_m3[2])+(co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                  ifelse(hypo_top_depth_m == 3.8, ((co2_3.8*vol_depths$Vol_m3[3])+(co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                         ifelse(hypo_top_depth_m == 5, ((co2_5*vol_depths$Vol_m3[4])+(co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                ifelse(hypo_top_depth_m == 6.2, ((co2_6.2*vol_depths$Vol_m3[5])+(co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                       ifelse(hypo_top_depth_m == 8, ((co2_8*vol_depths$Vol_m3[6])+(co2_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                              ifelse(hypo_top_depth_m == 9, co2_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Plot
ch4_umolL %>%  
  drop_na(epi_ch4,hypo_ch4) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_ch4,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_ch4,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_ch4,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_ch4,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~CH[4]~(~mu*mol~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_umolL %>%  
  drop_na(epi_co2,hypo_co2) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_co2,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_co2,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_co2,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_co2,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~CO[2]~(~mu*mol~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_ch4_umolL <- ch4_umolL %>% 
  select(DateTime,epi_ch4,hypo_ch4) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_CH4_umolL") %>% 
  mutate(Depth = ifelse(Depth == "epi_ch4", "Epi",
                        ifelse(Depth == "hypo_ch4", "Hypo", NA)))

final_co2_umolL <- co2_umolL %>% 
  select(DateTime,epi_co2,hypo_co2) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_co2_umolL") %>% 
  mutate(Depth = ifelse(Depth == "epi_co2", "Epi",
                        ifelse(Depth == "hypo_co2", "Hypo", NA)))

###############################################################################

## Organize data for ARIMA modeling - following EddyFlux, 3_Rev_EnvAnalysis.R
# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, SW Radiation, CO2, and CH4 for Epi and Hypo
arima_epi <- join_all(list(final_doc_mgL,final_temp_c,final_do_mgL,final_chla_ugL,final_ch4_umolL,final_co2_umolL),by=c("DateTime","Depth"),type="left") 

arima_epi <- join_all(list(arima_epi,met_daily,final_inflow_m3s),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Epi") %>% 
  drop_na(VW_DOC_mgL)

# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, SW Radiation, CO2, and CH4 for Epi and Hypo
arima_hypo <- join_all(list(final_doc_mgL,final_temp_c,final_do_mgL,final_chla_ugL,final_ch4_umolL,final_co2_umolL),by=c("DateTime","Depth"),type="left")

arima_hypo <- join_all(list(arima_hypo,met_daily,final_inflow_m3s),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Hypo") %>% 
  drop_na(VW_DOC_mgL)

###############################################################################

## Check correlations among environmental variables - what needs to be removed?
# Epi
epi_cor = as.data.frame(cor(arima_epi[,3:11],use = "complete.obs"),method=c("pearson"))
chart.Correlation(arima_epi[,3:11],histogram = TRUE,method=c("pearson"))
# No correlations!

# Hypo
hypo_cor = as.data.frame(cor(arima_hypo[,3:11],use = "complete.obs"),method=c("pearson"))
chart.Correlation(arima_hypo[,3:11],histogram = TRUE,method=c("pearson"))
# DO and CO2 correlated at -0.73...but going to keep it in and use a cut-off of 0.75?

###############################################################################

## Check for skewness following MEL script!
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# Epi
for (i in 3:11){
  print(colnames(arima_epi)[i])
  var <- arima_epi[,i]
  hist(as.matrix(var), main = colnames(arima_epi)[i])
  print(skewness(arima_epi[,i], na.rm = TRUE))
  print(skewness(log(arima_epi[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(arima_epi[,i]), na.rm = TRUE))
  print(skewness(arima_epi[,i]^2),na.rm=TRUE)
  var <- log(arima_epi[,i])
  hist(as.matrix(var), main = c("Log",colnames(arima_epi)[i]))
  var <- Math.cbrt(arima_epi[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(arima_epi)[i]))
  var <- (arima_epi[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(arima_epi)[i]))
}
# Nothing: Temp, DO, CH4, CO2, SW Radiation, Rain
# Log: DOC, Chla, Inflow

# Transform and scale data
arima_epi_scale <- arima_epi %>% 
  mutate(VW_DOC_mgL = log(VW_DOC_mgL),
         VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s))

arima_epi_scale[,3:11] <- scale(arima_epi_scale[,3:11])

# Hypo
for (i in 3:11){
  print(colnames(arima_hypo)[i])
  var <- arima_hypo[,i]
  hist(as.matrix(var), main = colnames(arima_hypo)[i])
  print(skewness(arima_hypo[,i], na.rm = TRUE))
  print(skewness(log(arima_hypo[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(arima_hypo[,i]), na.rm = TRUE))
  print(skewness(arima_hypo[,i]^2),na.rm=TRUE)
  var <- log(arima_hypo[,i])
  hist(as.matrix(var), main = c("Log",colnames(arima_hypo)[i]))
  var <- Math.cbrt(arima_hypo[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(arima_hypo)[i]))
  var <- (arima_hypo[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(arima_hypo)[i]))
}
# Nothing: DOC, Temp, DO, CH4, CO2, Rain, SW Radiation
# Log: Chla, Inflow

# Transform and scale data
arima_hypo_scale <- arima_hypo %>% 
  mutate(VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s))

arima_hypo_scale[,3:11] <- scale(arima_hypo_scale[,3:11])

###############################################################################

## ARIMA modeling
# Following MEL code : )
# THINGS TO CHANGE: 'cols' (change to the environmental variables); best fit!

###############################################################################

# Epi
colnames(arima_epi_scale)

cols <- c(4:11) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_epi_scale[,3] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "epi"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "epi"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for epi",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "epi"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_epi_scale)[combn(cols,4)[,22]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,22] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_epi_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_epi_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_epi_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_epi_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################

# Epi
colnames(arima_hypo_scale)

cols <- c(4:7,9:11) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,3] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "hypo"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "hypo"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for hypo",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "hypo"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_hypo_scale)[combn(cols,5)[,3]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,5)[,3] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_hypo_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################