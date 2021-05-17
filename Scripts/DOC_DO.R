### Script to explore long-term DOC patterns and loading to FCR ###
# 29 Mar 2021, A Hounshell
# Updated: 5 Apr 2021, A Hounshell
# Following CCC thoughts: constrain hypo to 9 m (assume 8m = outflow concentration)
# Don't worry about 'entrainment' - should be included in the outflow variable

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate,zoo)

# Load data: long-term DOC; weir discharge; temperature/DO casts
# DOC data from EDI - NOTE: May want to update for 2020 data??
# Updated: 19 Apr 2021 with 2020 data!!!
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98" 
#infile1 <- paste0(getwd(),"/Data/chem.csv")
#download.file(inUrl1,infile1,method="curl")

chem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

### Separate and calculate hypo mass (and change in mass/time)
chem_inputs <- chem %>% 
  filter(Site==50 & Depth_m %in% c(8.0,9.0) | Site==100) %>% 
  filter(DateTime >= as.POSIXct("2015-03-31")) %>% 
  rename(time = DateTime,depth = Depth_m) %>% 
  mutate(loc = ifelse(Site==50 & depth==8.0, "50 8m",
                      ifelse(Site==50 & depth==9.0, "50 9m",
                             ifelse(Site==100, "100",NA)))) %>% 
  mutate(DOC_mgL = ifelse(DOC_mgL > 12.0, NA, DOC_mgL)) # Remove outlier at 8 m (12 mg/L!!!)

chem_inputs <- chem_inputs[!is.na(chem_inputs$DOC_mgL),]

# Figure 3 - [DOC] timeseries for 8m, 9m, 100 from 2015-03-31 to 2020-12-02
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29"), xmax = as.POSIXct("2020-09-11"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-09-25"), xmax = as.POSIXct("2020-12-02"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2016-10-07"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2017-10-30"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(chem_inputs,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=1)+
  geom_point(chem_inputs,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=3)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Time")+
  scale_color_manual(breaks=c('50 8m','50 9m','100'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Site 50 8m","Site 50 9m","Inflow"))+
  xlim(as.POSIXct("2015-03-31"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggsave("./Fig_Output/Fig3.jpg",width=10,height=4,units="in",dpi=320)

chem_50_hypo <- chem %>%
  filter(Site==50) %>% 
  filter(Depth_m == 9.0) %>% 
  rename(time = DateTime,depth = Depth_m) %>% 
  select(time,depth,DOC_mgL) %>% 
  group_by(time,depth) %>% 
  summarize_all(funs(mean),na.rm=TRUE)
chem_50_hypo <- chem_50_hypo[!is.na(chem_50_hypo$DOC_mgL),]

# Create dataframe with daily time step from 2015-01-01 to 2019-12-31
chem_50_hypo_full <- as.data.frame(seq(as.POSIXct("2014-01-01",tz="EST"),as.POSIXct("2020-12-31",tz="EST"),by="days"))
chem_50_hypo_full <- chem_50_hypo_full %>% 
  rename(time = `seq(as.POSIXct("2014-01-01", tz = "EST"), as.POSIXct("2020-12-31", tz = "EST"), by = "days")`)
chem_50_hypo_full <- left_join(chem_50_hypo_full,chem_50_hypo,by="time")

# Extrapolate and constrain to 2015-2019 data
chem_50_hypo_full <- chem_50_hypo_full %>% 
  mutate(DOC_mgL = na.fill(na.approx(DOC_mgL,na.rm=FALSE),"extend")) %>% 
  select(-depth)

# Calculate hypo mass: 9m x vol
# Use same volumes as in L&O-L MS
chem_50_hypo_full <- chem_50_hypo_full %>% 
  mutate(hypo_mg = (DOC_mgL*1.95*10^3*1000)) %>% 
  mutate(dMdt_mgs = NA)

# Calculate dM/dt (mg/s)
for(i in 2:length(chem_50_hypo_full$time)){
  chem_50_hypo_full$dMdt_mgs[i] = (chem_50_hypo_full$hypo_mg[i]-chem_50_hypo_full$hypo_mg[i-1])/(24*60*60)
}

# Plot to check - change in DOC in mg/s
ggplot(chem_50_hypo_full,mapping=aes(time,dMdt_mgs))+
  geom_line()

# Start 'final' datasheet for calculations
box_data <- chem_50_hypo_full

### Separate out thermocline concentration for 'outflow' calculation
chem_50_therm <- chem %>% 
  filter(Site==50) %>% 
  filter(Depth_m %in% c(8.0)) %>% 
  rename(time = DateTime,depth = Depth_m) %>% 
  select(time,DOC_mgL) %>% 
  rename(DOC_mgL_therm = DOC_mgL) %>% 
  mutate(DOC_mgL_therm = ifelse(DOC_mgL_therm > 12.0, NA, DOC_mgL_therm)) # Remove outlier at 8 m (12 mg/L!!!)

# Add to box_data data sheet and interpolate
box_data <- left_join(box_data,chem_50_therm,by="time")
box_data <- box_data %>% 
  mutate(DOC_mgL_therm = na.fill(na.approx(DOC_mgL_therm,na.rm=FALSE),"extend"))

# Filter and add data from 100
chem_100 <- chem %>% 
  filter(Site==100) %>% 
  rename(time = DateTime, depth = Depth_m) %>% 
  select(time,DOC_mgL) %>% 
  rename(DOC_mgL_100 = DOC_mgL)

box_data <- left_join(box_data,chem_100,by="time")
box_data <- box_data %>% 
  mutate(DOC_mgL_100 = na.fill(na.approx(DOC_mgL_100,na.rm=FALSE),"extend"))

# Weir discharge/temperature
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/7/f5fa5de4b49bae8373f6e7c1773b026e" 
#infile1 <- paste0(getwd(),"/Data/inflow_for_EDI_2013_10Jan2021.csv")
#download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/inflow_for_EDI_2013_10Jan2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:VT_Temp_C)

# Average inflow by day
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_all(funs(mean))

### Add together inflow and box data
inflow_box <- inflow_daily %>% 
  select(DateTime,WVWA_Flow_cms) %>% 
  rename(time = DateTime)

box_data <- left_join(box_data,inflow_box,by="time")
box_data <- box_data %>% 
  mutate(WVWA_Flow_cms = na.fill(na.approx(WVWA_Flow_cms,na.rm=FALSE),"extend"))

# Plot flow data to check
ggplot(box_data,mapping=aes(x=time,y=WVWA_Flow_cms))+
  geom_line()

### Calculate 'sink/source' term following Gerling et al. 2016; Kruger et al. 2020
box_data <- box_data %>% 
  mutate(j_mgs = dMdt_mgs-(WVWA_Flow_cms*DOC_mgL_100*1000)+(WVWA_Flow_cms*DOC_mgL_therm*1000)) %>% 
  mutate(j_kgd = j_mgs*60*60*24/1000/1000)

# Select for 2015-2019 during the stratified period
box_data_sel <- box_data %>% 
  filter(time>as.POSIXct("2015-03-31")&time<as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(5,6,7,8,9,10))

# Plot - Figure 4: DOC Box model results
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29"), xmax = as.POSIXct("2020-09-11"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-09-25"), xmax = as.POSIXct("2020-12-02"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2016-10-07"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2017-10-30"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_line(box_data,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"),size=1)+
  geom_point(box_data,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"),size=2)+
  geom_line(box_data,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000,color="Inflow_kgd"),size=1)+
  geom_point(box_data,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000,color="Inflow_kgd"),size=2)+
  geom_line(box_data,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000,color="Outflow_kgd"),size=1)+
  geom_point(box_data,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000,color="Outflow_kgd"),size=2)+
  geom_line(box_data,mapping=aes(x=time,y=j_kgd,color="j_kgd"),size=1)+
  geom_point(box_data,mapping=aes(x=time,y=j_kgd,color="j_kgd"),size=2)+
  ylab(expression(paste("(kg d"^-1*")")))+
  xlab("Time")+
  scale_color_manual(breaks=c('Outflow_kgd','dMdt_kgd','Inflow_kgd','j_kgd'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Outflow","dM/dt","Inflow","Source/Sink"))+
  xlim(as.POSIXct("2015-03-31"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggsave("./Fig_Output/Fig4.jpg",width=10,height=4,units="in",dpi=320)

ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29"), xmax = as.POSIXct("2020-09-11"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-09-25"), xmax = as.POSIXct("2020-12-02"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2016-10-07"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2017-10-30"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(box_data_sel,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"))+
  geom_point(box_data_sel,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"))+
  geom_line(box_data_sel,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000,color="Inflow_kgd"))+
  geom_point(box_data_sel,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000,color="Inflow_kgd"))+
  geom_line(box_data_sel,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000,color="Outflow_kgd"))+
  geom_point(box_data_sel,mapping=aes(x=time,y=WVWA_Flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000,color="Outflow_kgd"))+
  geom_line(box_data_sel,mapping=aes(x=time,y=j_kgd,color="j_kgd"))+
  geom_point(box_data_sel,mapping=aes(x=time,y=j_kgd,color="j_kgd"))+
  xlim(as.POSIXct("2015-03-31"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size=15)

ggplot(box_data,mapping=aes(x=time,y=j_kgd))+
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-12-01"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29"), xmax = as.POSIXct("2020-09-11"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-09-25"), xmax = as.POSIXct("2020-12-02"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  geom_line(size=1)+
  geom_point(size=3)+
  geom_hline(yintercept = 0,linetype="dashed")+
  xlab("Date")+
  ylab(expression(paste("DOC Flux (kg d"^-1*")")))+
  xlim(as.POSIXct("2015-03-31"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())
  
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

### Create a dataframe for CTD/YSI parameters at each sampling depth
depths<- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

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
  layer2<- q[q[, "depth"] == closest(q$depth,1.6),][1,]
  layer3<- q[q[, "depth"] == closest(q$depth,3.8),][1,]
  layer4<- q[q[, "depth"] == closest(q$depth,5.0),][1,]
  layer5<- q[q[, "depth"] == closest(q$depth,6.2),][1,]
  layer6<- q[q[, "depth"] == closest(q$depth,8.0),][1,]
  layer7<- q[q[, "depth"] == closest(q$depth,9.0),][1,]
  
  temp<-rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
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

ctd_50_depths$Temp_C <- as.numeric(ctd_50_depths$Temp_C)
ctd_50_depths$depth <- as.numeric(ctd_50_depths$depth)

ctd_50_hypo <- ctd_50_depths %>% 
  filter(depth == 9.0)

ctd_50_hypo$DO_mgL <- as.numeric(ctd_50_hypo$DO_mgL)

# Combine box model data w/ DO data
hypo_do_box <- left_join(ctd_50_hypo,box_data,by="time")
hypo_do_box$DO_mgL <- as.numeric(hypo_do_box$DO_mgL)

hypo_do_box <- hypo_do_box %>% 
  filter(time >= as.POSIXct("2015-01-01") & time <= as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9,10)) %>% 
  mutate(oxy = ifelse(DO_mgL >= 1.0, 'Oxic',
                      ifelse(DO_mgL < 1.0, 'Anoxic',
                             NA))) %>% 
  mutate(oxy = ifelse(time == as.POSIXct("2017-07-06", tz="EST"), 'Oxic', # Time point in 2017 where it *shouldn't* be anoxic
                      ifelse(time == as.POSIXct("2015-06-08", tz="EST"), 'Oxic', # Time point in 2015 where it *shouldn't* be anoxic
                             oxy))) %>% 
  #filter(time != as.POSIXct("2017-07-06",tz="EST")) %>% # Time point in 2017 where it *shouldn't* be anoxic
  mutate(year = ifelse(time >= as.POSIXct("2015-01-01") & time <= as.POSIXct("2015-12-31"), "2015",
                       ifelse(time >= as.POSIXct("2016-01-01") & time <= as.POSIXct("2016-12-31"), "2016",
                              ifelse(time >= as.POSIXct("2017-01-01") & time <= as.POSIXct("2017-12-31"), "2017",
                                     ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                                            ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                                                   ifelse(time >= as.POSIXct("2020-01-01") & time <= as.POSIXct("2020-12-31"), "2020",
                                                   NA)))))))

hypo_do_box <- hypo_do_box[!is.na(hypo_do_box$DO_mgL),]

# Remove and interpolate DO values that are clearly wonky (aka: show anoxic when oxic pre and post)
# 2015-06-08
hypo_do_box$DO_mgL[2] <- NA
# 2017-07-06
hypo_do_box$DO_mgL[82] <- NA
# 2017-08-02
hypo_do_box$DO_mgL[91] <- NA

# Interpolate NA values
hypo_do_box <- hypo_do_box %>% 
  mutate(DO_mgL = na.fill(na.approx(DO_mgL,na.rm=FALSE),"extend"))

# Fig. 5: Plot DOC and DOC source/sink term by oxic vs. anoxic
# [DOC]
doc_full <- ggplot(hypo_do_box,mapping=aes(oxy,DOC_mgL,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

doc <- ggplot(hypo_do_box,mapping=aes(year,DOC_mgL,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Source/sink term
jterm_full <- ggplot(hypo_do_box,mapping=aes(oxy,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Source/sink term (kg d"^-1*")")))+
  xlab("")+
  ylim(-15,10)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

jterm <- ggplot(hypo_do_box,mapping=aes(year,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Source/sink term (kg d"^-1*")")))+
  xlab("Year")+
  ylim(-15,10)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Curious to plot dM/dt
ggplot(hypo_do_box,mapping=aes(oxy,dMdt_mgs*60*60*24/1000/1000,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("dM/dt (kg d"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(doc_full, jterm_full, doc,jterm,ncol=2,nrow=2,common.legend=TRUE)

ggsave("./Fig_Output/Fig5.jpg",width=10,height=8,units="in",dpi=320)

# Plot temp by year (as box plots?)
temp_box <- ggplot(hypo_do_box,mapping=aes(year,Temp_C,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Temp (C"^o*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Plot DO as box plots
do_box <- ggplot(hypo_do_box,mapping=aes(year,DO_mgL,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("DO (mg"^-1*"L"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(temp_box,do_box,ncol=2,nrow=1,common.legend=TRUE)

ggsave("./Fig_OutPut/DO_temp_Year.jpg",width=10,height=5,units="in",dpi=320)

all_med <- hypo_do_box %>% 
  select(-DO_pSat,-Cond_uScm,-Chla_ugL,-Turb_NTU,-pH,-ORP_mV,-PAR_umolm2s) %>% 
  group_by(year,oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE)

year_med <- hypo_do_box %>% 
  select(-DO_pSat,-Cond_uScm,-Chla_ugL,-Turb_NTU,-pH,-ORP_mV,-PAR_umolm2s,-oxy) %>% 
  group_by(year) %>% 
  summarize_all(funs(median),na.rm=TRUE)

oxy_med <- hypo_do_box %>% 
  select(-DO_pSat,-Cond_uScm,-Chla_ugL,-Turb_NTU,-pH,-ORP_mV,-PAR_umolm2s,-year) %>% 
  group_by(oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE)

#Plot Inflow and dM/dt by year
dm_dt <- ggplot(hypo_do_box,mapping=aes(year,dMdt_mgs*60*60*24/1000/1000,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("dM/dt (kg d"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

inflow <- ggplot(hypo_do_box,mapping=aes(year,WVWA_Flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Inflow (kg d"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

outflow <- ggplot(hypo_do_box,mapping=aes(year,WVWA_Flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Outflow (kg d"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

jterm <- ggplot(hypo_do_box,mapping=aes(year,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Source/sink term (kg d"^-1*")")))+
  xlab("Year")+
  ylim(-15,10)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(inflow,outflow,dm_dt,jterm,nrow=2,ncol=2,common.legend = TRUE)

ggsave("./Fig_Output/dmtodt_inflow.jpg",width=10,height=8,units="in",dpi=320)

### Let's start thinking about DOM quality - what metrics to use?
# a254? HIX? BIX? Peak C? Peak T? PARAFAC?
# Load in data
fdom <- read_csv("./Data/20210210_ResultsFiles_ResEEMs2019_RAW.csv") %>% 
  filter(Reservoir == "FCR" & Dilution %in% c(1,2)) %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%m/%d/%Y", tz="EST")))

fdom_hypo <- fdom %>% 
  filter(Depth == 9.0) %>% 
  group_by(Date,Depth) %>% 
  summarize_all(funs(mean))

fdom_hypo_sd <- fdom %>% 
  filter(Depth == 9.0) %>% 
  group_by(Date,Depth) %>% 
  summarize_all(funs(sd))
  
cdom <- read_csv("./Data/20200930_ResultsFiles_Abs2019.csv") %>% 
  filter(Reservoir == "FCR" & Dilution %in% c(1)) %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%m/%d/%Y", tz="EST")))

cdom_hypo <- cdom %>% 
  filter(Depth == 9.0) %>% 
  group_by(Date,Depth) %>% 
  summarize_all(funs(mean))

cdom_hypo_sd <- cdom %>% 
  filter(Depth == 9.0) %>% 
  group_by(Date,Depth) %>% 
  summarize_all(funs(sd))

# Plot various parameters for Hypo
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(cdom_hypo,mapping=aes(x=Date,y=a254))+
  geom_point(cdom_hypo,mapping=aes(x=Date,y=a254))+
  geom_errorbar(cdom_hypo,mapping=aes(x=Date,y=a254,ymin=a254-cdom_hypo_sd$a254,ymax=a254+cdom_hypo_sd$a254))+
  ylim(0,80)+
  theme_classic(base_size=15)

ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(cdom_hypo,mapping=aes(x=Date,y=a350))+
  geom_point(cdom_hypo,mapping=aes(x=Date,y=a350))+
  geom_errorbar(cdom_hypo,mapping=aes(x=Date,y=a350,ymin=a350-cdom_hypo_sd$a350,ymax=a350+cdom_hypo_sd$a350))+
  theme_classic(base_size=15)

ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(cdom_hypo,mapping=aes(x=Date,y=Sr))+
  geom_point(cdom_hypo,mapping=aes(x=Date,y=Sr))+
  geom_errorbar(cdom_hypo,mapping=aes(x=Date,y=Sr,ymin=Sr-cdom_hypo_sd$Sr,ymax=Sr+cdom_hypo_sd$Sr))+
  theme_classic(base_size=15)

fdom_hypo_2 <- fdom_hypo %>% 
  rename(time = Date,PeakT=T,PeakA=A)

cdom_hypo_2 <- cdom_hypo %>% 
  rename(time = Date)

# Combine with DO data
hypo_do_box_2 <- left_join(hypo_do_box,cdom_hypo_2,by="time")
hypo_do_box_2 <- left_join(hypo_do_box_2,fdom_hypo_2,by="time")

hypo_do_box_2 <- hypo_do_box_2 %>% 
  filter(time >= as.POSIXct("2019-01-01"))

# Plot boxplots of various FDOM/CDOM parameters by oxic/anoxic
peakt_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=Date,y=T))+
  geom_point(fdom_hypo,mapping=aes(x=Date,y=T))+
  geom_errorbar(fdom_hypo,mapping=aes(x=Date,y=T,ymin=T-fdom_hypo_sd$T,ymax=T+fdom_hypo_sd$T))+
  ylim(0,0.25)+
  theme_classic(base_size=15)

peakt_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=PeakT,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("Peak T")+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(peakt_time,peakt_box)

ggsave("./Fig_Output/PeakT_oxy.jpg",width=10,height=4,units="in",dpi=320)

peaka_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=Date,y=A))+
  geom_point(fdom_hypo,mapping=aes(x=Date,y=A))+
  geom_errorbar(fdom_hypo,mapping=aes(x=Date,y=A,ymin=A-fdom_hypo_sd$A,ymax=A+fdom_hypo_sd$A))+
  ylim(0,0.35)+
  theme_classic(base_size=15)

peaka_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=PeakA,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("Peak A")+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(peaka_time,peaka_box)

ggsave("./Fig_Output/PeakA_oxy.jpg",width=10,height=4,units="in",dpi=320)

hix_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=Date,y=HIX))+
  geom_point(fdom_hypo,mapping=aes(x=Date,y=HIX))+
  geom_errorbar(fdom_hypo,mapping=aes(x=Date,y=HIX,ymin=HIX-fdom_hypo_sd$HIX,ymax=HIX+fdom_hypo_sd$HIX))+
  ylim(0,6)+
  theme_classic(base_size=15)

hix_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=HIX,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("HIX")+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(hix_time,hix_box)

ggsave("./Fig_Output/HIX_oxy.jpg",width=10,height=4,units="in",dpi=320)

bix_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=Date,y=BIX))+
  geom_point(fdom_hypo,mapping=aes(x=Date,y=BIX))+
  geom_errorbar(fdom_hypo,mapping=aes(x=Date,y=BIX,ymin=BIX-fdom_hypo_sd$BIX,ymax=BIX+fdom_hypo_sd$BIX))+
  ylim(0,0.9)+
  theme_classic(base_size=15)

bix_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=BIX,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("BIX")+
  xlab("Year")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(bix_time,bix_box)

ggsave("./Fig_Output/BIX_oxy.jpg",width=10,height=4,units="in",dpi=320)

### OLD CODE ###
# Plot DO and DOC concentration together w/ color blocking
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2015-10-05"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2016-10-07"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2017-10-30"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(ctd_50_hypo,mapping=aes(x=time,y=DO_mgL,color="DO"))+
  geom_point(ctd_50_hypo,mapping=aes(x=time,y=DO_mgL,color="DO"))+
  xlim(as.POSIXct("2015-03-31"),as.POSIXct("2019-11-20"))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Plot DO and DOC by year (2015-2019)
chem_50_hypo_time <- chem_50_hypo %>% 
  filter(time >= as.POSIXct("2015-01-01")) %>% 
  mutate(year = ifelse(time >= as.POSIXct("2015-01-01") & time <= as.POSIXct("2015-12-31"), "2015",
                       ifelse(time >= as.POSIXct("2016-01-01") & time <= as.POSIXct("2016-12-31"), "2016",
                              ifelse(time >= as.POSIXct("2017-01-01") & time <= as.POSIXct("2017-12-31"), "2017",
                                     ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                                            ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                                                   NA))))))

ggplot(chem_50_hypo_time,mapping=aes(x=time,y=DOC_mgL))+
  geom_point()+
  geom_line()+
  facet_grid(~year)

ctd_50_hypo_time <- ctd_50_hypo %>% 
  filter(time >= as.POSIXct("2015-01-01")) %>% 
    mutate(year = ifelse(time >= as.POSIXct("2015-01-01") & time <= as.POSIXct("2015-12-31"), "2015",
                       ifelse(time >= as.POSIXct("2016-01-01") & time <= as.POSIXct("2016-12-31"), "2016",
                              ifelse(time >= as.POSIXct("2017-01-01") & time <= as.POSIXct("2017-12-31"), "2017",
                                     ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                                            ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                                                   NA))))))

ggplot(ctd_50_hypo_time,mapping=aes(x=time,y=DO_mgL))+
  geom_point()+
  geom_line()+
  facet_grid(~year)

# Separate code source/sink term for different oxygen conditions/year
box_data_full <- left_join(box_data,data_hypo,by="time")

box_data_full <- box_data_full %>% 
  mutate(hox = ifelse(time > as.POSIXct("2015-06-01") & time <= as.POSIXct("2015-06-08"), "2015_off_1",
                      ifelse(time > as.POSIXct("2015-06-08") & time <= as.POSIXct("2015-10-05"), "2015_on_1",
                             ifelse(time >= as.POSIXct("2016-06-01") & time <= as.POSIXct("2016-09-16"), "2016_on_1",
                                    ifelse(time > as.POSIXct("2016-09-16") & time <= as.POSIXct("2016-10-07"),"2016_off_1",
                                           ifelse(time >= as.POSIXct("2017-06-01") & time <= as.POSIXct("2017-10-30"),"2017_on_1",
                                                  ifelse(time >= as.POSIXct("2018-06-01") & time <= as.POSIXct("2018-07-30"),"2018_on_1",
                                                         ifelse(time > as.POSIXct("2018-07-30") & time <= as.POSIXct("2018-10-21"),"2018_off_1",
                                                                ifelse(time >= as.POSIXct("2019-06-01") & time <= as.POSIXct("2019-06-03"), "2019_off_1",
                                                                       ifelse(time > as.POSIXct("2019-06-03") & time <=as.POSIXct("2019-06-17"),"2019_on_1",
                                                                              ifelse(time > as.POSIXct("2019-06-17") & time <= as.POSIXct("2019-07-08"),"2019_off_1",
                                                                                     ifelse(time > as.POSIXct("2019-07-08") & time <= as.POSIXct("2019-07-19"),"2019_on_2",
                                                                                            ifelse(time > as.POSIXct("2019-07-19") & time <= as.POSIXct("2019-08-05"),"2019_off_2",
                                                                                                   ifelse(time > as.POSIXct("2019-08-05") & time <= as.POSIXct("2019-08-19"),"2019_on_3",
                                                                                                          ifelse(time > as.POSIXct("2019-08-19") & time <= as.POSIXct("2019-09-02"),"2019_off_3",
                                                                                                                 ifelse(time > as.POSIXct("2019-09-02") & time <= as.POSIXct("2019-10-23"),"2019_on_4",
                                                                                                                        NA))))))))))))))))

# Plot by on/off schedule?
ggplot(box_data_full,mapping=aes(x=hox,y=DOC_mgL.y))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2)+
  theme_classic(base_size=15)

ggplot(box_data_full,mapping=aes(x=hox,y=DO_mgL))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.2)+
  theme_classic(base_size=15)

# Select for days with DOC data
data_hypo_doc <- data_hypo[!is.na(data_hypo$DOC_mgL),]

# Categorize data as: oxic (>3 mg/L); hypoxic (3>DO>0); anoxic (DO = 0)
data_hypo_doc <- data_hypo_doc %>% 
  mutate(oxy = ifelse(DO_mgL >= 0.7, 'Oxic',
                             ifelse(DO_mgL < 0.7, 'Anoxic', NA))) %>% 
  mutate(year = ifelse(time >= as.POSIXct("2015-01-01") & time <= as.POSIXct("2015-12-31"), "2015",
                       ifelse(time >= as.POSIXct("2016-01-01") & time <= as.POSIXct("2016-12-31"), "2016",
                              ifelse(time >= as.POSIXct("2017-01-01") & time <= as.POSIXct("2017-12-31"), "2017",
                                     ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                                            ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                                                   NA)))))) %>% 
  filter(time >= as.POSIXct("2015-01-01"))

# Also should select for summer stratified period: for now, Jun 01 to Oct 31
data_hypo_doc <- data_hypo_doc %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9))

#############
# Plot temp at depth at 50 + inflow temp
# Where does the inflow go?
# Not convinced it's going to the hypo - looks like temp is between 1.6 and 3.8 m (meaning the
# inflow water likely stays in the epi?)
ggplot()+
  geom_line(ctd_50_depths,mapping=aes(x=time,y=Temp_C,color=as.factor(depth)))+
  geom_line(inflow_daily,mapping=aes(x=DateTime,y=WVWA_Temp_C))+
  theme_classic(base_size=15)

# Plot DOC vs. DO for all depths
ggplot(data_depth,mapping=aes(x=DO_mgL,y=DOC_mgL,color=as.factor(depth)))+
  geom_point()+
  theme_classic(base_size=15)

ggplot(data_hypo,mapping=aes(x=DO_mgL,y=DOC_mgL,color=as.factor(depth)))+
  geom_point()+
  theme_classic(base_size=15)

# Plot long-term DOC (add in oxygenation schedule?)
ggplot(data_hypo_doc,mapping=aes(x=time,y=DOC_mgL,color=as.factor(depth)))+
  annotate(geom="rect",xmin = as.POSIXct("2014-05-06"), xmax = as.POSIXct("2014-06-03"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2014-06-29"), xmax = as.POSIXct("2014-07-29"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2014-08-18"), xmax = as.POSIXct("2014-11-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-05-05"), xmax = as.POSIXct("2015-06-01"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2015-06-08"), xmax = as.POSIXct("2015-10-20"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2016-04-18"), xmax = as.POSIXct("2016-09-16"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2017-04-18"), xmax = as.POSIXct("2017-12-28"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-12-01"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  geom_line(size=1)+
  geom_point(size=3)+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  ylim(0,8)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

# Plot boxplots of DOC in oxic, hypoxic, anoxic conditions for stratified period
ggplot(data_hypo_doc,mapping=aes(year,DOC_mgL,color=oxy))+
  geom_boxplot()+
  theme_classic(base_size=15)
