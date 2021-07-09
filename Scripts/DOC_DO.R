### Script to explore long-term DOC patterns and loading to FCR ###
# 29 Mar 2021, A Hounshell
# Updated: 5 Apr 2021, A Hounshell
# Following CCC thoughts: constrain hypo to 9 m (assume 8m = outflow concentration)
# Don't worry about 'entrainment' - should be included in the outflow variable

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,rMR,lme4,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics)

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
  filter(Site==50 & Depth_m %in% c(6.2,8.0,9.0) | Site==100) %>% 
  filter(DateTime >= as.POSIXct("2018-01-01")) %>% 
  rename(time = DateTime,depth = Depth_m) %>% 
  mutate(loc = ifelse(Site == 50 & depth==6.2, "50 6.2m",
           ifelse(Site==50 & depth==8.0, "50 8m",
                      ifelse(Site==50 & depth==9.0, "50 9m",
                             ifelse(Site==100, "100",NA))))) %>% 
  mutate(DOC_mgL = ifelse(DOC_mgL > 12.0, NA, DOC_mgL)) # Remove outlier at 9 m (12 mg/L!!!)

chem_inputs <- chem_inputs[!is.na(chem_inputs$DOC_mgL),]

# Separate out hypo samples (8-9m) and calculate VW average
chem_hypo <- chem_inputs %>% 
  filter(depth %in% c(8.0,9.0)) %>% 
  select(time,depth,DOC_mgL) %>% 
  pivot_wider(names_from = depth, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_") %>% 
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend")) %>% 
  mutate(VW_Hypo_DOC_mgL = ((DOC_8*1.41*10^4)+(DOC_9*1.95*10^3))/((1.41*10^4)+(1.95*10^3)))

chem_inflow <- chem_inputs %>% 
  filter(Site==50 & depth == c(6.2) | Site==100) 

chem_50_hypo <- chem %>%
  filter(Site==50) %>% 
  filter(Depth_m %in% c(8.0,9.0)) %>% 
  filter(DateTime >= as.POSIXct("2018-01-01")) %>%
  rename(time = DateTime,depth = Depth_m) %>% 
  select(time,depth,DOC_mgL) %>% 
  pivot_wider(names_from = depth, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_")

# Create dataframe with daily time step from 2018-01-01 to 2020-12-31
chem_50_hypo_full <- as.data.frame(seq(as.POSIXct("2018-01-01",tz="EST"),as.POSIXct("2020-12-31",tz="EST"),by="days"))
chem_50_hypo_full <- chem_50_hypo_full %>% 
  rename(time = `seq(as.POSIXct("2018-01-01", tz = "EST"), as.POSIXct("2020-12-31", tz = "EST"), by = "days")`)
chem_50_hypo_full <- left_join(chem_50_hypo_full,chem_50_hypo,by="time")

# Extrapolate and constrain to 2015-2019 data
chem_50_hypo_full <- chem_50_hypo_full %>% 
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend"))

# Calculate VW Hypo mass
# Use same volumes as in L&O-L MS
chem_50_hypo_full <- chem_50_hypo_full %>% 
  mutate(hypo_mg = (DOC_8*1.41*10^4*1000)+(DOC_9*1.95*10^3*1000)) %>% 
  mutate(hypo_vw_mgL = ((DOC_8*1.41*10^4)+(DOC_9*1.95*10^3))/((1.41*10^4)+(1.95*10^3))) %>% 
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
  filter(Depth_m %in% c(6.2)) %>% 
  filter(DateTime >= as.POSIXct("2018-01-01")) %>%
  rename(time = DateTime,depth = Depth_m) %>% 
  select(time,DOC_mgL) %>% 
  rename(DOC_mgL_therm = DOC_mgL) %>% 
  mutate(DOC_mgL_therm = ifelse(DOC_mgL_therm > 12.0, NA, DOC_mgL_therm)) # Remove outlier at 8 m (12 mg/L!!!)

# Add to box_data data sheet and interpolate
box_data <- left_join(box_data,chem_50_therm,by="time")
box_data <- box_data %>% 
  mutate(DOC_mgL_therm = na.fill(na.approx(DOC_mgL_therm,na.rm=FALSE),"extend")) %>%
  group_by(time) %>% 
  summarise_all(mean,na.rm=TRUE)

# Filter and add data from 100
chem_100 <- chem %>% 
  filter(Site==100) %>% 
  filter(DateTime >= as.POSIXct("2018-01-01")) %>%
  rename(time = DateTime, depth = Depth_m) %>% 
  select(time,DOC_mgL) %>% 
  rename(DOC_mgL_100 = DOC_mgL)

box_data <- left_join(box_data,chem_100,by="time")
box_data <- box_data %>% 
  mutate(DOC_mgL_100 = na.fill(na.approx(DOC_mgL_100,na.rm=FALSE),"extend")) %>% 
  group_by(time) %>% 
  summarise_all(mean,na.rm=TRUE)

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
  summarize_all(funs(mean(.,na.rm=TRUE)))

### Add together inflow and box data
inflow_box <- inflow_daily %>% 
  rename(time = DateTime) %>% 
  filter(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2020-12-31")) %>% 
  mutate(year = ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                              ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                                     ifelse(time >= as.POSIXct("2020-01-01") & time <= as.POSIXct("2020-12-31"), "2020",
                                            NA)))) %>% 
  mutate(flow_cms = ifelse(is.na(WVWA_Flow_cms), VT_Flow_cms, WVWA_Flow_cms))

# Check flows
ggplot()+
  geom_line(inflow_box,mapping=aes(x=time,y=WVWA_Flow_cms,color="WVWA"))+
  geom_line(inflow_box,mapping=aes(x=time,y=VT_Flow_cms,color="VT"))+
  theme_classic(base_size = 17)

# Calculate median flow for each year
year_inflow <- inflow_box %>% 
  group_by(year) %>% 
  summarize_all(funs(median),na.rm=TRUE)

# Select flow_cms
inflow_box <- inflow_box %>% 
  select(time,flow_cms)

# Plot all inflow data for SI
ggplot()+
  geom_segment(aes(x = as.POSIXct("2018-01-01"),xend = as.POSIXct("2018-12-31"),y = year_inflow$flow_cms[1],yend = year_inflow$flow_cms[1]), linetype="dashed", color="firebrick", size=1)+
  geom_segment(aes(x = as.POSIXct("2019-01-01"),xend = as.POSIXct("2019-12-31"),y = year_inflow$flow_cms[2],yend = year_inflow$flow_cms[2]), linetype="dashed", color="firebrick", size=1)+
  geom_segment(aes(x = as.POSIXct("2020-01-01"),xend = as.POSIXct("2020-12-02"),y = year_inflow$flow_cms[3],yend = year_inflow$flow_cms[3]), linetype="dashed", color="firebrick", size=1)+
  geom_line(inflow_box,mapping=aes(x=time,y=flow_cms),size=1)+
  ylab(expression(paste("Discharge (m"^3*"s"^-1*")")))+
  xlab("Time")+
  xlim(as.POSIXct("2018-01-01"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggsave("./Fig_Output/Discharge_v3.jpg",width=10,height=4,units="in",dpi=320)

# Calculate WRT for FCR for all years
inflow_box <- inflow_box %>% 
  mutate(wrt_d = 3.1*10^5/flow_cms/60/60/24)

ggplot()+
  geom_line(inflow_box,mapping=aes(x=time,y=wrt_d),size=1)+
  ylab("WRT (d)")+
  xlab("Time")+
  xlim(as.POSIXct("2018-01-01"),as.POSIXct("2020-12-02"))+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

box_data <- left_join(box_data,inflow_box,by="time")
box_data <- box_data %>% 
  mutate(flow_cms = na.fill(na.approx(flow_cms,na.rm=FALSE),"extend"))

### Thinking about Flora data ----
# As a potential 'source' of DOC to the hypo

# Load from EDI: downloaded on 10 June 2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/5/a24f4dbe9f0d856688f257547069d0a3" 
#infile1 <- paste0(getwd(),"/Data/FluoroProbe.csv")
#download.file(inUrl1,infile1,method="curl")

flora <- read.csv("./Data/FluoroProbe.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2017-12-31")) %>% 
  filter(Site == 50)

flora_1 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=0 & Depth_m<9.0) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=0.1)

flora_2 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=4.5 & Depth_m<5.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=5.0)

flora_3 <- flora %>% select(DateTime,Depth_m,TotalConc_ugL) %>% 
  filter(Depth_m>=8.5 & Depth_m<9.5) %>% 
  group_by(DateTime) %>% summarize_all(funs(mean)) %>% arrange(DateTime) %>% 
  mutate(Depth=9.0)

flora_all <- rbind(flora_1,flora_2,flora_3)

flora_all <- flora_all %>% select(-c(Depth_m)) %>% mutate(Station = 50) %>% rename(Date = DateTime, Flora_ugL = TotalConc_ugL) %>% 
  arrange(Date,Station,Depth)

flora_all <- flora_all %>% relocate("Station",.after = "Date") %>% relocate("Depth",.after = "Station")

flora_wc <- flora_1 %>% 
  select(DateTime,TotalConc_ugL) %>% 
  rename(time = DateTime)

# Plot Flora for above the thermocline (> 6.2 m) for each summer stratified period
flora_summer <- flora_all %>% 
  filter(Date>as.POSIXct("2018-01-15")&Date<as.POSIXct("2020-12-02")) %>% 
  filter(Depth == 0.1) %>% 
  mutate(month = month(Date)) %>% 
  filter(month %in% c(6,7,8,9,10))

# Plot!
chla_18 <- ggplot(flora_summer,mapping=aes(x=Date,y=Flora_ugL))+
  geom_line()+
  geom_point()+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  ylab(expression(paste("Chla ("*mu*"g L"^-1*")")))+
  xlab("2018")+
  theme_classic(base_size = 15)

chla_19 <- ggplot(flora_summer,mapping=aes(x=Date,y=Flora_ugL))+
  geom_line()+
  geom_point()+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  ylab(expression(paste("Chla ("*mu*"g L"^-1*")")))+
  xlab("2019")+
  theme_classic(base_size = 15)

chla_20 <- ggplot(flora_summer,mapping=aes(x=Date,y=Flora_ugL))+
  geom_line()+
  geom_point()+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  ylab(expression(paste("Chla ("*mu*"g L"^-1*")")))+
  xlab("2020")+
  theme_classic(base_size = 15)

ggarrange(chla_18,chla_19,chla_20,ncol=1,nrow=3,common.legend = TRUE, labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Epi_Flora.jpg",width=10,height=12,units="in",dpi=320)

### Calculate 'sink/source' term following Gerling et al. 2016; Kruger et al. 2020
box_data <- box_data %>% 
  mutate(j_mgs = dMdt_mgs-(flow_cms*DOC_mgL_100*1000*0.26)+(flow_cms*DOC_mgL_therm*1000*0.26)) %>% 
  mutate(j_kgd = j_mgs*60*60*24/1000/1000)

# Select for 2015-2019 during the stratified period
box_data_sel <- box_data %>% 
  filter(time>as.POSIXct("2018-01-15")&time<as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9,10))

# Plot - Figure 4: DOC Box model results
box_18 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2018-08-06"), xmax = as.POSIXct("2018-08-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic based on CTD casts!
  annotate(geom="rect",xmin = as.POSIXct("2018-08-16"), xmax = as.POSIXct("2018-10-21"), ymin=-Inf, ymax=Inf,alpha = 0.3)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_line(box_data_sel,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26,color="Inflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26,color="Outflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=j_kgd,color="j_kgd"),size=2)+
  ylab(expression(paste("(kg C d"^-1*")")))+
  xlab("2018")+
  scale_color_manual(breaks=c('Outflow_kgd','dMdt_kgd','Inflow_kgd','j_kgd'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Outflow","dM/dt","Inflow","Source/Sink"))+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-11-01"))+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

box_19 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_line(box_data_sel,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26,color="Inflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26,color="Outflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=j_kgd,color="j_kgd"),size=2)+
  xlab("2019")+
  ylab(expression(paste("(kg C d"^-1*")")))+
  scale_color_manual(breaks=c('Outflow_kgd','dMdt_kgd','Inflow_kgd','j_kgd'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Outflow","dM/dt","Inflow","Source/Sink"))+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-11-01"))+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

box_20 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-01"), xmax = as.POSIXct("2020-07-06"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2020-07-08"), xmax = as.POSIXct("2020-07-17"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-23"), xmax = as.POSIXct("2020-08-06"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-09-15"), xmax = as.POSIXct("2020-09-25"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_line(box_data_sel,mapping=aes(x=time,y=dMdt_mgs*60*60*24/1000/1000,color="dMdt_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26,color="Inflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26,color="Outflow_kgd"),size=2)+
  geom_line(box_data_sel,mapping=aes(x=time,y=j_kgd,color="j_kgd"),size=2)+
  xlab("2020")+
  ylab(expression(paste("(kg C d"^-1*")")))+
  scale_color_manual(breaks=c('Outflow_kgd','dMdt_kgd','Inflow_kgd','j_kgd'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Outflow","dM/dt","Inflow","Source/Sink"))+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-11-01"))+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(box_18,box_19,box_20,ncol=1,nrow=3,common.legend = TRUE, labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig4_scaledFlows_v2.jpg",width=10,height=12,units="in",dpi=320)

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
ctd_50_depths$Cond_uScm <- as.numeric(ctd_50_depths$Cond_uScm)
ctd_50_depths$DO_mgL <- as.numeric(ctd_50_depths$DO_mgL)

# Plot DO for each depth
ggplot(ctd_50_depths,mapping=aes(x=time,y=DO_mgL,color=as.factor(depth)))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 2)+
  xlim(as.POSIXct("2018-01-01"),as.POSIXct("2020-12-31"))+
  theme_classic(base_size=15)

ggplot(chem,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Depth_m)))+
  geom_line()+
  geom_point()+
  xlim(as.POSIXct("2018-01-01"),as.POSIXct("2020-12-31"))+
  theme_classic(base_size=15)

# Calculate %DO Saturation using rMR
ctd_do_calc <- ctd_50_depths %>% 
  select(time,Temp_C,DO_mgL,Cond_uScm,depth)

ctd_do_calc <- na.omit(ctd_do_calc)
  
ctd_do_calc <- ctd_do_calc %>% 
  mutate(calc_DO_pSat = DO.saturation(DO.mgl = ctd_do_calc$DO_mgL,temp.C = ctd_do_calc$Temp_C,
                                      elevation.m = 509, bar.press = NULL, bar.units = NULL,
                                      ctd_do_calc$Cond_uScm,salinity.units = "uS")*100)

ctd_50_depths_2 <- left_join(ctd_50_depths,ctd_do_calc,by=c("time","depth","Temp_C","DO_mgL","Cond_uScm"))

ctd_50_depths_2 <- ctd_50_depths_2 %>% 
  mutate(calc_DO_pSat = ifelse(is.na(calc_DO_pSat),DO_pSat,calc_DO_pSat))

ctd_50_depths_2$calc_DO_pSat <- as.numeric(ctd_50_depths_2$calc_DO_pSat)

ctd_50_hypo <- ctd_50_depths_2 %>% 
  filter(depth == 9.0 | depth == 8.0) %>% 
  filter(time >= as.POSIXct("2018-01-01")) %>% 
  select(time,Temp_C,DO_mgL,calc_DO_pSat,depth) %>% 
  pivot_wider(names_from = depth, values_from = c(Temp_C,DO_mgL,calc_DO_pSat), values_fil = NA, values_fn = mean)

ctd_50_hypo <- na.omit(ctd_50_hypo)

ctd_50_hypo <- ctd_50_hypo %>% 
  mutate(vw_temp_C = ((Temp_C_8*1.41*10^4)+(Temp_C_9*1.95*10^3))/((1.41*10^4)+(1.95*10^3))) %>% 
  mutate(vw_DO_mgL = ((DO_mgL_8*1.41*10^4)+(DO_mgL_9*1.95*10^3))/((1.41*10^4)+(1.95*10^3))) %>% 
  mutate(vw_pSat_DO = ((calc_DO_pSat_8*1.41*10^4)+(calc_DO_pSat_9*1.95*10^3))/((1.41*10^4)+(1.95*10^3)))

ctd_50_hypo$vw_DO_mgL <- as.numeric(ctd_50_hypo$vw_DO_mgL)

# Calculate DOC mineralization following Carey et al. 2018
# Need DOC concentration, temperature, amount O2 added
box_data_remin <- left_join(box_data,ctd_50_hypo, by = "time")

box_data_remin <- box_data_remin %>% 
  select(-Temp_C_8,-Temp_C_9,-DO_mgL_8,-DO_mgL_9,-calc_DO_pSat_8,-calc_DO_pSat_9) %>% 
  mutate(vw_temp_C = na.fill(na.approx(vw_temp_C,na.rm=FALSE),"extend"))

box_data_remin <- box_data_remin %>% 
  mutate(oxy_mgL = ifelse(time >= "2018-04-23" & time <= "2018-07-20", (15*1e6)/(1046476.8),
                          ifelse(time >= "2019-06-03" & time <= "2019-06-18", (25*1e6)/(1133683.2),
                                 ifelse(time >= "2019-07-08" & time <= "2019-07-20", (25*1e6)/(1133683.2),
                                      ifelse(time >= "2019-08-05" & time <= "2019-08-20", (25*1e6)/(1133683.2),
                                             ifelse(time >= "2019-09-02" & time <= "2019-12-02", (25*1e6)/(1133683.2),
                                                    ifelse(time >= "2020-06-29" & time <= "2020-07-06", (25*1e6)/(1133683.2),
                                                           ifelse(time >= "2020-07-06" & time <= "2020-07-12", (25*1e6)/(1133683.2),
                                                                  ifelse(time >= "2020-07-13" & time <= "2020-07-22", (25*1e6)/(1133683.2),
                                                                         ifelse(time >= "2020-07-23" & time <= "2020-12-02", (25*1e6)/(1133683.2),
                                                                                0))))))))))

# Plot to check oxygenation schedule
ggplot(box_data_remin,mapping=aes(x=time,y=oxy_mgL))+
  geom_line()

# Calculate DOC mineralization using equations in Carey et al. 2018
# Start by calculating Rh-o2
box_data_remin <- box_data_remin %>% 
  mutate(rh_oxy = 1+(0.46*(oxy_mgL/7.015*10^7))/(0.2+(oxy_mgL/7.015*10^7))) %>% 
  mutate(rh_doc = hypo_vw_mgL*0.0093*1.07^(vw_temp_C-20)*rh_oxy) %>% 
  mutate(rh_perc = rh_doc/1e6/j_kgd*100) %>% 
  select(-vw_temp_C,-vw_DO_mgL,-vw_pSat_DO)

# Combine box model data w/ DO data
hypo_do_box_full <- left_join(box_data_remin,ctd_50_hypo,by="time")

anoxia_time <- hypo_do_box_full %>% 
  select(time,vw_DO_mgL) %>% 
  mutate(vw_DO_mgL = na.fill(na.approx(vw_DO_mgL,na.rm=FALSE),"extend")) %>% 
  mutate(anoxia = ifelse(vw_DO_mgL< 1.0, 1, 0)) %>% 
  mutate(anoxia_time_d = 0) %>% 
  mutate(oxygenation = ifelse(time >= "2018-04-23" & time <= "2018-07-20", 1,
                          ifelse(time >= "2019-06-03" & time <= "2019-06-18", 1,
                                 ifelse(time >= "2019-07-08" & time <= "2019-07-20", 1,
                                        ifelse(time >= "2019-08-05" & time <= "2019-08-20", 1,
                                               ifelse(time >= "2019-09-02" & time <= "2019-12-02", 1,
                                                      ifelse(time >= "2020-06-29" & time <= "2020-07-06", 1,
                                                             ifelse(time >= "2020-07-06" & time <= "2020-07-12", 1,
                                                                    ifelse(time >= "2020-07-13" & time <= "2020-07-22", 1,
                                                                           ifelse(time >= "2020-07-23" & time <= "2020-12-02", 1,
                                                                                  0))))))))))

for (i in 1:length(anoxia_time$time)){
  if (anoxia_time$anoxia[i] == 1){
    anoxia_time$anoxia_time_d[i] = anoxia_time$anoxia_time_d[i-1]+1
  } else {
    anoxia_time$anoxia_time_d[i] = 0
  }
}

anoxia_time <- anoxia_time %>% 
  select(time,anoxia_time_d,oxygenation)

# Combine box model data w/ DO data and constrain to measured time points
hypo_do_box <- left_join(ctd_50_hypo,box_data_remin,by="time")
hypo_do_box <- left_join(hypo_do_box,anoxia_time,by="time")
hypo_do_box$vw_DO_mgL <- as.numeric(hypo_do_box$vw_DO_mgL)

hypo_do_box <- hypo_do_box %>% 
  filter(time >= as.POSIXct("2018-01-15") & time <= as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9,10)) %>% 
  mutate(oxy = ifelse(vw_DO_mgL >= 1.0, 'Oxic',
                      ifelse(vw_DO_mgL < 1.0, 'Anoxic',
                             NA))) %>% 
  mutate(year = ifelse(time >= as.POSIXct("2018-01-01") & time <= as.POSIXct("2018-12-31"), "2018",
                       ifelse(time >= as.POSIXct("2019-01-01") & time <= as.POSIXct("2019-12-31"), "2019",
                              ifelse(time >= as.POSIXct("2020-01-01") & time <= as.POSIXct("2020-12-31"), "2020",
                                     NA))))

hypo_do_box <- hypo_do_box[!is.na(hypo_do_box$vw_DO_mgL),]

# Interpolate NA values
hypo_do_box <- hypo_do_box %>% 
  mutate(vw_DO_mgL = na.fill(na.approx(vw_DO_mgL,na.rm=FALSE),"extend")) %>% 
  mutate(vw_temp_C = na.fill(na.approx(vw_temp_C,na.rm=FALSE),"extend"))

# Add Box plot for WRT?
ggplot(hypo_do_box,mapping=aes(oxy,wrt_d,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("WRT (d)")+
  xlab("")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

# Think about plotting DOC remineralization
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2018-04-23"), xmax = as.POSIXct("2018-07-30"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29"), xmax = as.POSIXct("2020-09-11"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2020-09-25"), xmax = as.POSIXct("2020-12-02"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(box_data_remin,mapping=aes(x=time,y=rh_doc/1e6),size=1)

# Plot discharge for each summer stratified period
q_18 <- ggplot()+
  geom_line(inflow_box,mapping=aes(x=time,y=flow_cms),size=0.8)+
  geom_point(hypo_do_box,mapping=aes(x=time,y=flow_cms),size=3)+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  ylab(expression(paste("Discharge (m"^3*"s"^-1*")")))+
  xlab("2018")+
  theme_classic(base_size = 15)

q_19 <- ggplot()+
  geom_line(inflow_box,mapping=aes(x=time,y=flow_cms),size=0.8)+
  geom_point(hypo_do_box,mapping=aes(x=time,y=flow_cms),size=3)+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  ylab(expression(paste("Discharge (m"^3*"s"^-1*")")))+
  xlab("2019")+
  theme_classic(base_size = 15)

q_20 <- ggplot()+
  geom_line(inflow_box,mapping=aes(x=time,y=flow_cms),size=0.8)+
  geom_point(hypo_do_box,mapping=aes(x=time,y=flow_cms),size=3)+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  ylab(expression(paste("Discharge (m"^3*"s"^-1*")")))+
  xlab("2020")+
  theme_classic(base_size = 15)

ggarrange(q_18,q_19,q_20,ncol=1,nrow=3,common.legend = TRUE, labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Discharge_v2.jpg",width=10,height=12,units="in",dpi=320)

# Thinking about time since anoxia as a driver?
# Calculate linear trend line for each year
# Remove values when duration of anoxia = 0
hypo_do_box_sel <- hypo_do_box %>% 
  filter(anoxia_time_d > 0)

fits <- lmList(j_kgd ~ anoxia_time_d | year, data = hypo_do_box_sel)
fits <- lmList(hypo_vw_mgL ~ anoxia_time_d | year, data = hypo_do_box_sel)

anoxia_doc_zero <- hypo_do_box %>% 
  filter(anoxia_time_d == 0) %>% 
  ggplot(mapping=aes(x=anoxia_time_d,y=hypo_vw_mgL,color=year))+
  geom_boxplot()+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  ylab(expression(paste("VW DOC (mg L"^-1*")")))+
  xlab("Duration of anoxia (d)")+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

anoxia_doc <- ggplot(hypo_do_box_sel,mapping=aes(x=anoxia_time_d,y=hypo_vw_mgL,color=year))+
  geom_point(size=3)+
  ylab(expression(paste("VW DOC (mg L"^-1*")")))+
  xlab("Duration of anoxia (d)")+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

anoxia_loading_zero <- hypo_do_box %>% 
  filter(anoxia_time_d == 0) %>% 
  ggplot(mapping=aes(x=anoxia_time_d,y=j_kgd,color=year))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_boxplot()+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  ylab(expression(paste("Internal DOC loading (kg C d"^-1*")")))+
  xlab("Duration of anoxia (d)")+
  ylim(-10,10)+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

anoxia_loading <- ggplot(hypo_do_box_sel,mapping=aes(anoxia_time_d,j_kgd,colour=year))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(size=3)+
  #geom_smooth(method="lm")+
  ylab(expression(paste("Internal DOC loading (kg C d"^-1*")")))+
  xlab("Duration of anoxia (d)")+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  ylim(-10,10)+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

ggarrange(anoxia_loading_zero, anoxia_loading,common.legend = TRUE,ncol=2,nrow=1,widths=c(1,2),
          labels = c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/DurationAnoxia_v3.jpg",width=10,height=4,units="in",dpi=320)

# Fig. 5: Plot DOC and DOC source/sink term by oxic vs. anoxic
# Figure 3 - [DOC] timeseries for 6.2 m, VW Hypo, 100 from 2018-01-01 to 2020-12-02
# Constrain to stratified period (June-Oct)
chem_hypo_strat <- chem_hypo %>% 
  filter(time >= as.POSIXct("2018-01-15") & time <= as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9,10))

chem_inflow_strat <- chem_inflow %>% 
  filter(time >= as.POSIXct("2018-01-15") & time <= as.POSIXct("2020-12-02")) %>% 
  mutate(month = month(time)) %>% 
  filter(month %in% c(6,7,8,9,10))

# Plot DOC by year
# 2018
doc_18 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2018-08-06"), xmax = as.POSIXct("2018-08-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic based on CTD casts!
  annotate(geom="rect",xmin = as.POSIXct("2018-08-16"), xmax = as.POSIXct("2018-10-21"), ymin=-Inf, ymax=Inf,alpha = 0.3)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_line(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=0.8)+
  geom_point(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=3)+
  geom_line(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=0.8)+
  geom_point(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=3)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Time")+
  scale_color_manual(breaks=c('50 6.2m','VW_Hypo_DOC','100'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Thermo","VW Hypo","Inflow"))+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  ylim(0,5)+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

# 2019
doc_19 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=0.8)+
  geom_point(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=3)+
  geom_line(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=0.8)+
  geom_point(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=3)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Time")+
  scale_color_manual(breaks=c('50 6.2m','VW_Hypo_DOC','100'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Thermo","VW Hypo","Inflow"))+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  ylim(0,5)+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

# 2020
doc_20 <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-01"), xmax = as.POSIXct("2020-07-06"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2020-07-08"), xmax = as.POSIXct("2020-07-17"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-23"), xmax = as.POSIXct("2020-08-06"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-09-15"), xmax = as.POSIXct("2020-09-25"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=0.8)+
  geom_point(chem_inflow_strat,mapping=aes(x=time,y=DOC_mgL,color=as.factor(loc)),size=3)+
  geom_line(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=0.8)+
  geom_point(chem_hypo_strat,mapping=aes(x=time,y=VW_Hypo_DOC_mgL,color="VW_Hypo_DOC"),size=3)+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Time")+
  scale_color_manual(breaks=c('50 6.2m','VW_Hypo_DOC','100'),values=c("#7EBDC2","#393E41","#F0B670"),labels=c("Thermo","VW Hypo","Inflow"))+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  ylim(0,5)+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

# [DOC]
doc_full <- ggplot(hypo_do_box,mapping=aes(oxy,hypo_vw_mgL,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("VW DOC (mg L"^-1*")")))+
  xlab("")+
  ylim(0,5)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

doc <- ggplot(hypo_do_box,mapping=aes(year,hypo_vw_mgL,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("VW DOC (mg L"^-1*")")))+
  xlab("Year")+
  ylim(0,5)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(ggarrange(doc_18,doc_19,doc_20,ncol=3,labels = c("A.","B.","C."),common.legend = TRUE,font.label = list(face="plain",size=15)),
          ggarrange(doc_full,doc,ncol=2,labels = c("D.","E."),common.legend = TRUE,font.label=list(face="plain",size=15)),
          nrow=2)

ggsave("./Fig_Output/Fig3_v5.jpg",width=10,height=8,units="in",dpi=320)

## Plot VW Hypo DOC vs. DO for each year/summer stratified period
do_doc <- ggplot(hypo_do_box,mapping=aes(vw_DO_mgL,hypo_vw_mgL,colour=year))+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_point(size=3)+
  ylab(expression(paste("VW DOC (mg L"^-1*")")))+
  xlab(expression(paste("VW DO (mg L"^-1*")")))+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

do_loading <- ggplot(hypo_do_box,mapping=aes(vw_DO_mgL,j_kgd,colour=year))+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_point(size=3)+
  ylab(expression(paste("Internal DOC loading (kg C d"^-1*")")))+
  xlab(expression(paste("VW DO (mg L"^-1*")")))+
  scale_color_manual(breaks=c('2018','2019','2020'),values=c("#7EBDC2","#393E41","#F0B670"))+
  theme_classic(base_size = 17)+
  theme(legend.title=element_blank())

ggarrange(do_doc,do_loading,common.legend = TRUE,ncol=2,nrow=1,labels = c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/DO_DOC.jpg",width=10,height=5,units="in",dpi=320)

# Source/sink term
jterm_full <- ggplot(hypo_do_box,mapping=aes(oxy,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Internal loading (kg C d"^-1*")")))+
  xlab("")+
  ylim(-10,11)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

jterm <- ggplot(hypo_do_box,mapping=aes(year,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Internal loading (kg C d"^-1*")")))+
  xlab("Year")+
  ylim(-10,11)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(jterm_full ,jterm,ncol=2,nrow=1,common.legend=TRUE,labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig6_HypoBypass_v2.jpg",width=10,height=5,units="in",dpi=320)

# Plot temp by year (as box plots?)
# Remove wonky casts: 2019-04-29, 2019-05-31
#ctd_50_hypo <- ctd_50_hypo %>% 
#  mutate(Temp_C = ifelse(Temp_C >= 16, NA, Temp_C)) %>% 
#  mutate(DO_mgL = ifelse(Temp_C >= 16, NA, DO_mgL))

# Also want to add in timeseries of DO for the full time period
#ctd_50_hypo <- ctd_50_hypo[!is.na(ctd_50_hypo$DO_mgL),]
temp_box <- ggplot(hypo_do_box,mapping=aes(year,vw_temp_C,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Temp (C"^o*")")))+
  xlab("Year")+
  ylim(0,15)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

# Plot DO as box plots
do_box <- ggplot(hypo_do_box,mapping=aes(year,vw_DO_mgL,colour=oxy))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("DO (mg L"^-1*")")))+
  xlab("Year")+
  ylim(0,8)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

# Separate timeseries by year and focus on stratified period (June-Oct)
# 2018
do_18 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_DO_mgL))+
  geom_hline(yintercept = 1, linetype="dashed")+
  annotate(geom="rect",xmin = as.POSIXct("2018-08-06"), xmax = as.POSIXct("2018-08-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic based on CTD casts!
  annotate(geom="rect",xmin = as.POSIXct("2018-08-16"), xmax = as.POSIXct("2018-10-21"), ymin=-Inf, ymax=Inf,alpha = 0.3)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("DO (mg L"^-1*")")))+
  xlab("Year")+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  ylim(0,8)+
  theme_classic(base_size = 17)

# 2019
do_19 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_DO_mgL))+
  geom_hline(yintercept = 1, linetype="dashed")+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("DO (mg L"^-1*")")))+
  xlab("Year")+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  ylim(0,8)+
  theme_classic(base_size = 17)

# 2020
do_20 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_DO_mgL))+
  geom_hline(yintercept = 1, linetype="dashed")+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-01"), xmax = as.POSIXct("2020-07-06"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2020-07-08"), xmax = as.POSIXct("2020-07-17"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-23"), xmax = as.POSIXct("2020-08-06"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-09-15"), xmax = as.POSIXct("2020-09-25"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("DO (mg L"^-1*")")))+
  xlab("Year")+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  ylim(0,8)+
  theme_classic(base_size = 17)

# Separate temp by year
# 2018
temp_18 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_temp_C))+
  annotate(geom="rect",xmin = as.POSIXct("2018-08-06"), xmax = as.POSIXct("2018-08-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic based on CTD casts!
  annotate(geom="rect",xmin = as.POSIXct("2018-08-16"), xmax = as.POSIXct("2018-10-21"), ymin=-Inf, ymax=Inf,alpha = 0.3)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("Temp ("^o*"C)")))+
  xlab("Year")+
  ylim(0,15)+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  theme_classic(base_size = 17)

# 2019
temp_19 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_temp_C))+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("Temp ("^o*"C)")))+
  xlab("Year")+
  ylim(0,15)+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  theme_classic(base_size = 17)

# 2020
temp_20 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_temp_C))+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-01"), xmax = as.POSIXct("2020-07-06"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2020-07-08"), xmax = as.POSIXct("2020-07-17"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-23"), xmax = as.POSIXct("2020-08-06"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-09-15"), xmax = as.POSIXct("2020-09-25"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab(expression(paste("Temp ("^o*"C)")))+
  xlab("Year")+
  ylim(0,15)+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  theme_classic(base_size = 17)

ggarrange(ggarrange(temp_18,temp_19,temp_20,do_18,do_19,do_20,nrow=2,ncol=3,labels=c("A.","B.","C.","D.","E.","F."),font.label = list(face="plain",size=15)),
          ggarrange(temp_box,do_box,ncol=3,common.legend=TRUE,labels=c("G.","H."),font.label=list(face="plain",size=15)),
          nrow=2,heights=c(2,1))
  
ggsave("./Fig_OutPut/DO_temp_Year_v5.jpg",width=12,height=11,units="in",dpi=320)

# Plot %DO sat for SI
hypo_do_box$vw_pSat_DO <- as.numeric(hypo_do_box$vw_pSat_DO)
ctd_50_hypo$vw_pSat_DO <- as.numeric(ctd_50_hypo$vw_pSat_DO)

dosat_box <- ggplot(hypo_do_box,mapping=aes(year,vw_pSat_DO,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("DO % Saturation")+
  xlab("Year")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

# Separate by year
# 2018
dosat_18 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_pSat_DO))+
  annotate(geom="rect",xmin = as.POSIXct("2018-08-06"), xmax = as.POSIXct("2018-08-10"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic based on CTD casts!
  annotate(geom="rect",xmin = as.POSIXct("2018-08-16"), xmax = as.POSIXct("2018-10-21"), ymin=-Inf, ymax=Inf,alpha = 0.3)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab("DO % Saturation")+
  xlab("Year")+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-10-31"))+
  theme_classic(base_size = 17)

# 2019
dosat_19 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_pSat_DO))+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab("DO % Saturation")+
  xlab("Year")+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  theme_classic(base_size = 17)

# 2020
dosat_20 <- ggplot(hypo_do_box,mapping=aes(x=time,y=vw_pSat_DO))+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-01"), xmax = as.POSIXct("2020-07-06"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2020-07-08"), xmax = as.POSIXct("2020-07-17"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-23"), xmax = as.POSIXct("2020-08-06"),ymin=-Inf,ymax=Inf,alpha = 0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-09-15"), xmax = as.POSIXct("2020-09-25"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed")+ #Turnover FCR; operationally defined
  geom_line(size=0.8)+
  geom_point(size=3)+
  ylab("DO % Saturation")+
  xlab("Year")+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-10-31"))+
  theme_classic(base_size = 17)

ggarrange(dosat_18,dosat_19,dosat_20,dosat_box,nrow=2,ncol=2,
          labels=c("A.","B.","C.","D."),font.label = list(face="plain",size=15))

ggsave("./Fig_OutPut/DOSat_v3.jpg",width=10,height=9,units="in",dpi=320)

# Calculate median DO for each year/oxygenation period
do_med <- hypo_do_box %>% 
  #select(-oxy) %>% 
  group_by(oxy,year) %>% 
  summarize_all(funs(max),na.rm=TRUE) %>% 
  select(oxy,year,j_kgd)

#Plot Inflow and dM/dt by year
dm_dt <- ggplot(hypo_do_box,mapping=aes(year,dMdt_mgs*60*60*24/1000/1000,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("dM/dt (kg C d"^-1*")")))+
  xlab("Year")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

inflow <- ggplot(hypo_do_box,mapping=aes(year,flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Inflow (kg C d"^-1*")")))+
  xlab("")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

outflow <- ggplot(hypo_do_box,mapping=aes(year,flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Outflow (kg C d"^-1*")")))+
  xlab("")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

jterm <- ggplot(hypo_do_box,mapping=aes(year,j_kgd,colour=oxy))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("Internal loading (kg C d"^-1*")")))+
  xlab("Year")+
  ylim(-15,11)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(inflow,outflow,dm_dt,jterm,nrow=2,ncol=2,common.legend = TRUE,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/dmtodt_inflow_hypobypass_v3.jpg",width=10,height=8,units="in",dpi=320)

### Generate AR model for DOC internal loading ----
# Including: AR term, DO, Temp, Anoxia Duration, and Oxygenation
# Following Howard et al. 2021; AGU GHG presentation

# First add in Flora data (full water column mean)
hypo_do_box <- left_join(hypo_do_box,flora_wc,by="time")

# Select data need for AR Model: 2018
ar_model_2018 <- hypo_do_box %>% 
  filter(year==2018) %>% 
  select(time,vw_temp_C,vw_DO_mgL,j_kgd,flow_cms,TotalConc_ugL)

# Create list of each Monday from June 2018 to end of October 2018
mon_2018 <- seq(as.Date("2018-06-04"),as.Date("2018-10-31"),"7 days")
mon_2018 <- as.POSIXct(strptime(mon_2018,"%Y-%m-%d", tz="EST"))
mon_2018 <- as.data.frame(mon_2018)
mon_2018 <- mon_2018 %>% 
  rename(dates_2018 = mon_2018)

#thurs_2018 <- seq(as.Date("2018-06-07"),as.Date("2018-10-31"),"7 days")
#thurs_2018 <- as.POSIXct(strptime(thurs_2018,"%Y-%m-%d", tz="EST"))
#thurs_2018 <- as.data.frame(thurs_2018)
#thurs_2018 <- thurs_2018 %>% 
#  rename(dates_2018 = thurs_2018)

#dates_2018 <- rbind(mon_2018,thurs_2018)

dates_2018 <- mon_2018

dates_2018 <- dates_2018 %>% 
  arrange(dates_2018) %>% 
  mutate(time = dates_2018)

ar_model_2018 <- full_join(ar_model_2018,dates_2018,by="time")

ar_model_2018 <- ar_model_2018 %>% 
  arrange(time) %>% 
  mutate(vw_temp_C = na.fill(na.approx(vw_temp_C,na.rm=FALSE),"extend")) %>% 
  mutate(vw_DO_mgL = na.fill(na.approx(vw_DO_mgL,na.rm=FALSE),"extend")) %>% 
  mutate(j_kgd = na.fill(na.approx(j_kgd,na.rm=FALSE),"extend")) %>% 
  mutate(flow_cms = na.fill(na.approx(flow_cms,na.rm=FALSE),"extend")) %>% 
  mutate(TotalConc_ugL = na.fill(na.approx(TotalConc_ugL,na.rm=FALSE),"extend")) %>% 
  group_by(time) %>% 
  summarise_all(mean,na.rm=TRUE)

ar_model_2018 <- left_join(ar_model_2018,anoxia_time,by="time")

ar_model_2018 <- ar_model_2018[!is.na(ar_model_2018$dates_2018),]

# Select data need for AR Model: 2019
ar_model_2019 <- hypo_do_box %>% 
  filter(year==2019) %>% 
  select(time,vw_temp_C,vw_DO_mgL,j_kgd,flow_cms,TotalConc_ugL)

# Create list of each Monday from June 2019 to end of October 2019
mon_2019 <- seq(as.Date("2019-06-03"),as.Date("2019-10-31"),"7 days")
mon_2019 <- as.POSIXct(strptime(mon_2019,"%Y-%m-%d", tz="EST"))
mon_2019 <- as.data.frame(mon_2019)
mon_2019 <- mon_2019 %>% 
  rename(dates_2019 = mon_2019)

#thurs_2019 <- seq(as.Date("2019-06-06"),as.Date("2019-10-31"),"7 days")
#thurs_2019 <- as.POSIXct(strptime(thurs_2019,"%Y-%m-%d", tz="EST"))
#thurs_2019 <- as.data.frame(thurs_2019)
#thurs_2019 <- thurs_2019 %>% 
#  rename(dates_2019 = thurs_2019)

#dates_2019 <- rbind(mon_2019,thurs_2019)

dates_2019 <- mon_2019

dates_2019 <- dates_2019 %>% 
  arrange(dates_2019) %>% 
  mutate(time = dates_2019)

ar_model_2019 <- full_join(ar_model_2019,dates_2019,by="time")

ar_model_2019 <- ar_model_2019 %>% 
  arrange(time) %>% 
  mutate(vw_temp_C = na.fill(na.approx(vw_temp_C,na.rm=FALSE),"extend")) %>% 
  mutate(vw_DO_mgL = na.fill(na.approx(vw_DO_mgL,na.rm=FALSE),"extend")) %>% 
  mutate(j_kgd = na.fill(na.approx(j_kgd,na.rm=FALSE),"extend"))%>% 
  mutate(flow_cms = na.fill(na.approx(flow_cms,na.rm=FALSE),"extend")) %>% 
  mutate(TotalConc_ugL = na.fill(na.approx(TotalConc_ugL,na.rm=FALSE),"extend")) %>% 
  group_by(time) %>% 
  summarise_all(mean,na.rm=TRUE)

ar_model_2019 <- left_join(ar_model_2019,anoxia_time,by="time")

ar_model_2019 <- ar_model_2019[!is.na(ar_model_2019$dates_2019),]

# Select data need for AR Model: 2020
ar_model_2020 <- hypo_do_box %>% 
  filter(year==2020) %>% 
  select(time,vw_temp_C,vw_DO_mgL,j_kgd,flow_cms,TotalConc_ugL)

# Create list of each Monday from June 2020 to end of October 2020
mon_2020 <- seq(as.Date("2020-06-01"),as.Date("2020-10-31"),"7 days")
mon_2020 <- as.POSIXct(strptime(mon_2020,"%Y-%m-%d", tz="EST"))
mon_2020 <- as.data.frame(mon_2020)
mon_2020 <- mon_2020 %>% 
  rename(dates_2020 = mon_2020)

#thurs_2020 <- seq(as.Date("2020-06-04"),as.Date("2020-10-31"),"7 days")
#thurs_2020 <- as.POSIXct(strptime(thurs_2020,"%Y-%m-%d", tz="EST"))
#thurs_2020 <- as.data.frame(thurs_2020)
#thurs_2020 <- thurs_2020 %>% 
#  rename(dates_2020 = thurs_2020)

#dates_2020 <- rbind(mon_2020,thurs_2020)

dates_2020 <- mon_2020

dates_2020 <- dates_2020 %>% 
  arrange(dates_2020) %>% 
  mutate(time = dates_2020)

ar_model_2020 <- full_join(ar_model_2020,dates_2020,by="time")

ar_model_2020 <- ar_model_2020 %>% 
  arrange(time) %>% 
  mutate(vw_temp_C = na.fill(na.approx(vw_temp_C,na.rm=FALSE),"extend")) %>% 
  mutate(vw_DO_mgL = na.fill(na.approx(vw_DO_mgL,na.rm=FALSE),"extend")) %>% 
  mutate(j_kgd = na.fill(na.approx(j_kgd,na.rm=FALSE),"extend"))%>% 
  mutate(flow_cms = na.fill(na.approx(flow_cms,na.rm=FALSE),"extend")) %>% 
  mutate(TotalConc_ugL = na.fill(na.approx(TotalConc_ugL,na.rm=FALSE),"extend")) %>% 
  group_by(time) %>% 
  summarise_all(mean,na.rm=TRUE)

ar_model_2020 <- left_join(ar_model_2020,anoxia_time,by="time")

ar_model_2020 <- ar_model_2020[!is.na(ar_model_2020$dates_2020),]

### Check for autocorrelation among DOC internal loading parameters (for each year)

# DOC internal loading, 2018: no lag term!
plot(ar_model_2018$j_kgd,type="b")
lag1.plot(ar_model_2018$j_kgd,10)
PlotACF(ar_model_2018$j_kgd)
acf2(ar_model_2018$j_kgd,na.action=na.pass)
pacf(ar_model_2018$j_kgd,na.action=na.pass)

# DOC internal loading 2019: first lag term!
plot(ar_model_2019$j_kgd,type="b")
lag1.plot(ar_model_2019$j_kgd,10)
PlotACF(ar_model_2019$j_kgd)
acf2(ar_model_2019$j_kgd,na.action=na.pass)
pacf(ar_model_2019$j_kgd,na.action=na.pass)

# Add lag term to 2019 data
xlag1 = lag(ar_model_2019$j_kgd,1)
y = cbind(ar_model_2019$j_kgd,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("j","j_kgd_ARLag1")
ar_model_2019 <- cbind(ar_model_2019,y)
ar_model_2019 <- ar_model_2019 %>% select(-j)

# DOC internal loading 2020: no lag term!
plot(ar_model_2020$j_kgd,type="b")
lag1.plot(ar_model_2020$j_kgd,10)
PlotACF(ar_model_2020$j_kgd)
acf2(ar_model_2020$j_kgd,na.action=na.pass)
pacf(ar_model_2020$j_kgd,na.action=na.pass)

# Plot ACF and PACF for each year
pdf("./Fig_Output/ACF_Plots.pdf", width=12, height=8)

acf2(ar_model_2018$j_kgd,na.action=na.pass)
acf2(ar_model_2019$j_kgd,na.action=na.pass)
acf2(ar_model_2020$j_kgd,na.action=na.pass)

dev.off()

# Check for collinearity among variables for each year and standardize to z-scores
# 2018
ar_model_2018_2 <- ar_model_2018 %>% 
  select(-time,-dates_2018,-oxygenation,-flow_cms,-TotalConc_ugL) %>% 
  scale()

# Keep variables with r<0.80 (pushing the limits here!)
# Removed oxygenation (keep: temp, DO, anoxia_time)
chart.Correlation(ar_model_2018_2,histogram=TRUE,method=c("spearman"))

# 2019
ar_model_2019_2 <- ar_model_2019 %>% 
  select(-time,-dates_2019,-flow_cms,-TotalConc_ugL) %>% 
  scale()

# Keep all variables??
chart.Correlation(ar_model_2019_2,histogram=TRUE,method=c("spearman"))

# 2020
ar_model_2020_2 <- ar_model_2020 %>% 
  select(-time,-dates_2020,-flow_cms,-TotalConc_ugL,-anoxia_time_d) %>% 
  scale()

# Remove anoxia_time (keep: temp, DO, oxygenation)
chart.Correlation(ar_model_2020_2,histogram=TRUE,method=c("spearman"))

## Generate AR models!
# 2018
ar_model_2018_2 <- as.data.frame(ar_model_2018_2)
model_2018 <- glm(j_kgd ~ vw_DO_mgL + vw_temp_C + anoxia_time_d, data = ar_model_2018_2, family = gaussian,
                  na.action = 'na.fail')

glm_2018 <- dredge(model_2018,rank="AICc")

select_glm_2018 <- subset(glm_2018,delta<2)

# 2019
ar_model_2019_2 <- as.data.frame(ar_model_2019_2)

ar_model_2019_2 <- ar_model_2019_2[complete.cases(ar_model_2019_2),]

model_2019 <- glm(j_kgd ~ j_kgd_ARLag1 + vw_DO_mgL + vw_temp_C + anoxia_time_d + oxygenation, data = ar_model_2019_2, 
                  family = gaussian,
                  na.action = 'na.fail')

glm_2019 <- dredge(model_2019,rank="AICc")

select_glm_2019 <- subset(glm_2019,delta<2)

# 2020
ar_model_2020_2 <- as.data.frame(ar_model_2020_2)
model_2020 <- glm(j_kgd ~ vw_DO_mgL + vw_temp_C + oxygenation, data = ar_model_2020_2, family = gaussian,
                  na.action = 'na.fail')

glm_2020 <- dredge(model_2020,rank="AICc")

select_glm_2020 <- subset(glm_2020,delta<2)


### Let's start thinking about DOM quality - what metrics to use? ----
# a254? HIX? BIX? Peak C? Peak T? PARAFAC?
# Load in data
fdom <- read_csv("./EDI_2021/20210511_OpticalData.csv") %>% 
  filter(Reservoir == "FCR" & Dilution %in% c(1,2)) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST"))) %>% 
  na.omit(fdom)

fdom_hypo <- fdom %>% 
  filter(Depth_m == 9.0) %>% 
  group_by(DateTime,Depth_m) %>% 
  summarize_all(funs(mean(.,na.rm=TRUE))) %>% 
  rename(time = DateTime,PeakT=T,PeakA=A)

fdom_hypo_sd <- fdom %>% 
  filter(Depth_m == 9.0) %>% 
  group_by(DateTime,Depth_m) %>% 
  summarize_all(funs(sd(.,na.rm=TRUE))) %>% 
  rename(time = DateTime,PeakT=T,PeakA=A)

fdom_inflow <- fdom %>% 
  filter(Site == 100) %>% 
  group_by(DateTime) %>% 
  summarize_all(funs(mean(.,na.rm=TRUE))) %>% 
  rename(time = DateTime,PeakT=T,PeakA=A)

fdom_inflow_sd <- fdom %>% 
  filter(Site == 100) %>% 
  group_by(DateTime) %>% 
  summarize_all(funs(sd(.,na.rm=TRUE))) %>% 
  rename(time = DateTime,PeakT=T,PeakA=A)

# Combine with DO data
hypo_do_box_2 <- left_join(hypo_do_box,fdom_hypo,by="time")

hypo_do_box_2 <- hypo_do_box_2 %>% 
  filter(time >= as.POSIXct("2019-01-01")&time<=as.POSIXct("2019-12-31")) %>% 
  mutate(SUVA = a254_m/2.303/hypo_vw_mgL)

suva_hypo <- hypo_do_box_2 %>% 
  select(time,vw_DO_mgL,a254_m,hypo_vw_mgL,SUVA,year,month,oxy)

suva_hypo <- suva_hypo[!is.na(suva_hypo$SUVA),]

# Plot boxplots of various FDOM/CDOM parameters by oxic/anoxic
peakt_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_inflow,mapping=aes(x=time,y=PeakT,color="Inflow"),size=0.8)+
  geom_point(fdom_inflow,mapping=aes(x=time,y=PeakT,color="Inflow"),size=3)+
  geom_errorbar(fdom_inflow,mapping=aes(x=time,y=PeakT,ymin=PeakT-fdom_inflow_sd$PeakT,ymax=PeakT+fdom_inflow_sd$PeakT,color="Inflow"))+
  geom_line(fdom_hypo,mapping=aes(x=time,y=PeakT,color="9 m"),size=0.8)+
  geom_point(fdom_hypo,mapping=aes(x=time,y=PeakT,color="9 m"),size=3)+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=PeakT,ymin=PeakT-fdom_hypo_sd$PeakT,ymax=PeakT+fdom_hypo_sd$PeakT,color="9 m"))+
  scale_color_manual(breaks=c('9 m','Inflow'),values=c("#393E41","#F0B670"))+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  ylim(0,0.25)+
  ylab("Peak T (R.F.U.)")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

peakt_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=PeakT,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("Peak T (R.F.U.)")+
  xlab("Year")+
  ylim(0,0.25)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

peaka_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-04"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"), xmax = as.POSIXct("2019-07-09"), ymin = -Inf, ymax=Inf, alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2019-07-29"), xmax = as.POSIXct("2019-08-05"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  annotate(geom="rect",xmin = as.POSIXct("2019-08-22"),xmax = as.POSIXct("2019-09-11"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Anoxic
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_inflow,mapping=aes(x=time,y=PeakA,color="Inflow"),size=0.8)+
  geom_point(fdom_inflow,mapping=aes(x=time,y=PeakA,color="Inflow"),size=3)+
  geom_errorbar(fdom_inflow,mapping=aes(x=time,y=PeakA,ymin=PeakA-fdom_inflow_sd$PeakA,ymax=PeakA+fdom_inflow_sd$PeakA,color="Inflow"))+
  geom_line(fdom_hypo,mapping=aes(x=time,y=PeakA,color="9 m"),size=0.8)+
  geom_point(fdom_hypo,mapping=aes(x=time,y=PeakA, color="9 m"),size=3)+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=PeakA,ymin=PeakA-fdom_hypo_sd$PeakA,ymax=PeakA+fdom_hypo_sd$PeakA, color="9 m"))+
  scale_color_manual(breaks=c('9 m','Inflow'),values=c("#393E41","#F0B670"))+
  ylim(0,0.35)+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-10-31"))+
  xlab("")+
  ylab("Peak A (R.F.U.)")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

peaka_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=PeakA,colour=oxy))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("Peak A (R.F.U.)")+
  xlab("")+
  ylim(0,0.35)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(peaka_time,peaka_box,peakt_time,peakt_box,nrow=2,ncol=2,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig6_v4.jpg",width=11,height=7,units="in",dpi=320)

# Try plotting some other parameters as well
a254_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=time,y=a254_m),size=0.8)+
  geom_point(fdom_hypo,mapping=aes(x=time,y=a254_m),size=3)+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=a254_m,ymin=a254_m-fdom_hypo_sd$a254_m,ymax=a254_m+fdom_hypo_sd$a254_m))+
  ylim(0,30)+
  xlab("")+
  ylab(expression(paste("a254 (m"^-1*")")))+
  theme_classic(base_size=17)

a254_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=a254_m,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("a254 (m"^-1*")")))+
  xlab("")+
  theme_classic(base_size=17)+
  ylim(0,30)+
  theme(legend.title=element_blank())

# Calculate SUVA
suva_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(suva_hypo,mapping=aes(x=time,y=SUVA),size=0.8)+
  geom_point(suva_hypo,mapping=aes(x=time,y=SUVA),size=3)+
  ylim(0,8)+
  xlab("")+
  ylab(expression(paste("SUVA (L mg"^-1*" m"^-1*")")))+
  theme_classic(base_size=17)

suva_box <- ggplot(suva_hypo,mapping=aes(x=year,y=SUVA,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab(expression(paste("SUVA (L mg"^-1*" m"^-1*")")))+
  xlab("")+
  theme_classic(base_size=17)+
  ylim(0,8)+
  theme(legend.title=element_blank())

ggarrange(a254_time,a254_box,suva_time,suva_box,common.legend=TRUE,nrow=2,ncol=2,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig6_abs.jpg",width=10,height=8,units="in",dpi=320)

# Plot indicators (HIX, BIX)
hix_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=time,y=HIX))+
  geom_point(fdom_hypo,mapping=aes(x=time,y=HIX))+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=HIX,ymin=HIX-fdom_hypo_sd$HIX,ymax=HIX+fdom_hypo_sd$HIX))+
  ylim(0,6)+
  xlab("")+
  ylab("HIX")+
  theme_classic(base_size=17)

hix_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=HIX,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("HIX")+
  xlab("")+
  ylim(0,6)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

bix_time <- ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=time,y=BIX))+
  geom_point(fdom_hypo,mapping=aes(x=time,y=BIX))+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=BIX,ymin=BIX-fdom_hypo_sd$BIX,ymax=BIX+fdom_hypo_sd$BIX))+
  ylim(0,0.9)+
  xlab("")+
  ylab("BIX")+
  theme_classic(base_size=17)

bix_box <- ggplot(hypo_do_box_2,mapping=aes(x=year,y=BIX,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("BIX")+
  xlab("")+
  ylim(0,0.9)+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(hix_time,hix_box,bix_time,bix_box,common.legend=TRUE,nrow=2,ncol=2,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig6_IX.jpg",width=10,height=8,units="in",dpi=320)

# Check Peak M?
ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = as.POSIXct("2019-11-20"),ymin=-Inf,ymax=Inf,alpha=0.3)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed")+ #Turnover FCR
  geom_line(fdom_hypo,mapping=aes(x=time,y=M))+
  geom_point(fdom_hypo,mapping=aes(x=time,y=M))+
  geom_errorbar(fdom_hypo,mapping=aes(x=time,y=M,ymin=M-fdom_hypo_sd$M,ymax=M+fdom_hypo_sd$M))+
  #ylim(0,0.25)+
  ylab("Peak M (R.F.U.)")+
  theme_classic(base_size=17)

ggplot(hypo_do_box_2,mapping=aes(x=year,y=M,colour=oxy))+
  geom_boxplot()+
  scale_color_manual(breaks=c('Anoxic','Oxic'),values=c("#CD5C5C","#598BAF"))+
  ylab("Peak M (R.F.U.)")+
  xlab("Year")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

# Find range for FDOM parameters
med_fdom <- hypo_do_box_2 %>% 
  select(oxy,PeakA,PeakT) %>% 
  group_by(oxy) %>% 
  summarize_all(funs(min,max),na.rm=TRUE)

# Calculate median for various groups of data ----
all_med <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,oxy,hypo_vw_mgL) %>% 
  group_by(year,oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE)

col_order <- c("year","oxy","vw_temp_C","vw_DO_mgL","hypo_vw_mgL","j_kgd","inflow_kgd","outflow_kgd","dMdt_kgd")

all_med <- all_med[,col_order]

all_med_fdom <- hypo_do_box_2 %>% 
  select(year,oxy,PeakA,PeakT) %>% 
  group_by(year,oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE)

all_med <- left_join(all_med,all_med_fdom,by=c("year","oxy"))

year_med <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,hypo_vw_mgL) %>%
  group_by(year) %>% 
  summarize_all(funs(median),na.rm=TRUE) %>% 
  mutate(oxy = "All")

year_med <- year_med[,col_order]

year_med_fdom <- hypo_do_box_2 %>% 
  select(year,PeakA,PeakT) %>% 
  group_by(year) %>% 
  summarize_all(funs(median),na.rm=TRUE) %>% 
  mutate(oxy = "All")

year_med <- left_join(year_med,year_med_fdom,by=c("year","oxy"))

oxy_med <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,oxy,hypo_vw_mgL) %>%
  group_by(oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE) %>% 
  mutate(year = "All")

oxy_med <- oxy_med[,col_order]

oxy_med_fdom <- hypo_do_box_2 %>% 
  select(oxy,PeakA,PeakT) %>% 
  group_by(oxy) %>% 
  summarize_all(funs(median),na.rm=TRUE) %>% 
  mutate(year = "All")

oxy_med <- left_join(oxy_med,oxy_med_fdom,by=c("year","oxy"))

all_year_med <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,hypo_vw_mgL) %>%
  summarize_all(funs(median),na.rm=TRUE) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_med <- all_year_med[,col_order]

all_year_med_fdom <- hypo_do_box_2 %>% 
  select(PeakA,PeakT) %>% 
  summarize_all(funs(median)) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_med <- left_join(all_year_med,all_year_med_fdom,by=c("year","oxy"))

all_all_med <- rbind(all_year_med,oxy_med,year_med,all_med)

all_all_med <- all_all_med %>% 
  arrange(year,oxy)

all_all_med$PeakA[11] <- NA
all_all_med$PeakT[11] <- NA
all_all_med$PeakA[12] <- NA
all_all_med$PeakT[12] <- NA

# Calculate 25th percentile for various groups of data
all_quan_25 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,oxy,hypo_vw_mgL) %>%
  group_by(year,oxy) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE)))

col_order <- c("year","oxy","vw_temp_C","vw_DO_mgL","hypo_vw_mgL","j_kgd","inflow_kgd","outflow_kgd","dMdt_kgd")

all_quan_25 <- all_quan_25[,col_order]

all_quan_25_fdom <- hypo_do_box_2 %>% 
  select(year,oxy,PeakA,PeakT) %>% 
  group_by(year,oxy) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE)))

all_quan_25 <- left_join(all_quan_25,all_quan_25_fdom,by=c("year","oxy"))

year_quan_25 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,hypo_vw_mgL) %>%
  group_by(year) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>% 
  mutate(oxy = "All")

year_quan_25 <- year_quan_25[,col_order]

year_quan_25_fdom <- hypo_do_box_2 %>% 
  select(year,PeakA,PeakT) %>% 
  group_by(year) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>%  
  mutate(oxy = "All")

year_quan_25 <- left_join(year_quan_25,year_quan_25_fdom,by=c("year","oxy"))

oxy_quan_25 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,oxy,hypo_vw_mgL) %>%
  group_by(oxy) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>% 
  mutate(year = "All")

oxy_quan_25 <- oxy_quan_25[,col_order]

oxy_quan_25_fdom <- hypo_do_box_2 %>% 
  select(oxy,PeakA,PeakT) %>% 
  group_by(oxy) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>% 
  mutate(year = "All")

oxy_quan_25 <- left_join(oxy_quan_25,oxy_quan_25_fdom,by=c("year","oxy"))

all_year_quan_25 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,hypo_vw_mgL) %>%
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_quan_25 <- all_year_quan_25[,col_order]

all_year_quan_25_fdom <- hypo_do_box_2 %>% 
  select(PeakA,PeakT) %>% 
  summarize_all(funs(quantile(.,.25,na.rm=TRUE))) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_quan_25 <- left_join(all_year_quan_25,all_year_quan_25_fdom,by=c("year","oxy"))

all_all_quan_25 <- rbind(all_year_quan_25,oxy_quan_25,year_quan_25,all_quan_25)

all_all_quan_25 <- all_all_quan_25 %>% 
  arrange(year,oxy)

all_all_quan_25$PeakA[11] <- NA
all_all_quan_25$PeakT[11] <- NA
all_all_quan_25$PeakA[12] <- NA
all_all_quan_25$PeakT[12] <- NA

# Calculate 75th percentile for various groups of data
all_quan_75 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,oxy,hypo_vw_mgL) %>%
  group_by(year,oxy) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE)))

col_order <- c("year","oxy","vw_temp_C","vw_DO_mgL","hypo_vw_mgL","j_kgd","inflow_kgd","outflow_kgd","dMdt_kgd")

all_quan_75 <- all_quan_75[,col_order]

all_quan_75_fdom <- hypo_do_box_2 %>% 
  select(year,oxy,PeakA,PeakT) %>% 
  group_by(year,oxy) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE)))

all_quan_75 <- left_join(all_quan_75,all_quan_75_fdom,by=c("year","oxy"))

year_quan_75 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,year,hypo_vw_mgL) %>%
  group_by(year) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>% 
  mutate(oxy = "All")

year_quan_75 <- year_quan_75[,col_order]

year_quan_75_fdom <- hypo_do_box_2 %>% 
  select(year,PeakA,PeakT) %>% 
  group_by(year) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>%  
  mutate(oxy = "All")

year_quan_75 <- left_join(year_quan_75,year_quan_75_fdom,by=c("year","oxy"))

oxy_quan_75 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,oxy,hypo_vw_mgL) %>%
  group_by(oxy) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>% 
  mutate(year = "All")

oxy_quan_75 <- oxy_quan_75[,col_order]

oxy_quan_75_fdom <- hypo_do_box_2 %>% 
  select(oxy,PeakA,PeakT) %>% 
  group_by(oxy) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>% 
  mutate(year = "All")

oxy_quan_75 <- left_join(oxy_quan_75,oxy_quan_75_fdom,by=c("year","oxy"))

all_year_quan_75 <- hypo_do_box %>% 
  mutate(outflow_kgd = flow_cms*DOC_mgL_therm*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(inflow_kgd = flow_cms*DOC_mgL_100*1000*60*60*24/1000/1000*0.26) %>% 
  mutate(dMdt_kgd = dMdt_mgs*60*60*24/1000/1000) %>% 
  select(vw_temp_C,vw_DO_mgL,j_kgd,inflow_kgd,outflow_kgd,dMdt_kgd,hypo_vw_mgL) %>%
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_quan_75 <- all_year_quan_75[,col_order]

all_year_quan_75_fdom <- hypo_do_box_2 %>% 
  select(PeakA,PeakT) %>% 
  summarize_all(funs(quantile(.,.75,na.rm=TRUE))) %>% 
  mutate(oxy = "All") %>% 
  mutate(year = "All")

all_year_quan_75 <- left_join(all_year_quan_75,all_year_quan_75_fdom,by=c("year","oxy"))

all_all_quan_75 <- rbind(all_year_quan_75,oxy_quan_75,year_quan_75,all_quan_75)

all_all_quan_75 <- all_all_quan_75 %>% 
  arrange(year,oxy)

all_all_quan_75$PeakA[11] <- NA
all_all_quan_75$PeakT[11] <- NA
all_all_quan_75$PeakA[12] <- NA
all_all_quan_75$PeakT[12] <- NA

really_all <- left_join(all_all_med,all_all_quan_25,by=c("year","oxy"))

really_all <- left_join(really_all,all_all_quan_75,by=c("year","oxy"))

really_all_round <- really_all %>% 
  select(-year,-oxy) %>% 
  round(.,digits=2)

really_all_round$Temp_med <- paste(really_all_round$vw_temp_C.x,really_all_round$vw_temp_C.y,really_all_round$vw_temp_C,sep=";")
really_all_round$DO_med <- paste(really_all_round$vw_DO_mgL.x,really_all_round$vw_DO_mgL.y,really_all_round$vw_DO_mgL,sep=";")
really_all_round$hypo_vw_mgL_med <- paste(really_all_round$hypo_vw_mgL.x,really_all_round$hypo_vw_mgL.y,really_all_round$hypo_vw_mgL,sep=";")
really_all_round$j_med <- paste(really_all_round$j_kgd.x,really_all_round$j_kgd.y,really_all_round$j_kgd,sep=";")
really_all_round$inflow_med <- paste(really_all_round$inflow_kgd.x,really_all_round$inflow_kgd.y,really_all_round$inflow_kgd,sep=";")
really_all_round$outflow_med <- paste(really_all_round$outflow_kgd.x,really_all_round$outflow_kgd.y,really_all_round$outflow_kgd,sep=";")
really_all_round$dMdt_med <- paste(really_all_round$dMdt_kgd.x,really_all_round$dMdt_kgd.y,really_all_round$dMdt_kgd,sep=";")
really_all_round$PeakA_med <- paste(really_all_round$PeakA.x,really_all_round$PeakA.y,really_all_round$PeakA,sep=";")
really_all_round$PeakT_med <- paste(really_all_round$PeakT.x,really_all_round$PeakT.y,really_all_round$PeakT,sep=";")

really_all <- really_all %>% 
  select(year,oxy)

really_all_round <- really_all_round %>% 
  select(Temp_med,DO_med,hypo_vw_mgL_med,j_med,inflow_med,outflow_med,dMdt_med,PeakA_med,PeakT_med)

really_all <- cbind(really_all,really_all_round)

# Export out and format as a table in excel
write_csv(really_all,"./Fig_Output/20210628_SITable_quantiles.csv")


### OLD CODE ----
# Loading in and plotting various FDOM/CDOM parameters (using non-EDI data!)
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

cdom_hypo_2 <- cdom_hypo %>% 
  rename(time = Date)

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
