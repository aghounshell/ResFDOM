### Script to make graphs for 2021 SFS Presentation
### A Hounshell, 22 Apr 2021

# Load in packages
pacman::p_load(tidyverse,ggplot2,ggpubr)

# DOC data: Downloaded from EDI 19 Apr 21 ----
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98" 
#infile1 <- paste0(getwd(),"/Data/chem.csv")
#download.file(inUrl1,infile1,method="curl")

doc <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2019-04-29") & DateTime <= as.POSIXct("2020-04-01")) %>% 
  filter(Site %in% c('100','200','20','30','45','50','1') & Depth_m == 0.1) %>% 
  mutate(Loc = paste(Reservoir,Site)) %>% 
  filter(Rep == 1)

doc_fcr <- doc %>% 
  filter(Reservoir == "FCR") %>% 
  filter(Site != '1')

doc_bvr <- doc %>% 
  filter(Reservoir == "BVR")

doc_all_res <- rbind(doc_fcr,doc_bvr)

# Designate BVR as left and right arm
doc_all_res <- doc_all_res %>% 
  mutate(Site2 = ifelse(Loc == "BVR 20", "20_L",
                        ifelse(Loc == "BVR 30", "20_R",
                               ifelse(Loc == "BVR 1", "30",
                                      Site))))

# Boxplots?
doc_all_res$Site2 <- factor(doc_all_res$Site2, levels=c("100","200","20_L","20_R","20","30","45","50"))

doc_box <- ggplot(doc_all_res,mapping=aes(as.factor(Site2),DOC_mgL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)

# What about for nitrate and srp?
ggplot(doc_all_res,mapping=aes(as.factor(Site2),NO3NO2_ugL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)

ggplot(doc_all_res,mapping=aes(as.factor(Site2),SRP_ugL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)

# Also load discharge and calculate loading ----
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/454/4/a18421fd2e95c15d6f97009d5fef3e59" 
infile1 <- paste0(getwd(),"/Data/2019_Continuum_Discharge.csv")
download.file(inUrl1,infile1,method="curl")
# Downloaded data on 29 Apr 2021

qcharge <- read_csv("./Data/2019_Continuum_Discharge.csv") %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%Y-%m-%d", tz="EST"))) %>% 
  filter(Site %in% c(100,200)) %>% 
  rename(DateTime = Date)

loads <- left_join(qcharge,doc,by=c("DateTime","Reservoir","Site"))

loads_2 <- loads %>% 
  mutate(N_load_gd = Flow_cms*NO3NO2_ugL*1000*60*60*24/(1*10^6)) %>% 
  mutate(P_load_gd = Flow_cms*SRP_ugL*1000*60*60*24/(1*10^6)) %>% 
  mutate(C_load_kgd = Flow_cms*DOC_mgL*60*60*24/1000)

# Boxplot
ggplot(loads_2,mapping=aes(as.factor(Site),N_load_gd,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("N load (g d"^-1*")")))+
  theme_classic(base_size=15)

ggplot(loads_2,mapping=aes(as.factor(Site),P_load_gd,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("P load (g d"^-1*")")))+
  theme_classic(base_size=15)

ggplot(loads_2,mapping=aes(as.factor(Site),C_load_kgd,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("C load (kg d"^-1*")")))+
  theme_classic(base_size=15)

# Fluorescence (HIX, BIX) ----
fl <- read_csv("./Data/20210210_ResultsFiles_ResEEMs2019_RAW.csv") %>% 
  filter(Dilution %in% c(1,2)) %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%m/%d/%Y", tz="EST"))) %>% 
  filter(Date >= as.POSIXct("2019-04-29") & Date <= as.POSIXct("2020-04-01")) %>% 
  filter(Station %in% c('100','200','20','30','45','50','1') & Depth == 0.1) %>% 
  mutate(Loc = paste(Reservoir,Station))

fl_fcr <- fl %>% 
  filter(Reservoir == "FCR") %>% 
  filter(Station != '1')

fl_bvr <- fl %>% 
  filter(Reservoir == "BVR")

fl_all_res <- rbind(fl_fcr,fl_bvr)

# Designate BVR as left and right arm
fl_all_res <- fl_all_res %>% 
  mutate(Site2 = ifelse(Loc == "BVR 20", "20_L",
                        ifelse(Loc == "BVR 30", "20_R",
                               ifelse(Loc == "BVR 1", "30",
                                      Station))))

# Boxplots?
fl_all_res$Site2 <- factor(fl_all_res$Site2, levels=c("100","200","20_L","20_R","20","30","45","50"))

ggplot(fl_all_res,mapping=aes(as.factor(Site2),HIX,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab("HIX")+
  theme_classic(base_size=15)

ggplot(fl_all_res,mapping=aes(as.factor(Site2),BIX,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab("BIX")+
  theme_classic(base_size=15)

peakt_box <- ggplot(fl_all_res,mapping=aes(as.factor(Site2),T,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab("Peak T")+
  theme_classic(base_size=15)

peaka_box <- ggplot(fl_all_res,mapping=aes(as.factor(Site2),A,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab("Peak A")+
  theme_classic(base_size=15)

peakc_box <- ggplot(fl_all_res,mapping=aes(as.factor(Site2),C,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab("Peak C")+
  theme_classic(base_size=15)

ggarrange(doc_box,peakc_box,peaka_box,peakt_box,ncol=2,nrow=2,common.legend = TRUE)

ggsave("./Fig_Output/DOC_spatial_RCDays.png",width=10,height=7,units="in",dpi=320)

# Think about adding in filtered Chla data ----
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/555/1/93e2c69f314809705ea21d9244eb368d" 
infile1 <- paste0(getwd(),"/Data/chla_master_df_dt.csv")
download.file(inUrl1,infile1,method="curl")
# Downloaded data on 29 Apr 2021

chla <- read.csv("./Data/chla_master_df_dt.csv", header=T) %>%
  select(ï..Reservoir:Chla_ugL) %>%
  rename(Reservoir = ï..Reservoir) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-04-01") & DateTime < as.POSIXct("2020-04-01")) %>% 
  filter(Site %in% c(100, 200, 20, 30, 45, 50, 1) & Depth_m == 0.1) %>% 
  mutate(Site2 = ifelse(Reservoir == "BVR" & Site == "20", "20_L",
                        ifelse(Reservoir == "BVR" & Site == "30", "20_R",
                               ifelse(Reservoir == "BVR" & Site == "1", "30",
                                      Site)))) %>% 
  filter(Site2 %in% c("100","200","20_L","20_R","20","30","45","50"))

chla$Site2 <- factor(chla$Site2, levels=c("100","200","20_L","20_R","20","30","45","50"))

ggplot(chla,mapping=aes(as.factor(Site2),Chla_ugL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("Chla (",mu,"g L"^-1*")")))+
  theme_classic(base_size=15)

# Load in GHGs? ----
# From Ryan's dropbox (RC day NOT published to EDI!)
ghg <- read_csv("./Data/RCDays_GHG.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-04-01") & DateTime < as.POSIXct("2019-10--0="))

ghg$ch4_umolL <- as.numeric(ghg$ch4_umolL)
ghg$co2_umolL <- as.numeric(ghg$co2_umolL)

ghg <- ghg %>% 
  group_by(DateTime,Depth_m,Reservoir) %>% 
  summarize_all(mean,na.rm=TRUE) %>% 
  filter(Depth_m %in% c("0.1","B100","B200","B20","B30","B45","B50","B01","F100","F200",
                        "F20","F30","F45","F50")) %>% 
  mutate(Site2 = ifelse(Reservoir == "FCR" & Depth_m == "0.1", "50",
                        ifelse(Reservoir == "BVR" & Depth_m=="0.1", "50",
                               ifelse(Depth_m == "B20", "20_L",
                                      ifelse(Depth_m == "B30", "20_R",
                                             ifelse(Depth_m == "B01", "30",
                                                    ifelse(Depth_m == "B100", "100",
                                                           ifelse(Depth_m == "B200", "200",
                                                                  ifelse(Depth_m == "B45", "45",
                                                                         ifelse(Depth_m == "F100", "100",
                                                                                ifelse(Depth_m == "F200", "200",
                                                                                       ifelse(Depth_m == "F20", "20",
                                                                                              ifelse(Depth_m == "F30", "30",
                                                                                                     ifelse(Depth_m == "F45","45",
                                                                                                            ifelse(Depth_m == "B50", "50",
                                                                                                            Depth_m))))))))))))))) %>% 
  
  filter(Site2 %in% c("100","200","20_L","20_R","20","30","45","50"))

ghg$Site2 <- factor(ghg$Site2, levels=c("100","200","20_L","20_R","20","30","45","50"))

# Plot
ggplot(ghg,mapping=aes(as.factor(Site2),ch4_umolL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("CH4 (",mu,"mol L"^-1*")")))+
  theme_classic(base_size=15)

ggplot(ghg,mapping=aes(as.factor(Site2),co2_umolL,color=Reservoir))+
  geom_boxplot()+
  scale_color_manual(breaks=c('BVR','FCR'),
                     values=c("#1B9E77","#1D1A31"))+
  xlab("Site")+
  ylab(expression(paste("CO2 (",mu,"mol L"^-1*")")))+
  theme_classic(base_size=15)

# Old code to look at temporal variability and explore synchronoy ----

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

doc_fcr <- completeFun(doc_fcr,"DOC_mgL") %>% 
  select(Site, DateTime, Depth_m, DOC_mgL) %>% 
  rename(DOC_mgL_FCR = DOC_mgL)

doc_bvr <- completeFun(doc_bvr,"DOC_mgL") %>% 
  select(Site, DateTime, Depth_m, DOC_mgL) %>% 
  rename(DOC_mgL_BVR = DOC_mgL)

doc_fcr_res <- doc_fcr %>% 
  filter(Site %in% c('20','30','45','50'))

doc_fcr_in <- doc_fcr %>% 
  filter(Site %in% c('100','200'))

doc_bvr_res <- doc_bvr %>% 
  filter(Site %in% c('20','30','45','50','1')) %>% 
  mutate(Site2 = ifelse(Site == '30', '20_R',
                        ifelse(Site == '1', '30',
                               Site)))

doc_bvr_left <- doc_bvr_res %>% 
  filter(Site2!="20_R") %>% 
  select(DateTime,Depth_m,DOC_mgL_BVR,Site2) %>% 
  rename(Site = Site2)

doc_bvr_left$Site <- as.double(doc_bvr_left$Site)

doc_bvr_right <- doc_bvr_res %>% 
  filter(Site2 == "20_R") %>% 
  mutate(Site = ifelse(Site2 == "20_R", "20", NA)) %>% 
  rename(DOC_mgL_BVR_r = DOC_mgL_BVR)

doc_bvr_right$Site <- as.double(doc_bvr_right$Site)

doc_bvr_in <- doc_bvr %>% 
  filter(Site %in% c('100','200'))

# Combine for synchrony (start with basics - don't change any names)
syn_all <- full_join(doc_fcr,doc_bvr,by=c("Site","DateTime","Depth_m"))

# Plot?
ggplot(syn_all,mapping=aes(DOC_mgL_FCR,DOC_mgL_BVR,colour=as.factor(Site)))+
  geom_point()+
  theme_classic(base_size=15)

# Probably makes more sense to focus on the in reservoir sites (from loctic to lentic)?
syn_res <- full_join(doc_fcr_res,doc_bvr_res,by=c("Site","DateTime","Depth_m"))

# Plot?
ggplot(syn_res,mapping=aes(DOC_mgL_FCR,DOC_mgL_BVR,colour=as.factor(Site)))+
  geom_point()+
  theme_classic(base_size=15)

# Combine for synchrony - separating into left and right arms of BVR
syn_arms <- full_join(doc_fcr_res,doc_bvr_left,by=c("Site","DateTime","Depth_m"))
syn_arms <- full_join(syn_arms,doc_bvr_right,by=c("Site","DateTime","Depth_m"))

cor.test(syn_arms$DOC_mgL_FCR,syn_arms$DOC_mgL_BVR,method="pearson")

# Plot?
ggplot(syn_arms,mapping=aes(DOC_mgL_FCR,DOC_mgL_BVR,colour=as.factor(Site)))+
  geom_point()+
  theme_classic(base_size=15)

# Plot?
# All from both reservoirs
ggplot()+
  geom_line(doc_fcr,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_fcr,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  geom_line(doc_bvr,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1,linetype="longdash")+
  geom_point(doc_bvr,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  theme_classic(base_size=15)

ggplot()+
  geom_line(doc_fcr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_fcr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  geom_line(doc_bvr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1,linetype="longdash")+
  geom_point(doc_bvr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  theme_classic(base_size=15)

fcr_res <- 
  ggplot()+
  geom_line(doc_fcr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_fcr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  scale_color_manual(breaks=c('20','30','45','50'),
                     values=c("#1B9E77","#CA5902","#7570B3","#1D1A31"))+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  labs(fill = "Site")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_res <- ggplot()+
  geom_line(doc_bvr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_bvr_res,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  scale_color_manual(breaks=c('20','1','30','45','50'),labels=c('20','30-L','30-R','45','50'),
                     values=c("#1B9E77","#FD8C35","#CA5902","#7570B3","#1D1A31"))+
  ylim(0,9)+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

fcr_in <- 
ggplot()+
  geom_line(doc_fcr_in,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_fcr_in,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0,9)+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_in <- ggplot()+
  geom_line(doc_bvr_in,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=1)+
  geom_point(doc_bvr_in,mapping=aes(x=DateTime,y=DOC_mgL,color=as.factor(Site)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0,9)+
  xlab("Date")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(fcr_in,bvr_in,fcr_res,bvr_res,nrow=2,ncol=2,common.legend = FALSE)

ggsave("./Fig_Output/DOC_RCDays.png",width=10,height=7,units="in",dpi=320)

# Fluorescence (HIX, BIX) ----
fl_mean <- fl %>% 
  group_by(Date,Reservoir,Station) %>% 
  summarize_all(funs(mean),na.rm=TRUE)

fl_sd <- fl %>% 
  group_by(Date,Reservoir,Station) %>% 
  summarize_all(funs(sd),na.rm=TRUE)

fl_fcr_res <- fl %>% 
  filter(Station %in% c('20','30','45','50') & Reservoir == "FCR")

fl_fcr_in <- fl %>% 
  filter(Station %in% c('100','200') & Reservoir == "FCR")

fl_bvr_res <- fl %>% 
  filter(Station %in% c('20','30','45','50','1') & Reservoir == "BVR")

fl_bvr_in <- fl %>% 
  filter(Station %in% c('100','200') & Reservoir == "BVR")

# Plot HIX
fcr_res <- 
  ggplot()+
  geom_hline(yintercept=4,linetype="dashed")+
  geom_hline(yintercept=6,linetype="dashed")+
  geom_line(fl_fcr_res,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=1)+
  geom_point(fl_fcr_res,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('20','30','45','50'),
                     values=c("#1B9E77","#CA5902","#7570B3","#1D1A31"))+
  ylim(0,8)+
  xlab("Date")+
  ylab("HIX")+
  labs(fill = "Site")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_res <- ggplot()+
  geom_hline(yintercept=4,linetype="dashed")+
  geom_hline(yintercept=6,linetype="dashed")+
  geom_line(fl_bvr_res,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=1)+
  geom_point(fl_bvr_res,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('20','1','30','45','50'),labels=c('20','30-L','30-R','45','50'),
                     values=c("#1B9E77","#FD8C35","#CA5902","#7570B3","#1D1A31"))+
  ylim(0,8)+
  xlab("Date")+
  ylab("HIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

fcr_in <- 
  ggplot()+
  geom_hline(yintercept=4,linetype="dashed")+
  geom_hline(yintercept=6,linetype="dashed")+
  geom_line(fl_fcr_in,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=1)+
  geom_point(fl_fcr_in,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0,8)+
  xlab("Date")+
  ylab("HIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_in <- ggplot()+
  geom_hline(yintercept=4,linetype="dashed")+
  geom_hline(yintercept=6,linetype="dashed")+
  geom_line(fl_bvr_in,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=1)+
  geom_point(fl_bvr_in,mapping=aes(x=Date,y=HIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0,8)+
  xlab("Date")+
  ylab("HIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(fcr_in,bvr_in,fcr_res,bvr_res,nrow=2,ncol=2,common.legend = FALSE)

ggsave("./Fig_Output/HIX_RCDays.png",width=10,height=7,units="in",dpi=320)

# BIX
fcr_res <- 
  ggplot()+
  geom_hline(yintercept=0.7,linetype="dashed")+
  geom_line(fl_fcr_res,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=1)+
  geom_point(fl_fcr_res,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('20','30','45','50'),
                     values=c("#1B9E77","#CA5902","#7570B3","#1D1A31"))+
  ylim(0.5,0.8)+
  xlab("Date")+
  ylab("BIX")+
  labs(fill = "Site")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_res <- ggplot()+
  geom_hline(yintercept=0.7,linetype="dashed")+
  geom_line(fl_bvr_res,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=1)+
  geom_point(fl_bvr_res,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('20','1','30','45','50'),labels=c('20','30-L','30-R','45','50'),
                     values=c("#1B9E77","#FD8C35","#CA5902","#7570B3","#1D1A31"))+
  ylim(0.5,0.8)+
  xlab("Date")+
  ylab("BIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

fcr_in <- 
  ggplot()+
  geom_hline(yintercept=0.7,linetype="dashed")+
  geom_line(fl_fcr_in,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=1)+
  geom_point(fl_fcr_in,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0.5,0.8)+
  xlab("Date")+
  ylab("BIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

bvr_in <- ggplot()+
  geom_hline(yintercept=0.7,linetype="dashed")+
  geom_line(fl_bvr_in,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=1)+
  geom_point(fl_bvr_in,mapping=aes(x=Date,y=BIX,color=as.factor(Station)),size=2)+
  scale_color_manual(breaks=c('100','200'),
                     values=c("#384D48","#2274A5"))+ 
  ylim(0.5,0.8)+
  xlab("Date")+
  ylab("BIX")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(fcr_in,bvr_in,fcr_res,bvr_res,nrow=2,ncol=2,common.legend = FALSE)

ggsave("./Fig_Output/BIX_RCDays.png",width=10,height=7,units="in",dpi=320)
