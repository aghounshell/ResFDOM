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

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

doc_fcr <- completeFun(doc_fcr,"DOC_mgL")

doc_bvr <- completeFun(doc_bvr,"DOC_mgL")

doc_fcr_res <- doc_fcr %>% 
  filter(Site %in% c('20','30','45','50'))

doc_fcr_in <- doc_fcr %>% 
  filter(Site %in% c('100','200'))

doc_bvr_res <- doc_bvr %>% 
  filter(Site %in% c('20','30','45','50','1'))

doc_bvr_in <- doc_bvr %>% 
  filter(Site %in% c('100','200'))

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

#####

# Fluorescence (HIX, BIX) ----
fl <- read_csv("./Data/20210210_ResultsFiles_ResEEMs2019_RAW.csv") %>% 
  filter(Dilution %in% c(1,2)) %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%m/%d/%Y", tz="EST"))) %>% 
  filter(Date >= as.POSIXct("2019-04-29") & Date <= as.POSIXct("2020-04-01")) %>% 
  filter(Station %in% c('100','200','20','30','45','50','1') & Depth == 0.1) %>% 
  mutate(Loc = paste(Reservoir,Station))

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
