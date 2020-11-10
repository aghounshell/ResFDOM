### Script to conduct data analyses on all compiled data (see All_Data.R)
### A Hounshell, 04 Nov 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn)

# Load compiled data
data <- read_csv("./Data/20201103_All_Data.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%Y-%m-%d", tz = "EST"))

# Epi data
epi <- data %>% filter(Station == 50 & Depth == 0.1)

# Meta data
meta <- data %>% filter(Station == 50 & Depth == 5)

# Hypo data
hypo <- data %>% filter(Station == 50 & Depth == 9)

# Inflow (weir) data
weir <- data %>% filter(Station == 100)

# Missing some inflow data
inflow <- read_csv('./Data/Inflow.csv') %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d",tz='EST'))) %>% 
  filter(DateTime > as.POSIXct("2019-01-01") & DateTime < as.POSIXct("2020-01-01")) %>% 
  select(Site,DateTime,WVWA_Flow_cms) %>% 
  rename(Station = Site, Date = DateTime, Flow_cms = WVWA_Flow_cms)

inflow <- inflow %>%  group_by(Date) %>% summarise_all(funs(mean(.,na.rm=TRUE)))

# Combine w/ Epi, Meta, and Hypo data
epi <- left_join(epi,inflow,by="Date")
epi <- epi %>% 
  mutate(Flow_cms.x = Flow_cms.y) %>% 
  select(-c("Station.y","Flow_cms.y")) %>% 
  rename(Flow_cms = Flow_cms.x)

meta <- left_join(meta,inflow,by="Date")
meta <- meta %>% 
  mutate(Flow_cms.x = Flow_cms.y) %>% 
  select(-c("Station.y","Flow_cms.y")) %>% 
  rename(Flow_cms = Flow_cms.x)

hypo <- left_join(hypo,inflow,by="Date")
hypo <- hypo %>% 
  mutate(Flow_cms.x = Flow_cms.y) %>% 
  select(-c("Station.y","Flow_cms.y")) %>% 
  rename(Flow_cms = Flow_cms.x)

## Separate into GHG and DOM datasets and filter for complete cases
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Load in 'weekly' dates for AR model
ar_dates <- read.csv("./Data/AR_Dates.csv") %>% 
  mutate(AR.Model.Dates = as.POSIXct(strptime(AR.Model.Dates, "%m/%d/%Y",tz='EST'))) %>% 
  rename(Date = AR.Model.Dates)

# Included: DOC, temp, DO, Flora, Flow, Rain, SW
epi_ghg <- epi %>% select(Date,BIX:HIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:ShortwaveRadiationUp_Average_W_m2)
epi_ghg <- completeFun(epi_ghg,c("ch4_umolL","co2_umolL"))
epi_ghg <- epi_ghg %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))
epi_ghg <- left_join(ar_dates,epi_ghg)

#################### NEED TO UPDATE FOR NEAR WEEKLY DATA BELOW ###########################

# Included: DOC, temp, DO, Flora, Flow, Rain, SW
epi_dom <- epi %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:ShortwaveRadiationUp_Average_W_m2)
epi_dom <- completeFun(epi_dom,"HIX")
epi_dom <- epi_dom %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))

# Included: DOC, temp, DO, Flora, Flow
meta_ghg <- meta %>% select(Date,BIX:HIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
meta_ghg <- completeFun(meta_ghg,c("ch4_umolL","co2_umolL"))
meta_ghg <- meta_ghg %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))

# Included: DOC, temp, DO, Flora, Flow
meta_dom <- meta %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
meta_dom <- completeFun(meta_dom,"HIX")
meta_dom <- meta_dom %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))

# Included: DOC, temp, DO, Flora, Flow
hypo_ghg <- hypo %>% select(Date,HIX:BIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
hypo_ghg <- completeFun(hypo_ghg,c("ch4_umolL","co2_umolL"))
hypo_ghg <- hypo_ghg %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))

# Included: DOC, temp, DO, Flora, Flow
hypo_dom <- hypo %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
hypo_dom <- completeFun(hypo_dom,"HIX")
hypo_dom <- hypo_dom %>% filter(Date > as.POSIXct("2019-05-19") & Date < as.POSIXct("2019-11-21"))

############################################################################################
### Look at frequency of sampling
plot(epi_dom$Date,epi_dom$Fmax1)
plot(epi_ghg$Date,epi_ghg$ch4_umolL)

# Extrapolate data to daily then select weekly data from 5-06-19 to 11-20-19
# Create daily timepoints
epi_ghg_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
epi_ghg_ts <- full_join(epi_ghg_ts,epi_ghg)
epi_ghg_ts <- epi_ghg_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend")) %>% 
  mutate(Rain_Total_mm=na.fill(na.approx(Rain_Total_mm),"extend")) %>% 
  mutate(ShortwaveRadiationUp_Average_W_m2=na.fill(na.approx(ShortwaveRadiationUp_Average_W_m2),"extend")) 
  
epi_ghg_weekly <- epi_ghg_ts[seq(1,nrow(epi_ghg_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(epi_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_point(epi_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_line(epi_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  geom_point(epi_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  theme_classic(base_size=15)

# Meta ghg data
meta_ghg_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
meta_ghg_ts <- full_join(meta_ghg_ts,meta_ghg)
meta_ghg_ts <- meta_ghg_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend"))

meta_ghg_weekly <- meta_ghg_ts[seq(1,nrow(meta_ghg_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(meta_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_point(meta_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_line(meta_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  geom_point(meta_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  theme_classic(base_size=15)

# Hypo ghg data
hypo_ghg_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
hypo_ghg_ts <- full_join(hypo_ghg_ts,hypo_ghg)
hypo_ghg_ts <- hypo_ghg_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend"))

hypo_ghg_weekly <- hypo_ghg_ts[seq(1,nrow(hypo_ghg_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(hypo_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_point(hypo_ghg_weekly,mapping=aes(x=Date,y=ch4_umolL,color="Weekly"))+
  geom_line(hypo_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  geom_point(hypo_ghg,mapping=aes(x=Date,y=ch4_umolL,color="Sampling"))+
  theme_classic(base_size=15)

###########################################################################
# Check autocorrelation among different parameters

# Methane, epi - 1 lag is important
plot(epi_ghg_weekly$ch4_umolL,type="b")
lag1.plot(epi_ghg_weekly$ch4_umolL,10)
PlotACF(epi_ghg_weekly$ch4_umolL)
acf2(epi_ghg_weekly$ch4_umolL,na.action=na.pass)
pacf(epi_ghg_weekly$ch4_umolL,na.action=na.pass)
xlag1 = lag(epi_ghg_weekly$ch4_umolL,1)
y = cbind(epi_ghg_weekly$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
epi_ghg_weekly <- cbind(epi_ghg_weekly,y)
epi_ghg_weekly <- epi_ghg_weekly %>% select(-ch4)

# CO2, epi - 1 lag is important
plot(epi_ghg_weekly$co2_umolL,type="b")
lag1.plot(epi_ghg_weekly$co2_umolL,10)
PlotACF(epi_ghg_weekly$co2_umolL)
acf2(epi_ghg_weekly$co2_umolL,na.action=na.pass)
pacf(epi_ghg_weekly$co2_umolL,na.action=na.pass)
xlag1 = lag(epi_ghg_weekly$co2_umolL,1)
y = cbind(epi_ghg_weekly$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
epi_ghg_weekly <- cbind(epi_ghg_weekly,y)
epi_ghg_weekly <- epi_ghg_weekly %>% select(-co2)

# Methane, meta - 1 lag is important
plot(meta_ghg_weekly$ch4_umolL,type="b")
lag1.plot(meta_ghg_weekly$ch4_umolL,10)
PlotACF(meta_ghg_weekly$ch4_umolL)
acf2(meta_ghg_weekly$ch4_umolL,na.action=na.pass)
pacf(meta_ghg_weekly$ch4_umolL,na.action=na.pass)
xlag1 = lag(meta_ghg_weekly$ch4_umolL,1)
y = cbind(meta_ghg_weekly$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
meta_ghg_weekly <- cbind(meta_ghg_weekly,y)
meta_ghg_weekly <- meta_ghg_weekly %>% select(-ch4)

# CO2, meta - 1 lag is important
plot(meta_ghg_weekly$co2_umolL,type="b")
lag1.plot(meta_ghg_weekly$co2_umolL,10)
PlotACF(meta_ghg_weekly$co2_umolL)
acf2(meta_ghg_weekly$co2_umolL,na.action=na.pass)
pacf(meta_ghg_weekly$co2_umolL,na.action=na.pass)
xlag1 = lag(meta_ghg_weekly$co2_umolL,1)
y = cbind(meta_ghg_weekly$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
meta_ghg_weekly <- cbind(meta_ghg_weekly,y)
meta_ghg_weekly <- meta_ghg_weekly %>% select(-co2)

# Methane, hypo - 1 lag is important
plot(hypo_ghg_weekly$ch4_umolL,type="b")
lag1.plot(hypo_ghg_weekly$ch4_umolL,10)
PlotACF(hypo_ghg_weekly$ch4_umolL)
acf2(hypo_ghg_weekly$ch4_umolL,na.action=na.pass)
pacf(hypo_ghg_weekly$ch4_umolL,na.action=na.pass)
xlag1 = lag(hypo_ghg_weekly$ch4_umolL,1)
y = cbind(hypo_ghg_weekly$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
hypo_ghg_weekly <- cbind(hypo_ghg_weekly,y)
hypo_ghg_weekly <- hypo_ghg_weekly %>% select(-ch4)

# CO2, hypo - 1 lag is important
plot(hypo_ghg_weekly$co2_umolL,type="b")
lag1.plot(hypo_ghg_weekly$co2_umolL,10)
PlotACF(hypo_ghg_weekly$co2_umolL)
acf2(hypo_ghg_weekly$co2_umolL,na.action=na.pass)
pacf(hypo_ghg_weekly$co2_umolL,na.action=na.pass)
xlag1 = lag(hypo_ghg_weekly$co2_umolL,1)
y = cbind(hypo_ghg_weekly$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
hypo_ghg_weekly <- cbind(hypo_ghg_weekly,y)
hypo_ghg_weekly <- hypo_ghg_weekly %>% select(-co2)

######################################################################################################
# Check for correlations among variables
epi_ghg_corr <- epi_ghg_weekly %>% select(-Date)
chart.Correlation(epi_ghg_corr,histogram=TRUE,method=c("spearman"))

meta_ghg_corr <- meta_ghg_weekly %>% select(-Date)
chart.Correlation(meta_ghg_corr,histogram=TRUE,method=c("spearman"))

hypo_ghg_corr <- hypo_ghg_weekly %>% select(-Date)
chart.Correlation(hypo_ghg_corr,histogram=TRUE,method=c("spearman"))

####################################################################################################
# Develop AR Models for GHG data
# Build a global model

############### Epi
epi_ghg_weekly_2 <- epi_ghg_weekly[complete.cases(epi_ghg_weekly),]

# Remove: HIX, CO2
model_epi_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + BIX + DOC_mgL + temp + DO + Flora_ugL + Flow_cms + 
                       Rain_Total_mm + ShortwaveRadiationUp_Average_W_m2, data = epi_ghg_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_ch4 <- dredge(model_epi_ch4,rank = "AICc", fixed = "ch4_umolL_ARLag1")

select_glm_epi_ch4 <- subset(glm_epi_ch4,delta<2)

# Remove: HIX, CH4
model_epi_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + BIX + DOC_mgL + temp + DO + Flora_ugL + Flow_cms + 
                       Rain_Total_mm + ShortwaveRadiationUp_Average_W_m2, data = epi_ghg_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_co2 <- dredge(model_epi_co2,rank="AICc",fixed="co2_umolL_ARLag1")

select_glm_epi_co2 <- subset(glm_epi_co2,delta<2)

############ Meta
meta_ghg_weekly_2 <- meta_ghg_weekly[complete.cases(meta_ghg_weekly),]

# Remove: CO2
model_meta_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + BIX + HIX + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, data = meta_ghg_weekly_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_ch4 <- dredge(model_meta_ch4,rank="AICc",fixed="ch4_umolL_ARLag1")

select_glm_meta_ch4 <- subset(glm_meta_ch4,delta<2)

# Keep all parameters
model_meta_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + BIX + HIX + DOC_mgL + ch4_umolL + temp + DO + Flora_ugL + Flow_cms, 
                      data = meta_ghg_weekly_2, family = gaussian, na.action = 'na.fail')

glm_meta_co2 <- dredge(model_meta_co2,rank="AICc",fixed="co2_umolL_ARLag1")

select_glm_meta_co2 <- subset(glm_meta_co2,delta<2)

############### Hypo
hypo_ghg_weekly_2 <- hypo_ghg_weekly[complete.cases(hypo_ghg_weekly),]

# Remove: CO2
model_hypo_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + BIX + HIX + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, 
                      data = hypo_ghg_weekly_2, family = gaussian, na.action = 'na.fail')

glm_hypo_ch4 <- dredge(model_hypo_ch4,rank="AICc",fixed="ch4_umolL_ARLag1")

select_glm_hypo_ch4 <- subset(glm_hypo_ch4,delta<2)

# Keep all parameters (r2 < 0.8)
model_hypo_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + BIX + HIX + DOC_mgL + ch4_umolL + temp + DO + Flora_ugL + Flow_cms, 
                      data = hypo_ghg_weekly_2, family = gaussian, na.action = 'na.fail')

glm_hypo_co2 <- dredge(model_hypo_co2,rank="AICc",fixed="co2_umolL_ARLag1")

select_glm_hypo_co2 <- subset(glm_hypo_co2,delta<2)

############ NOW DO IT ALL FOR DOM DATA - USE HIX AND BIX AS RESPONSE VARIABLES ###########
### Extrapolate data to daily
# Extrapolate data to daily then select weekly data from 5-20-19 to 11-20-19
# Create daily timepoints
epi_dom_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
epi_dom_ts <- full_join(epi_dom_ts,epi_dom)
epi_dom_ts <- epi_dom_ts %>% select(-c(Fmax1:B,'A/T':'M/C'))
epi_dom_ts <- epi_dom_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend")) %>% 
  mutate(Rain_Total_mm=na.fill(na.approx(Rain_Total_mm),"extend")) %>% 
  mutate(ShortwaveRadiationUp_Average_W_m2=na.fill(na.approx(ShortwaveRadiationUp_Average_W_m2),"extend")) 

epi_dom_weekly <- epi_dom_ts[seq(1,nrow(epi_dom_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(epi_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_point(epi_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_line(epi_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  geom_point(epi_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  theme_classic(base_size=15)

# Meta
meta_dom_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
meta_dom_ts <- full_join(meta_dom_ts,meta_dom)
meta_dom_ts <- meta_dom_ts %>% select(-c(Fmax1:B,'A/T':'M/C'))
meta_dom_ts <- meta_dom_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend")) 

meta_dom_weekly <- meta_dom_ts[seq(1,nrow(meta_dom_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(meta_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_point(meta_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_line(meta_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  geom_point(meta_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  theme_classic(base_size=15)

# Hypo
hypo_dom_ts <- inflow %>% 
  filter(Date > as.POSIXct("2019-05-20") & Date < as.POSIXct("2019-11-21")) %>% 
  select(Date)
hypo_dom_ts <- full_join(hypo_dom_ts,hypo_dom)
hypo_dom_ts <- hypo_dom_ts %>% select(-c(Fmax1:B,'A/T':'M/C'))
hypo_dom_ts <- hypo_dom_ts %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend")) 

hypo_dom_weekly <- hypo_dom_ts[seq(1,nrow(hypo_dom_ts),7), ]

# Check extrapolation
ggplot()+
  geom_line(hypo_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_point(hypo_dom_weekly,mapping=aes(x=Date,y=BIX,color="Weekly"))+
  geom_line(hypo_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  geom_point(hypo_dom,mapping=aes(x=Date,y=BIX,color="Sampling"))+
  theme_classic(base_size=15)

###########################################################################
# Check autocorrelation among different parameters

# HIX, epi - 1 lag important
plot(epi_dom_weekly$HIX,type="b")
lag1.plot(epi_dom_weekly$HIX,10)
PlotACF(epi_dom_weekly$HIX)
acf2(epi_dom_weekly$HIX,na.action=na.pass)
pacf(epi_dom_weekly$HIX,na.action=na.pass)
xlag1 = lag(epi_dom_weekly$HIX,1)
y = cbind(epi_dom_weekly$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
epi_dom_weekly <- cbind(epi_dom_weekly,y)
epi_dom_weekly <- epi_dom_weekly %>% select(-hix_ar)

# BIX, epi - 1 lag important
plot(epi_dom_weekly$BIX,type="b")
lag1.plot(epi_dom_weekly$BIX,10)
PlotACF(epi_dom_weekly$BIX)
acf2(epi_dom_weekly$BIX,na.action=na.pass)
pacf(epi_dom_weekly$BIX,na.action=na.pass)
xlag1 = lag(epi_dom_weekly$BIX,1)
y = cbind(epi_dom_weekly$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
epi_dom_weekly <- cbind(epi_dom_weekly,y)
epi_dom_weekly <- epi_dom_weekly %>% select(-bix_ar)

# HIX, meta - NO LAG IMPORTANT???
plot(meta_dom_weekly$HIX,type="b")
lag1.plot(meta_dom_weekly$HIX,10)
PlotACF(meta_dom_weekly$HIX)
acf2(meta_dom_weekly$HIX,na.action=na.pass)
pacf(meta_dom_weekly$HIX,na.action=na.pass)
xlag1 = lag(meta_dom_weekly$HIX,1)
y = cbind(meta_dom_weekly$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
meta_dom_weekly <- cbind(meta_dom_weekly,y)
meta_dom_weekly <- meta_dom_weekly %>% select(-hix_ar)

# BIX, meta - 1 lag important
plot(meta_dom_weekly$BIX,type="b")
lag1.plot(meta_dom_weekly$BIX,10)
PlotACF(meta_dom_weekly$BIX)
acf2(meta_dom_weekly$BIX,na.action=na.pass)
pacf(meta_dom_weekly$BIX,na.action=na.pass)
xlag1 = lag(meta_dom_weekly$BIX,1)
y = cbind(meta_dom_weekly$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
meta_dom_weekly <- cbind(meta_dom_weekly,y)
meta_dom_weekly <- meta_dom_weekly %>% select(-bix_ar)

# HIX, hypo - NO LAG IMPORTANT???
plot(hypo_dom_weekly$HIX,type="b")
lag1.plot(hypo_dom_weekly$HIX,10)
PlotACF(hypo_dom_weekly$HIX)
acf2(hypo_dom_weekly$HIX,na.action=na.pass)
pacf(hypo_dom_weekly$HIX,na.action=na.pass)
xlag1 = lag(hypo_dom_weekly$HIX,1)
y = cbind(hypo_dom_weekly$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
hypo_dom_weekly <- cbind(hypo_dom_weekly,y)
hypo_dom_weekly <- hypo_dom_weekly %>% select(-hix_ar)

# BIX, hypo - 1 lag important
plot(hypo_dom_weekly$BIX,type="b")
lag1.plot(hypo_dom_weekly$BIX,10)
PlotACF(hypo_dom_weekly$BIX)
acf2(hypo_dom_weekly$BIX,na.action=na.pass)
pacf(hypo_dom_weekly$BIX,na.action=na.pass)
xlag1 = lag(hypo_dom_weekly$BIX,1)
y = cbind(hypo_dom_weekly$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
hypo_dom_weekly <- cbind(hypo_dom_weekly,y)
hypo_dom_weekly <- hypo_dom_weekly %>% select(-bix_ar)

######################################################################################################
# Check for correlations among variables
epi_dom_corr <- epi_dom_weekly %>% select(-Date)
chart.Correlation(epi_dom_corr,histogram=TRUE,method=c("spearman"))

meta_dom_corr <- meta_dom_weekly %>% select(-Date)
chart.Correlation(meta_dom_corr,histogram=TRUE,method=c("spearman"))

hypo_dom_corr <- hypo_dom_weekly %>% select(-Date)
chart.Correlation(hypo_dom_corr,histogram=TRUE,method=c("spearman"))

####################################################################################################
# Develop AR Models for GHG data
# Build a global model

############### Epi
epi_dom_weekly_2 <- epi_dom_weekly[complete.cases(epi_dom_weekly),]

# Remove CO2, CH4
model_epi_BIX <- glm(BIX ~ BIX_ARLag1 + HIX + DOC_mgL + temp + DO + Flora_ugL + Flow_cms + 
                       Rain_Total_mm + ShortwaveRadiationUp_Average_W_m2, data = epi_dom_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_BIX <- dredge(model_epi_BIX,rank = "AICc", fixed = "BIX_ARLag1")

select_glm_epi_BIX <- subset(glm_epi_BIX,delta<2)

# Remove BIX (correlated w/ SW), CO2, CH4
model_epi_HIX <- glm(HIX ~ HIX_ARLag1 + DOC_mgL + temp + DO + Flora_ugL + Flow_cms + 
                       Rain_Total_mm + ShortwaveRadiationUp_Average_W_m2, data = epi_dom_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_HIX <- dredge(model_epi_HIX,rank = "AICc", fixed = "HIX_ARLag1")

select_glm_epi_HIX <- subset(glm_epi_HIX,delta<2)

######################## Meta
meta_dom_weekly_2 <- meta_dom_weekly[complete.cases(meta_dom_weekly),]

# Remove CO2
model_meta_BIX <- glm(BIX ~ BIX_ARLag1 + HIX + ch4_umolL + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, data = meta_dom_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_meta_BIX <- dredge(model_meta_BIX,rank = "AICc", fixed = "BIX_ARLag1")

select_glm_meta_BIX <- subset(glm_meta_BIX,delta<2)

# Remove CO2
model_meta_HIX <- glm(HIX ~ BIX + ch4_umolL + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, data = meta_dom_weekly_2, 
                     family = gaussian, na.action = 'na.fail')

glm_meta_HIX <- dredge(model_meta_HIX,rank = "AICc")

select_glm_meta_HIX <- subset(glm_meta_HIX,delta<2)

########################## Hypo
hypo_dom_weekly_2 <- hypo_dom_weekly[complete.cases(hypo_dom_weekly),]

# Remove CO2
model_hypo_BIX <- glm(BIX ~ BIX_ARLag1 + HIX + ch4_umolL + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, 
                      data = hypo_dom_weekly_2, family = gaussian, na.action = 'na.fail')

glm_hypo_BIX <- dredge(model_hypo_BIX,rank = "AICc", fixed = "BIX_ARLag1")

select_glm_hypo_BIX <- subset(glm_hypo_BIX,delta<2)

# Remove CO2
model_hypo_HIX <- glm(HIX ~ BIX + ch4_umolL + DOC_mgL + temp + DO + Flora_ugL + Flow_cms, data = hypo_dom_weekly_2, 
                      family = gaussian, na.action = 'na.fail')

glm_hypo_HIX <- dredge(model_hypo_HIX,rank = "AICc")

select_glm_hypo_HIX <- subset(glm_hypo_HIX,delta<2)

##############################################################
# Plot ACF and PACF for GHG and DOM parameters
pdf("./Fig_Output/ACF_Plots.pdf", width=12, height=8)

acf2(epi_ghg_weekly$ch4_umolL,na.action=na.pass)
acf2(epi_ghg_weekly$co2_umolL,na.action=na.pass)
acf2(epi_dom_weekly$HIX,na.action=na.pass)
acf2(epi_dom_weekly$BIX,na.action=na.pass)

acf2(meta_ghg_weekly$ch4_umolL,na.action=na.pass)
acf2(meta_ghg_weekly$co2_umolL,na.action=na.pass)
acf2(meta_dom_weekly$HIX,na.action=na.pass)
acf2(meta_dom_weekly$BIX,na.action=na.pass)

acf2(hypo_ghg_weekly$ch4_umolL,na.action=na.pass)
acf2(hypo_ghg_weekly$co2_umolL,na.action=na.pass)
acf2(hypo_dom_weekly$HIX,na.action=na.pass)
acf2(hypo_dom_weekly$BIX,na.action=na.pass)

dev.off()
