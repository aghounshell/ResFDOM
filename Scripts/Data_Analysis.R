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
epi_ghg <- left_join(ar_dates,epi_ghg)

# Included: DOC, temp, DO, Flora, Flow, Rain, SW
epi_dom <- epi %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:ShortwaveRadiationUp_Average_W_m2)
epi_dom <- left_join(ar_dates,epi_dom)

# Included: DOC, temp, DO, Flora, Flow
meta_ghg <- meta %>% select(Date,BIX:HIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
meta_ghg <- left_join(ar_dates,meta_ghg)

# Included: DOC, temp, DO, Flora, Flow
meta_dom <- meta %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
meta_dom <- left_join(ar_dates,meta_dom)

# Included: DOC, temp, DO, Flora, Flow
hypo_ghg <- hypo %>% select(Date,HIX:BIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
hypo_ghg <- left_join(ar_dates,hypo_ghg)

# Included: DOC, temp, DO, Flora, Flow
hypo_dom <- hypo %>% select(Date,Fmax1:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
hypo_dom <- left_join(ar_dates,hypo_dom)

############################################################################################
### Look at frequency of sampling
plot(epi_dom$Date,epi_dom$Fmax1)
plot(epi_ghg$Date,epi_ghg$ch4_umolL)

#################### NEED TO UPDATE FOR NEAR WEEKLY DATA BELOW ###########################

# Extrapolate data to daily then select weekly data from 5-06-19 to 11-20-19
# Create daily timepoints
epi_data <- epi_ghg %>% 
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

# Meta ghg data
meta_data <- meta_ghg %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend"))

# Hypo ghg data
hypo_data <- hypo_ghg %>% 
  mutate(HIX=na.fill(na.approx(HIX),"extend")) %>% 
  mutate(BIX=na.fill(na.approx(BIX),"extend")) %>% 
  mutate(ch4_umolL=na.fill(na.approx(ch4_umolL),"extend")) %>% 
  mutate(co2_umolL=na.fill(na.approx(co2_umolL),"extend")) %>% 
  mutate(DOC_mgL=na.fill(na.approx(DOC_mgL),"extend")) %>% 
  mutate(temp=na.fill(na.approx(temp),"extend")) %>% 
  mutate(DO=na.fill(na.approx(DO),"extend")) %>% 
  mutate(Flora_ugL=na.fill(na.approx(Flora_ugL),"extend")) %>% 
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend"))

###########################################################################
# Check autocorrelation among different parameters

# Methane, epi - 1 lag is important
plot(epi_data$ch4_umolL,type="b")
lag1.plot(epi_data$ch4_umolL,10)
PlotACF(epi_data$ch4_umolL)
acf2(epi_data$ch4_umolL,na.action=na.pass)
pacf(epi_data$ch4_umolL,na.action=na.pass)
xlag1 = lag(epi_data$ch4_umolL,1)
y = cbind(epi_data$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
epi_data <- cbind(epi_data,y)
epi_data <- epi_data %>% select(-ch4)

# CO2, epi - 1 lag is important
plot(epi_data$co2_umolL,type="b")
lag1.plot(epi_data$co2_umolL,10)
PlotACF(epi_data$co2_umolL)
acf2(epi_data$co2_umolL,na.action=na.pass)
pacf(epi_data$co2_umolL,na.action=na.pass)
xlag1 = lag(epi_data$co2_umolL,1)
y = cbind(epi_data$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
epi_data <- cbind(epi_data,y)
epi_data <- epi_data %>% select(-co2)

# HIX, epi - 1 lag important
plot(epi_data$HIX,type="b")
lag1.plot(epi_data$HIX,10)
PlotACF(epi_data$HIX)
acf2(epi_data$HIX,na.action=na.pass)
pacf(epi_data$HIX,na.action=na.pass)
xlag1 = lag(epi_data$HIX,1)
y = cbind(epi_data$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
epi_data <- cbind(epi_data,y)
epi_data <- epi_data %>% select(-hix_ar)

# BIX, epi - 1 lag important
plot(epi_data$BIX,type="b")
lag1.plot(epi_data$BIX,10)
PlotACF(epi_data$BIX)
acf2(epi_data$BIX,na.action=na.pass)
pacf(epi_data$BIX,na.action=na.pass)
xlag1 = lag(epi_data$BIX,1)
y = cbind(epi_data$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
epi_data <- cbind(epi_data,y)
epi_data <- epi_data %>% select(-bix_ar)

# Methane, meta - 1 lag is important
plot(meta_data$ch4_umolL,type="b")
lag1.plot(meta_data$ch4_umolL,10)
PlotACF(meta_data$ch4_umolL)
acf2(meta_data$ch4_umolL,na.action=na.pass)
pacf(meta_data$ch4_umolL,na.action=na.pass)
xlag1 = lag(meta_data$ch4_umolL,1)
y = cbind(meta_data$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
meta_data <- cbind(meta_data,y)
meta_data <- meta_data %>% select(-ch4)

# CO2, meta - 1 lag is important
plot(meta_data$co2_umolL,type="b")
lag1.plot(meta_data$co2_umolL,10)
PlotACF(meta_data$co2_umolL)
acf2(meta_data$co2_umolL,na.action=na.pass)
pacf(meta_data$co2_umolL,na.action=na.pass)
xlag1 = lag(meta_data$co2_umolL,1)
y = cbind(meta_data$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
meta_data <- cbind(meta_data,y)
meta_data <- meta_data %>% select(-co2)

# HIX, meta - NO LAG IMPORTANT???
plot(meta_data$HIX,type="b")
lag1.plot(meta_data$HIX,10)
PlotACF(meta_data$HIX)
acf2(meta_data$HIX,na.action=na.pass)
pacf(meta_data$HIX,na.action=na.pass)
xlag1 = lag(meta_data$HIX,1)
y = cbind(meta_data$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
meta_data <- cbind(meta_data,y)
meta_data <- meta_data %>% select(-hix_ar)

# BIX, meta - 1 lag important
plot(meta_data$BIX,type="b")
lag1.plot(meta_data$BIX,10)
PlotACF(meta_data$BIX)
acf2(meta_data$BIX,na.action=na.pass)
pacf(meta_data$BIX,na.action=na.pass)
xlag1 = lag(meta_data$BIX,1)
y = cbind(meta_data$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
meta_data <- cbind(meta_data,y)
meta_data <- meta_data %>% select(-bix_ar)

# Methane, hypo - NO LAG IS IMPORTANT???
plot(hypo_data$ch4_umolL,type="b")
lag1.plot(hypo_data$ch4_umolL,10)
PlotACF(hypo_data$ch4_umolL)
acf2(hypo_data$ch4_umolL,na.action=na.pass)
pacf(hypo_data$ch4_umolL,na.action=na.pass)
xlag1 = lag(hypo_data$ch4_umolL,1)
y = cbind(hypo_data$ch4_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("ch4","ch4_umolL_ARLag1")
hypo_data <- cbind(hypo_data,y)
hypo_data <- hypo_data %>% select(-ch4)

# CO2, hypo - 1 lag is important
plot(hypo_data$co2_umolL,type="b")
lag1.plot(hypo_data$co2_umolL,10)
PlotACF(hypo_data$co2_umolL)
acf2(hypo_data$co2_umolL,na.action=na.pass)
pacf(hypo_data$co2_umolL,na.action=na.pass)
xlag1 = lag(hypo_data$co2_umolL,1)
y = cbind(hypo_data$co2_umolL,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to methane epi data
colnames(y) <- c("co2","co2_umolL_ARLag1")
hypo_data <- cbind(hypo_data,y)
hypo_data <- hypo_data %>% select(-co2)

# HIX, hypo - NO LAG IMPORTANT???
plot(hypo_data$HIX,type="b")
lag1.plot(hypo_data$HIX,10)
PlotACF(hypo_data$HIX)
acf2(hypo_data$HIX,na.action=na.pass)
pacf(hypo_data$HIX,na.action=na.pass)
xlag1 = lag(hypo_data$HIX,1)
y = cbind(hypo_data$HIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("hix_ar","HIX_ARLag1")
hypo_data <- cbind(hypo_data,y)
hypo_data <- hypo_data %>% select(-hix_ar)

# BIX, hypo - 1 lag important
plot(hypo_data$BIX,type="b")
lag1.plot(hypo_data$BIX,10)
PlotACF(hypo_data$BIX)
acf2(hypo_data$BIX,na.action=na.pass)
pacf(hypo_data$BIX,na.action=na.pass)
xlag1 = lag(hypo_data$BIX,1)
y = cbind(hypo_data$BIX,xlag1)
ar1fit = lm(y[,1]~y[,2])
summary(ar1fit)
plot(ar1fit$fit,ar1fit$residuals)
acf(ar1fit$residuals)
# Add AR lag to HIX epi data
colnames(y) <- c("bix_ar","BIX_ARLag1")
hypo_data <- cbind(hypo_data,y)
hypo_data <- hypo_data %>% select(-bix_ar)

##############################################################
# Plot ACF and PACF for GHG and DOM parameters
pdf("./Fig_Output/ACF_Plots.pdf", width=12, height=8)

acf2(epi_data$ch4_umolL,na.action=na.pass)
acf2(epi_data$co2_umolL,na.action=na.pass)
acf2(epi_data$HIX,na.action=na.pass)
acf2(epi_data$BIX,na.action=na.pass)

acf2(meta_data$ch4_umolL,na.action=na.pass)
acf2(meta_data$co2_umolL,na.action=na.pass)
acf2(meta_data$HIX,na.action=na.pass)
acf2(meta_data$BIX,na.action=na.pass)

acf2(hypo_data$ch4_umolL,na.action=na.pass)
acf2(hypo_data$co2_umolL,na.action=na.pass)
acf2(hypo_data$HIX,na.action=na.pass)
acf2(hypo_data$BIX,na.action=na.pass)

dev.off()

######################################################################################################
# Check for correlations among variables
# Remove variables that are collinear with r2 > 0.60
epi_ghg_corr <- epi_data %>% select(-c(Date))
chart.Correlation(epi_ghg_corr,histogram=TRUE,method=c("spearman"))

# Epi CH4 = ch4_lag + Rain + DO + DOC + Flora
# Remove: BIX_lag; HIX_lag; CO2_lag; SW Radiation; Flow; HIX; BIX; temp; CO2
epi_ch4_corr <- epi_data %>% select(ch4_umolL,ch4_umolL_ARLag1,Rain_Total_mm,DO,DOC_mgL,Flora_ugL)
chart.Correlation(epi_ch4_corr,histogram=TRUE,method=c("spearman"))

epi_data_2 <- epi_data[complete.cases(epi_data),]

model_epi_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + DOC_mgL + DO + Flora_ugL +
                       Rain_Total_mm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_ch4 <- dredge(model_epi_ch4,rank = "AICc")

select_glm_epi_ch4 <- subset(glm_epi_ch4,delta<2)

# Epi CO2 = co2_lag + rain + flora + DO + DOC
# Remove: BIX_lag, HIX_lag, ch4_lag, SW radiation, Flow, CH4, HIX, BIX, temp, Ch4
epi_co2_corr <- epi_data %>% select(co2_umolL,co2_umolL_ARLag1,Rain_Total_mm,DO,DOC_mgL,Flora_ugL)
chart.Correlation(epi_co2_corr,histogram=TRUE,method=c("spearman"))

model_epi_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + DOC_mgL + DO + Flora_ugL +
                       Rain_Total_mm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_co2 <- dredge(model_epi_co2,rank = "AICc")

select_glm_epi_co2 <- subset(glm_epi_co2,delta<2)

# Epi BIX = bix_lag  + rain + Flora + temp + DOC
# Remove: HIX_lag, CO2_lag, CH4_lag, flow, CO2, CH4, HIX, SW_Radiation, DO
epi_bix_corr <- epi_data %>% select(BIX,BIX_ARLag1,Rain_Total_mm,DOC_mgL,Flora_ugL,temp)
chart.Correlation(epi_bix_corr,histogram=TRUE,method=c("spearman"))

model_epi_bix <- glm(BIX ~ BIX_ARLag1 + DOC_mgL + Flora_ugL + temp +
                       Rain_Total_mm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_bix <- dredge(model_epi_bix,rank = "AICc")

select_glm_epi_bix <- subset(glm_epi_bix,delta<2)

# Epi HIX = hix_lag + rain + flora + temp + DOC
# Remove: SW Radiation, DO, BIX, BIX_lag, CO2_lag, CH4_lag, flow, CO2, C4
epi_hix_corr <- epi_data %>% select(HIX,HIX_ARLag1,Rain_Total_mm,DOC_mgL,Flora_ugL,temp)
chart.Correlation(epi_hix_corr,histogram=TRUE,method=c("spearman"))

model_epi_hix <- glm(HIX ~ HIX_ARLag1 + DOC_mgL + Flora_ugL + temp +
                       Rain_Total_mm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_hix <- dredge(model_epi_hix,rank = "AICc")

select_glm_epi_hix <- subset(glm_epi_hix,delta<2)

############ Meta
meta_data_corr <- meta_data %>% select(-Date)
chart.Correlation(meta_data_corr,histogram=TRUE,method=c("spearman"))

meta_data_2 <- meta_data[complete.cases(meta_data),]

# Meta CH4 = ch4_lag + flora + DO + DOC + HIX + BIX
# Remove: flow, temp, co2, co2_lag, bix_lag, hix_lag
meta_ch4_corr <- meta_data %>% select(ch4_umolL,ch4_umolL_ARLag1,Flora_ugL,DO,DOC_mgL,HIX,BIX)
chart.Correlation(meta_ch4_corr,histogram=TRUE,method=c("spearman"))

model_meta_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + Flora_ugL + DO + DOC_mgL + HIX + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_ch4 <- dredge(model_meta_ch4,rank="AICc")

select_glm_meta_ch4 <- subset(glm_meta_ch4,delta<2)

# Meta CO2 = co2_lag + flora + DOC + HIX + BIX
# Remove: hix_lag, bix_lag, ch4_lag, flow, DO, temp, ch4
meta_co2_corr <- meta_data %>% select(co2_umolL,co2_umolL_ARLag1,Flora_ugL,DOC_mgL,HIX,BIX)
chart.Correlation(meta_co2_corr,histogram=TRUE,method=c("spearman"))

model_meta_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + Flora_ugL + DOC_mgL + HIX + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_co2 <- dredge(model_meta_co2,rank="AICc")

select_glm_meta_co2 <- subset(glm_meta_co2,delta<2)

# Meta BIX = bix_lag + flora + temp + DOC + CH4 + HIX
# Remove: co2_lag, ch4_lag, hix_lag, flow, DO, co2
meta_bix_corr <- meta_data %>% select(BIX,BIX_ARLag1,Flora_ugL,temp,DOC_mgL,ch4_umolL,HIX)
chart.Correlation(meta_bix_corr,histogram=TRUE,method=c("spearman"))

model_meta_bix <- glm(BIX ~ BIX_ARLag1 + Flora_ugL + temp + DOC_mgL + ch4_umolL + HIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_bix <- dredge(model_meta_bix,rank="AICc")

select_glm_meta_bix <- subset(glm_meta_bix,delta<2)

# Meta HIX = hix_lag + flora + temp + DOC + CH4 + BIX
# Remove: co2_lag, ch4_lag, bix_lag, flow, DO, CO2
meta_hix_corr <- meta_data %>% select(HIX,HIX_ARLag1,Flora_ugL,temp,DOC_mgL,ch4_umolL,BIX)
chart.Correlation(meta_hix_corr,histogram=TRUE,method=c("spearman"))

model_meta_hix <- glm(HIX ~ HIX_ARLag1 + Flora_ugL + temp + DOC_mgL + ch4_umolL + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_hix <- dredge(model_meta_hix,rank="AICc")

select_glm_meta_hix <- subset(glm_meta_hix,delta<2)

########################## Hypo
hypo_data_corr <- hypo_data %>% select(-Date)
chart.Correlation(hypo_data_corr,histogram=TRUE,method=c("spearman"))

hypo_data_2 <- hypo_data[complete.cases(hypo_data),]

# Hypo CH4 = ch4_lag + DO + DOC + BIX + HIX
# Remove: co2_lag, hix_lag, bix_lag, flow, flora, co2, temp
hypo_ch4_corr <- hypo_data %>% select(ch4_umolL,ch4_umolL_ARLag1,DO,DOC_mgL,BIX,HIX)
chart.Correlation(hypo_ch4_corr,histogram=TRUE,method=c("spearman"))

model_hypo_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + BIX + HIX + DOC_mgL + DO, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_ch4 <- dredge(model_hypo_ch4,rank="AICc")

select_glm_hypo_ch4 <- subset(glm_hypo_ch4,delta<2)

# Hypo Co2 = co2_lag + DO + DOC + HIX
# Remove: co2_lag, hix_lag, bix_lag, temp, BIX, flora
hypo_co2_corr <- hypo_data %>% select(co2_umolL,co2_umolL_ARLag1,DO,DOC_mgL,HIX)
chart.Correlation(hypo_co2_corr,histogram=TRUE,method=c("spearman"))

model_hypo_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + HIX + DOC_mgL + DO, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_co2 <- dredge(model_hypo_co2,rank="AICc")

select_glm_hypo_co2 <- subset(glm_hypo_co2,delta<2)

# Hypo BIX = bix_lag + DO + DOC + CO2 + HIX
# Remove: ch4_lag, co2_lag, hix_lag, temp, flow, flora, ch4
hypo_bix_corr <- hypo_data %>% select(BIX,BIX_ARLag1,DO,DOC_mgL,co2_umolL,HIX)
chart.Correlation(hypo_bix_corr,histogram=TRUE,method=c("spearman"))

model_hypo_bix <- glm(BIX ~ BIX_ARLag1 + HIX + DOC_mgL + DO + co2_umolL, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_bix <- dredge(model_hypo_bix,rank="AICc")

select_glm_hypo_bix <- subset(glm_hypo_bix,delta<2)

# Hypo HIX = hix_lag + DO + temp + DOC
# Remove: ch4_lag, co2_lag, bix_lag, flow, flora, co2, ch4, BIX
hypo_hix_corr <- hypo_data %>% select(HIX,HIX_ARLag1,DO,DOC_mgL,temp)
chart.Correlation(hypo_hix_corr,histogram=TRUE,method=c("spearman"))

model_hypo_hix <- glm(HIX ~ HIX_ARLag1 + DOC_mgL + DO + temp, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_hix <- dredge(model_hypo_hix,rank="AICc")

select_glm_hypo_hix <- subset(glm_hypo_hix,delta<2)