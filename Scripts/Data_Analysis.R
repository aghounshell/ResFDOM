### Script to conduct data analyses on all compiled data (see All_Data.R)
### A Hounshell, 04 Nov 2020

setwd("C:/Users/ahoun/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics)

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

###################################################################################
## Separate into GHG and DOM datasets and filter for complete cases
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Load in 'weekly' dates for AR model
ar_dates <- read.csv("./Data/AR_Dates.csv") %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%m/%d/%Y",tz='EST')))

# Included: DOC, temp, DO, Flora, Flow, Rain, SW
epi_ghg <- epi %>% select(Date,BIX:HIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
epi_ghg <- left_join(ar_dates,epi_ghg)

# Included: DOC, temp, DO, Flora, Flow
meta_ghg <- meta %>% select(Date,BIX:HIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
meta_ghg <- left_join(ar_dates,meta_ghg)

# Included: DOC, temp, DO, Flora, Flow
hypo_ghg <- hypo %>% select(Date,HIX:BIX,ch4_umolL:co2_umolL,DOC_mgL,temp:DO,Flora_ugL:Flow_cms)
hypo_ghg <- left_join(ar_dates,hypo_ghg)

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
  mutate(Flow_cms=na.fill(na.approx(Flow_cms),"extend"))

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
####################### For AGU - Focus on GHG variables!!!! ##############################
# Check for correlations among variables
# Remove variables that are collinear with r2 > 0.60
epi_data_2 <- epi_data[complete.cases(epi_data),]
epi_data_2 <- epi_data_2 %>% select(-c(Date,HIX))
epi_data_2 <- as.data.frame(epi_data_2)
epi_data_2 <- epi_data_2 %>% 
  mutate(Flow_cms_norm = log(Flow_cms)) %>%
  mutate(co2_umolL_norm = log(co2_umolL)) %>% 
  mutate(co2_umolL_ARLag1_norm = log(co2_umolL_ARLag1)) %>% 
  select(-c(co2_umolL,Flow_cms,co2_umolL_ARLag1))
epi_data_3 <- scale(epi_data_2)
epi_data_3 <- as.data.frame(epi_data_3)
epi_data_3 <- epi_data_3 %>% rename(ch4_umolL_zt = ch4_umolL,co2_umolL_norm_zt = co2_umolL_norm)
epi_data_2 <- cbind(epi_data_2$ch4_umolL,epi_data_2$co2_umolL_norm,epi_data_3)
epi_data_2 <- epi_data_2 %>% rename(ch4_umolL = `epi_data_2$ch4_umolL`,co2_umolL_norm = `epi_data_2$co2_umolL_norm`)

chart.Correlation(epi_data_2,histogram=TRUE,method=c("spearman"))
epi_data_2 <- as.data.frame(epi_data_2)

# Epi CH4 = ch4_lag + flow + DOC + flora + BIX
# Remove: co2_lag, temp, co2, do
epi_ch4_corr <- epi_data_2 %>% select(ch4_umolL,ch4_umolL_ARLag1,Flow_cms_norm,DOC_mgL,Flora_ugL,BIX)
chart.Correlation(epi_ch4_corr,histogram=TRUE,method=c("spearman"))

model_epi_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + DOC_mgL + BIX + Flora_ugL +
                       Flow_cms_norm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_ch4 <- dredge(model_epi_ch4,rank = "AICc")

select_glm_epi_ch4 <- subset(glm_epi_ch4,delta<2)

# Epi CO2 = co2_lag + flow + DO + DOC + flora
# Remove: ch4_lag, temp, ch4, bix
epi_co2_corr <- epi_data_2 %>% select(co2_umolL_norm,co2_umolL_ARLag1_norm,Flow_cms_norm,DO,DOC_mgL,Flora_ugL)
chart.Correlation(epi_co2_corr,histogram=TRUE,method=c("spearman"))

model_epi_co2 <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm + DOC_mgL + DO + Flora_ugL +
                       Flow_cms_norm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_co2 <- dredge(model_epi_co2,rank = "AICc")

select_glm_epi_co2 <- subset(glm_epi_co2,delta<2)

########################## Hypo
hypo_data_2 <- hypo_data[complete.cases(hypo_data),]
hypo_data_2 <- hypo_data_2 %>% select(-c(Date,HIX))
hypo_data_2 <- as.data.frame(hypo_data_2)
hypo_data_2 <- hypo_data_2 %>% 
  mutate(Flow_cms_norm = log(Flow_cms)) %>% 
  mutate(Flora_ugL_norm = log(Flora_ugL)) %>% 
  mutate(ch4_umolL_norm = log(ch4_umolL)) %>% 
  mutate(ch4_umolL_ARLag1_norm = log(ch4_umolL_ARLag1)) %>% 
  select(-c(Flow_cms,ch4_umolL,ch4_umolL_ARLag1,Flora_ugL))
hypo_data_2 <- scale(hypo_data_2)
hypo_data_3 <- scale(hypo_data_2)
hypo_data_3 <- as.data.frame(hypo_data_3)
hypo_data_3 <- hypo_data_3 %>% rename(ch4_umolL_norm_zt = ch4_umolL_norm,co2_umolL_zt = co2_umolL)
hypo_data_2 <- as.data.frame(hypo_data_2)
hypo_data_2 <- cbind(hypo_data_2$ch4_umolL_norm,hypo_data_2$co2_umolL,hypo_data_3)
hypo_data_2 <- hypo_data_2 %>% rename(ch4_umolL_norm = `hypo_data_2$ch4_umolL_norm`,co2_umolL = `hypo_data_2$co2_umolL`)


chart.Correlation(hypo_data_2,histogram=TRUE,method=c("spearman"))
hypo_data_2 <- as.data.frame(hypo_data_2)

# Hypo CH4 = ch4_lag + DO + Flow + DOC + CO2 + BIX
# Remove: co2_lag, temp, flora
hypo_ch4_corr <- hypo_data_2 %>% select(ch4_umolL_norm,ch4_umolL_ARLag1_norm,Flow_cms_norm,DO,DOC_mgL,co2_umolL_zt,BIX)
chart.Correlation(hypo_ch4_corr,histogram=TRUE,method=c("spearman"))

model_hypo_ch4 <- glm(ch4_umolL_norm ~ ch4_umolL_ARLag1_norm + BIX + DOC_mgL + DO + Flow_cms_norm + co2_umolL_zt, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_ch4 <- dredge(model_hypo_ch4,rank="AICc")

select_glm_hypo_ch4 <- subset(glm_hypo_ch4,delta<2)

# Hypo Co2 = co2_lag + Flow + DO + DOC
# Remove: ch4_lag, temp, BIX, flora,ch4
hypo_co2_corr <- hypo_data_2 %>% select(co2_umolL,co2_umolL_ARLag1,Flow_cms_norm,DO,DOC_mgL)
chart.Correlation(hypo_co2_corr,histogram=TRUE,method=c("spearman"))

model_hypo_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + DOC_mgL + DO + Flow_cms_norm, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_co2 <- dredge(model_hypo_co2,rank="AICc")

select_glm_hypo_co2 <- subset(glm_hypo_co2,delta<2)

############################## 
# Calculate R2 for each developed model
mod1_epi_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1, data = epi_data_2, 
                    family = gaussian, na.action = 'na.fail')

mod1_epi_co2_lag <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm, data = epi_data_2,
                        family = gaussian, na.action = 'na.fail')
glm_epi_co2_lag <- dredge(mod1_epi_co2_lag)

mod1_epi_co2 <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm + DO + Flora_ugL +
                      Flow_cms_norm, data = epi_data_2, 
                    family = gaussian, na.action = 'na.fail')
mod2_epi_co2 <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm + DO, data = epi_data_2)
mod3_epi_co2 <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm + DO + Flora_ugL, 
                    data = epi_data_2)
mod4_epi_co2 <- glm(co2_umolL_norm ~ co2_umolL_ARLag1_norm + DO + Flow_cms_norm,
                    data = epi_data_2)

mod1_hypo_ch4_lag <- glm(ch4_umolL_norm ~ ch4_umolL_ARLag1_norm, 
                     data = hypo_data_2, family = gaussian, na.action = 'na.fail')
glm_hypo_ch4_lag <- dredge(mod1_hypo_ch4_lag)

mod1_hypo_ch4 <- glm(ch4_umolL_norm ~ BIX + DOC_mgL + DO, 
                     data = hypo_data_2, family = gaussian, na.action = 'na.fail')
mod2_hypo_ch4 <- glm(ch4_umolL_norm ~ BIX + DOC_mgL + DO + co2_umolL_zt, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')
mod1_hypo_ch4_do <- glm(ch4_umolL_norm ~ DO, 
                     data = hypo_data_2, family = gaussian, na.action = 'na.fail')

mod1_hypo_co2_lag <- glm(co2_umolL ~ co2_umolL_ARLag1, 
                         data = hypo_data_2, family = gaussian, na.action = 'na.fail')
glm_hypo_co2_lag <- dredge(mod1_hypo_co2_lag)

mod1_hypo_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + DOC_mgL + Flow_cms_norm, 
                     data = hypo_data_2, family = gaussian, na.action = 'na.fail')
mod2_hypo_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + DOC_mgL + DO + Flow_cms_norm, 
                     data = hypo_data_2, family = gaussian, na.action = 'na.fail')


round((rsq(mod1_epi_ch4, type = 'sse')), digits = 2)

round((rsq(mod1_epi_co2_lag, type = 'sse')), digits = 2)

round((rsq(mod1_epi_co2, type = 'sse')), digits = 2)
round((rsq(mod2_epi_co2, type = 'sse')), digits = 2)
round((rsq(mod3_epi_co2, type = 'sse')), digits = 2)
round((rsq(mod4_epi_co2, type = 'sse')), digits = 2)

round((rsq(mod1_hypo_ch4_lag, type = 'sse')), digits = 2)

round((rsq(mod1_hypo_ch4, type = 'sse')), digits = 2)
round((rsq(mod2_hypo_ch4, type = 'sse')), digits = 2)
round((rsq(mod1_hypo_ch4_do, type = 'sse')), digits = 2)

round((rsq(mod1_hypo_co2_lag, type = 'sse')), digits = 2)

round((rsq(mod1_hypo_co2, type = 'sse')), digits = 2)
round((rsq(mod2_hypo_co2, type = 'sse')), digits = 2)


pred1_mod1_epi_ch4 <- predict(mod1_epi_ch4,newdata = epi_data_2)

pred1_mod1_epi_co2_lag <- predict(mod1_epi_co2_lag,newdata=epi_data_2)

pred1_mod1_epi_co2 <- predict(mod1_epi_co2,newdata = epi_data_2)
pred1_mod2_epi_co2 <- predict(mod2_epi_co2,newdata = epi_data_2)
pred1_mod3_epi_co2 <- predict(mod3_epi_co2,newdata = epi_data_2)
pred1_mod4_epi_co2 <- predict(mod4_epi_co2,newdata = epi_data_2)

pred1_mod1_hypo_ch4_lag <- predict(mod1_hypo_ch4_lag,newdata = hypo_data_2)

pred1_mod1_hypo_ch4 <- predict(mod1_hypo_ch4,newdata = hypo_data_2)
pred1_mod2_hypo_ch4 <- predict(mod2_hypo_ch4,newdata = hypo_data_2)

pred1_mod1_hypo_co2_lag <- predict(mod1_hypo_co2_lag,newdata = hypo_data_2)

pred1_mod1_hypo_co2 <- predict(mod1_hypo_co2,newdata = hypo_data_2)
pred1_mod2_hypo_co2 <- predict(mod2_hypo_co2,newdata = hypo_data_2)


round(rmse(pred1_mod1_epi_ch4, epi_data_2$ch4_umolL), digits = 1)

round(rmse(pred1_mod1_epi_co2_lag, epi_data_2$co2_umolL_norm),digits = 1)

round(rmse(pred1_mod1_epi_co2, epi_data_2$co2_umolL_norm), digits = 1)
round(rmse(pred1_mod2_epi_co2, epi_data_2$co2_umolL_norm), digits = 1)
round(rmse(pred1_mod3_epi_co2, epi_data_2$co2_umolL_norm), digits = 1)
round(rmse(pred1_mod4_epi_co2, epi_data_2$co2_umolL_norm), digits = 1)

round(rmse(pred1_mod1_hypo_ch4_lag, hypo_data_2$ch4_umolL_norm), digits = 1)

round(rmse(pred1_mod1_hypo_ch4, hypo_data_2$ch4_umolL_norm), digits = 1)
round(rmse(pred1_mod2_hypo_ch4, hypo_data_2$ch4_umolL_norm), digits = 1)

round(rmse(pred1_mod1_hypo_co2_lag, hypo_data_2$co2_umolL), digits = 1)

round(rmse(pred1_mod1_hypo_co2, hypo_data_2$co2_umolL), digits = 1)
round(rmse(pred1_mod2_hypo_co2, hypo_data_2$co2_umolL), digits = 1)


# Then need to 'un-transform' variables
epi_data_uncorr <- epi_data[complete.cases(epi_data),]
epi_data_uncorr <- epi_data_uncorr %>% 
  mutate(co2_umolL_norm = log(co2_umolL)) %>% 
  mutate(Flow_cms_norm = log(Flow_cms)) %>% 
  mutate(co2_umolL_ARLag1_norm = log(co2_umolL_ARLag1)) %>% 
  select(-c(co2_umolL,Flow_cms,co2_umolL_ARLag1))
# CH4 lag term
(select_glm_epi_ch4$ch4_umolL_ARLag1*sd(epi_data_uncorr$ch4_umolL_ARLag1))+mean(epi_data_uncorr$ch4_umolL_ARLag1)
# CO2 Epi
(select_glm_epi_co2$DO*sd(epi_data_uncorr$DO))+mean(epi_data_uncorr$DO)
(select_glm_epi_co2$Flora_ugL*sd(epi_data_uncorr$Flora_ugL))+mean(epi_data_uncorr$Flora_ugL)
(select_glm_epi_co2$Flow_cms_norm*sd(epi_data_uncorr$Flow_cms_norm))+mean(epi_data_uncorr$Flow_cms_norm)
# CH4 Hypo
hypo_data_uncorr <- hypo_data[complete.cases(hypo_data),]
hypo_data_uncorr <- hypo_data_uncorr %>% 
  mutate(Flow_cms_norm = log(Flow_cms)) %>% 
  mutate(Flora_ugL_norm = log(Flora_ugL)) %>% 
  mutate(ch4_umolL_norm = log(ch4_umolL)) %>% 
  mutate(ch4_umolL_ARLag1_norm = log(ch4_umolL_ARLag1)) %>% 
  select(-c(Flow_cms,Flora_ugL,ch4_umolL,ch4_umolL_ARLag1))

(select_glm_hypo_ch4$BIX*sd(hypo_data_uncorr$BIX))+mean(hypo_data_uncorr$BIX)
(select_glm_hypo_ch4$co2_umolL*sd(hypo_data_uncorr$co2_umolL))+mean(hypo_data_uncorr$co2_umolL)
(select_glm_hypo_ch4$DO*sd(hypo_data_uncorr$DO))+mean(hypo_data_uncorr$DO)
(select_glm_hypo_ch4$DOC_mgL*sd(hypo_data_uncorr$DOC_mgL))+mean(hypo_data_uncorr$DOC_mgL)
# CO2 Hypo
(select_glm_hypo_co2$co2_umolL_ARLag1*sd(hypo_data_uncorr$co2_umolL_ARLag1))+mean(hypo_data_uncorr$co2_umolL_ARLag1)
(select_glm_hypo_co2$DO*sd(hypo_data_uncorr$DO))+mean(hypo_data_uncorr$DO)
(select_glm_hypo_co2$DOC_mgL*sd(hypo_data_uncorr$DOC_mgL))+mean(hypo_data_uncorr$DOC_mgL)

############ Meta - NEED TO NORMALIZE GHG DATA
meta_data_2 <- meta_data[complete.cases(meta_data),]
meta_data_2 <- meta_data_2 %>% select(-c(Date,BIX_ARLag1,HIX_ARLag1,HIX))
meta_data_2 <- as.data.frame(meta_data_2)
meta_data_2 <- meta_data_2 %>% 
  mutate(Flow_cms_norm = log(Flow_cms)) %>% 
  select(-c(Flow_cms))
meta_data_2 <- scale(meta_data_2)

chart.Correlation(meta_data_2,histogram=TRUE,method=c("spearman"))
meta_data_2 <- as.data.frame(meta_data_2)

# Meta CH4 = ch4_lag + Flora + DO + DOC + BIX
# Remove: co2_lag, temp, co2, flow
meta_ch4_corr <- meta_data_2 %>% select(ch4_umolL,ch4_umolL_ARLag1,Flora_ugL,DO,DOC_mgL,BIX)
chart.Correlation(meta_ch4_corr,histogram=TRUE,method=c("spearman"))

model_meta_ch4 <- glm(ch4_umolL ~ ch4_umolL_ARLag1 + Flora_ugL + DO + DOC_mgL + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_ch4 <- dredge(model_meta_ch4,rank="AICc")

select_glm_meta_ch4 <- subset(glm_meta_ch4,delta<2)

# Meta CO2 = co2_lag + Flora + DOC + BIX
# Remove: ch4_lag, flow, DO, temp, CH4
meta_co2_corr <- meta_data_2 %>% select(co2_umolL,co2_umolL_ARLag1,Flora_ugL,DOC_mgL,BIX)
chart.Correlation(meta_co2_corr,histogram=TRUE,method=c("spearman"))

model_meta_co2 <- glm(co2_umolL ~ co2_umolL_ARLag1 + Flora_ugL + DOC_mgL + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_co2 <- dredge(model_meta_co2,rank="AICc")

select_glm_meta_co2 <- subset(glm_meta_co2,delta<2)

############################## DOM data
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

# Epi BIX = bix_lag + Flow, DOC, Flora, Temp
# Remove: HIX_lag, CO2_lag, CH4_lag, DO, HIX, CO2, ch4
epi_bix_corr <- epi_data_2 %>% select(BIX,BIX_ARLag1,Flow_cms_norm,DOC_mgL,Flora_ugL,temp)
chart.Correlation(epi_bix_corr,histogram=TRUE,method=c("spearman"))

model_epi_bix <- glm(BIX ~ BIX_ARLag1 + DOC_mgL + Flora_ugL + temp +
                       Flow_cms_norm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_bix <- dredge(model_epi_bix,rank = "AICc")

select_glm_epi_bix <- subset(glm_epi_bix,delta<2)

# Epi HIX = hix_lag + flow + DOC + Flora + CO2
# Remove: BIX_lag, co2_lag, ch4_lag, DO, BIX, temp, ch4
epi_hix_corr <- epi_data_2 %>% select(HIX,HIX_ARLag1,Flow_cms_norm,DOC_mgL,Flora_ugL,co2_umolL_norm)
chart.Correlation(epi_hix_corr,histogram=TRUE,method=c("spearman"))

model_epi_hix <- glm(HIX ~ HIX_ARLag1 + DOC_mgL + Flora_ugL + temp + co2_umolL_norm +
                       Flow_cms_norm, data = epi_data_2, 
                     family = gaussian, na.action = 'na.fail')

glm_epi_hix <- dredge(model_epi_hix,rank = "AICc")

select_glm_epi_hix <- subset(glm_epi_hix,delta<2)

# Meta BIX = bix_lag + temp + HIX + DOC
# Remove: co2_lag, ch4_lag, hix_lag, DO, CO2, Flora, Ch4
meta_bix_corr <- meta_data_2 %>% select(BIX,BIX_ARLag1,Flow_cms_norm,temp,DOC_mgL,HIX)
chart.Correlation(meta_bix_corr,histogram=TRUE,method=c("spearman"))

model_meta_bix <- glm(BIX ~ BIX_ARLag1 + temp + DOC_mgL + HIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_bix <- dredge(model_meta_bix,rank="AICc")

select_glm_meta_bix <- subset(glm_meta_bix,delta<2)

# Meta HIX = hix_lag + BIX + temp + flow + DOC + BIX
# Remove: co2_lag, ch4_lag, bix_lag, DO, CO2, flora, ch4
meta_hix_corr <- meta_data_2 %>% select(HIX,HIX_ARLag1,Flow_cms_norm,temp,DOC_mgL,BIX)
chart.Correlation(meta_hix_corr,histogram=TRUE,method=c("spearman"))

model_meta_hix <- glm(HIX ~ HIX_ARLag1 + temp + DOC_mgL + Flow_cms_norm + BIX, data = meta_data_2, 
                      family = gaussian, na.action = 'na.fail')

glm_meta_hix <- dredge(model_meta_hix,rank="AICc")

select_glm_meta_hix <- subset(glm_meta_hix,delta<2)

# Hypo BIX = bix_lag + DOC + CO2 + ch4 + HIX
# Remove: ch4_lag, co2_lag, hix_lag, temp, flow, do, flora
hypo_bix_corr <- hypo_data_2 %>% select(BIX,BIX_ARLag1,DOC_mgL,co2_umolL,ch4_umolL_norm,HIX)
chart.Correlation(hypo_bix_corr,histogram=TRUE,method=c("spearman"))

model_hypo_bix <- glm(BIX ~ BIX_ARLag1 + HIX + DOC_mgL + ch4_umolL_norm + co2_umolL, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_bix <- dredge(model_hypo_bix,rank="AICc")

select_glm_hypo_bix <- subset(glm_hypo_bix,delta<2)

# Hypo HIX = hix_lag + Flow + DO + DOC + CO2 + BIX
# Remove: ch4_lag, co2_lag, bix_lag, temp, flora, ch4
hypo_hix_corr <- hypo_data_2 %>% select(HIX,HIX_ARLag1,Flow_cms_norm,DO,DOC_mgL,co2_umolL,BIX)
chart.Correlation(hypo_hix_corr,histogram=TRUE,method=c("spearman"))

model_hypo_hix <- glm(HIX ~ HIX_ARLag1 + DOC_mgL + DO + Flow_cms_norm + co2_umolL + BIX, 
                      data = hypo_data_2, family = gaussian, na.action = 'na.fail')

glm_hypo_hix <- dredge(model_hypo_hix,rank="AICc")

select_glm_hypo_hix <- subset(glm_hypo_hix,delta<2)
