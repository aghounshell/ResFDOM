### Script to conduct data analyses on all compiled data (see All_Data.R)
### A Hounshell, 03 Nov 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,PerformanceAnalytics)

# Load compiled data
data <- read_csv("./Data/20201103_All_Data.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%Y-%m-%d", tz = "EST"))

# Data for station 50 ONLY
data_50 <- data %>% 
  filter(Station == 50) %>% 
  select(Depth,Fmax1:DO,Flora_ugL)
data_50 <- data_50[complete.cases(data_50),]

# Check autocorrelation
chart.Correlation(data_50,histogram=TRUE,method=c("pearson"))

# Separate by depth (for now)
data_epi <- data %>% 
  filter(Station == 50 & Depth == 0.1) %>% 
  select(Date,Fmax1:DO,Flora_ugL,Rain_Total_mm,ShortwaveRadiationUp_Average_W_m2)
data_epi <- data_epi[complete.cases(data_epi),]

data_meta <- data %>% 
  filter(Station == 50 & Depth == 5) %>% 
  select(Date,Fmax1:DO,Flora_ugL)
data_meta <- data_meta[complete.cases(data_meta),]

data_hypo <- data %>% 
  filter(Station == 50 & Depth == 9) %>% 
  select(Date,Fmax1:DO,Flora_ugL)
data_hypo <- data_hypo[complete.cases(data_hypo),]

data_inflow <- data %>% 
  filter(Station == 100 | Station == 200) %>% 
  select(Date,Fmax1:DN_mgL,Flow_cms)
data_inflow <- data_inflow[complete.cases(data_inflow),]
