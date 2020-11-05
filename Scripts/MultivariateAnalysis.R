### Script to conduct PCA/RDA on 2019 data
### 1. Determine a 'universal' metric for DOM
### 2. Initially assess relationships between parameters
### A Hounshell, 05 Nov 2020

setwd("C:/Users/ahoun/OneDrive/Desktop/ResFDOM")

# Load packages
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc,ggplot2,tidyverse,vegan3d,
               scatterplot3d,rgl,car)

# Load data
data <- read_csv("./Data/20201105_All_Data.csv")
data$Date <- as.POSIXct(strptime(data$Date, "%Y-%m-%d", tz = "EST"))

# Separate DOM data
dom <- data %>% select(Date:'M/C')
dom <- dom[complete.cases(dom),]

# Double check correlations
dom_corr <- dom %>% select(-c(Date,Station,Depth))
chart.Correlation(dom_corr,histogram=TRUE,method=c("spearman"))

