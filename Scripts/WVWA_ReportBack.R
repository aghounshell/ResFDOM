### Getting data together for 2020 WVWA Report back
### DOC, Abs, and FL data for DBP formation
### Focus on: long-term DOC data; SUVA and Humic C data for 2019

setwd("C:/Users/ahoun/Desktop/ResFDOM")

# Load in libraries
pacman::p_load(tidyverse,ggplot2)

install.packages('tidyverse', INSTALL_opts = c('--no-lock'))
install.packages('ggplot2', INSTALL_opts = c('--no-lock'))
install.packages('ggpubr', INSTALL_opts = c('--no-lock'))

# Load long-term DOC data
fcrchem <- read.csv("./Data/chem.csv", header=T) %>%
  select(Reservoir:DN_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  filter(Site == 50) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  #filter(DateTime > as.POSIXct("2018-12-31") & DateTime < as.POSIXct("2020-01-01")) %>% 
  filter(Depth_m == 0.1|Depth_m == 5.0|Depth_m == 9.0)

fcrchem <- fcrchem %>% select(-c("Reservoir")) %>% 
  rename(Date = DateTime, Depth = Depth_m, Station = Site)

# Plot long-term DOC
# ggplot(fcrchem )