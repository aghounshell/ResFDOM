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

# Epi data
epi <- data %>% filter(Station == 50 & Depth == 0.1)

# Meta data
meta <- data %>% filter(Station == 50 & Depth == 5)

# Hypo data
hypo <- data %>% filter(Station == 50 & Depth == 9)

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

data_2 <- rbind(epi,meta,hypo) %>% arrange(Date)

############## PCA at station 50 only - All data used for AR Model ######################
data_50 <- data_2 %>% select(Date,Depth,BIX,HIX,ch4_umolL,co2_umolL,DOC_mgL,temp,DO,Flora_ugL,Flow_cms)
data_50 <- data_50[complete.cases(data_50),]

# Check correlations
corr_50 <- data_50 %>% select(-c(Date,Depth))
chart.Correlation(corr_50,histogram=TRUE,method=c("spearman"))

# Correct data for right-skew
skew_50 <- (corr_50)^(1/3)
chart.Correlation(skew_50,histogram=TRUE,method=c("spearman"))

# Normalize data using scale
scale_50 <- scale(skew_50)
chart.Correlation(scale_50,histogram=TRUE,method=c("spearman"))

# Conduct PCA on DOM data - normalized and transformed; no outliers removed
cpca_50 <- rda(scale_50)
summary(cpca_50,axes=0)
plot(cpca_50)
text(cpca_50)
screeplot(cpca_50, bstick = TRUE)

# Extract species scores for scaling 2
c50_spe_sc2 <- scores(cpca_50, choices=1:3, display="sp", scaling=2)

# Plot results on two axes
par(mar=c(5.1,5.1,4.1,2.1))

data_50$Depth <- factor(data_50$Depth, levels=c("0.1", "5", "9"))
with(data_50,levels(Depth))
colvec<-c("#7FC6A4","#7EBDC2","#393E41")
with(data_50,levels(Depth))
sq<-c(21,22,23)

plot(cpca_50,type="n",scaling=2,xlab="PC1 (37% var. explained)",ylab="PC2 (29% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(data_50,points(cpca_50,display="sites",col=c("black","black","black"),scaling=2,pch=sq[Depth],
                      bg=colvec[Depth],cex=1.5))
with(data_50,legend("topright",legend=levels(Depth),bty="n",col=c("black","black","black"),
                   pch=c(21,22,23),pt.bg=colvec,cex=1.3))
arrows(0, 0, c50_spe_sc2[,1], c50_spe_sc2[,2], angle=20, col="black")
text(cpca_50, display = "species", scaling=2, cex = 1, col = "black")

# Separate DOM data
dom <- data %>% select(Date:'M/C') %>% filter(Station == 50)
dom <- dom[complete.cases(dom),]

# Double check correlations
dom_corr <- dom %>% select(-c(Date,Station,Depth))
chart.Correlation(dom_corr,histogram=TRUE,method=c("spearman"))

# Remove Fmax4, Fmax2, M/C
dom_corr <- dom %>% select(-c(Date,Station,Depth,Fmax4,Fmax2,'M/C'))
chart.Correlation(dom_corr,histogram=TRUE,method=c("spearman"))

# Normalize data using scale
dom_scale <- scale(dom_corr)

# Conduct PCA on DOM data - normalized but not transformed; no outliers removed
dom_pca <- rda(dom_scale)
summary(dom_pca,axes=0)
plot(dom_pca)
text(dom_pca)
screeplot(dom_pca, bstick = TRUE)

# Extract species scores for scaling 2
domspe_sc2 <- scores(dom_pca, choices=1:3, display="sp", scaling=2)

# Plot results on two axes
dom$Depth <- factor(dom$Depth, levels=c("0.1", "5", "9"))
with(dom,levels(Depth))
colvec<-c("#7FC6A4","#7EBDC2","#393E41")
with(dom,levels(Depth))
sq<-c(21,22,23)

plot(dom_pca,type="n",scaling=2,xlab="PC1 (41% var. explained)",ylab="PC2 (26% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(dom,points(dom_pca,display="sites",col=c("black","black","black"),scaling=2,pch=sq[Depth],
                bg=colvec[Depth],cex=1.5))
with(dom,legend("topright",legend=levels(Depth),bty="n",col=c("black","black","black"),
                  pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, display = "species", scaling=2, cex = 0.8, col = "black")

############################### All C data ###############################################
c_data <- data %>% select(Date:co2_umolL,DOC_mgL)
c_data <- c_data[complete.cases(c_data),]

c_data$location <- paste(c_data$Station,c_data$Depth)

# Check correlations
c_corr <- c_data %>% select(-c(Date,Station,Depth,location))
chart.Correlation(c_corr,histogram=TRUE,method=c("spearman"))

c_corr <- c_data %>% select(-c(Date,Station,Depth,location,Fmax4,'M/C'))
chart.Correlation(c_corr,histogram=TRUE,method=c("spearman"))

# Correct data for right-skew
c_skew <- (c_corr)^(1/3)
chart.Correlation(c_skew,histogram=TRUE,method=c("spearman"))

# Normalize data using scale
c_scale <- scale(c_skew)
chart.Correlation(c_scale,histogram=TRUE,method=c("spearman"))

# Conduct PCA on DOM data - normalized and transformed; no outliers removed
c_pca <- rda(c_scale)
summary(c_pca,axes=0)
plot(c_pca)
text(c_pca)
screeplot(c_pca, bstick = TRUE)
# NOTE! First three axes are important!!!

# Extract species scores for scaling 2
cspe_sc2 <- scores(c_pca, choices=1:3, display="sp", scaling=2)

# Plot results
c_data$location<-factor(c_data$location, levels=c("50 0.1", "50 5", "50 9", "100 0.1", "200 0.1"))
with(c_data,levels(location))
colvec<-c("#7FC6A4","#7EBDC2","#393E41","#F0B670","#FE5F55")
with(c_data,levels(location))
sq<-c(21,22,23,24,25)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(c_pca,type="n",scaling=2,xlab="PC1 (33% var. explained)",ylab="PC2 (23% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(c_data,points(dom_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[location],
                  bg=colvec[location],cex=1.5))
with(c_data,legend("topright",legend=levels(location),bty="n",col=c("black","black","black","black","black"),
                  pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, cspe_sc2[,1], cspe_sc2[,2], angle=20, col="black")
text(c_pca, display = "species", scaling=2, cex = 1, col = "black")

plot(c_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (33% var. explained)",ylab="PC3 (18% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(c_data,points(c_pca,choices=c(1,3),display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[location],
                     bg=colvec[location],cex=1.5))
with(c_data,legend("bottomleft",legend=levels(location),bty="n",col=c("black","black","black","black","black"),
                     pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, cspe_sc2[,1], cspe_sc2[,3], angle=20, col="black")
text(c_pca,choices=c(1,3), display = "species", scaling=2, cex = 1, col = "black")

plot(c_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (23% var. explained)",ylab="PC3 (18% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(c_data,points(c_pca,choices=c(2,3),display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[location],
                   bg=colvec[location],cex=1.5))
with(c_data,legend("bottomleft",legend=levels(location),bty="n",col=c("black","black","black","black","black"),
                   pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, cspe_sc2[,2], cspe_sc2[,3], angle=20, col="black")
text(c_pca,choices=c(2,3), display = "species", scaling=2, cex = 1, col = "black")

############## PCA at station 50 only? All C data ######################
c_data_50 <- data %>% select(Date:co2_umolL) %>% filter(Station == "50")
c_data_50 <- c_data_50[complete.cases(c_data_50),]

# Check correlations
c_corr_50 <- c_data_50 %>% select(-c(Date,Station,Depth))
chart.Correlation(c_corr_50,histogram=TRUE,method=c("spearman"))

c_corr_50 <- c_data_50 %>% select(-c(Date,Station,Depth,Fmax2,Fmax4,'M/C'))
chart.Correlation(c_corr_50,histogram=TRUE,method=c("spearman"))

# Correct data for right-skew
c_skew_50 <- (c_corr_50)^(1/3)
chart.Correlation(c_skew_50,histogram=TRUE,method=c("spearman"))

# Normalize data using scale
c_scale_50 <- scale(c_skew_50)
chart.Correlation(c_scale_50,histogram=TRUE,method=c("spearman"))

# Conduct PCA on DOM data - normalized and transformed; no outliers removed
c_pca_50 <- rda(c_scale_50)
summary(c_pca_50,axes=0)
plot(c_pca_50)
text(c_pca_50)
screeplot(c_pca_50, bstick = TRUE)

# Extract species scores for scaling 2
c50_spe_sc2 <- scores(c_pca_50, choices=1:3, display="sp", scaling=2)

# Plot results on two axes
c_data_50$Depth <- factor(c_data_50$Depth, levels=c("0.1", "5", "9"))
with(c_data_50,levels(Depth))
colvec<-c("#7FC6A4","#7EBDC2","#393E41")
with(c_data_50,levels(Depth))
sq<-c(21,22,23)

plot(c_pca_50,type="n",scaling=2,xlab="PC1 (33% var. explained)",ylab="PC2 (30% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(c_data_50,points(c_pca_50,display="sites",col=c("black","black","black"),scaling=2,pch=sq[Depth],
                bg=colvec[Depth],cex=1.5))
with(c_data_50,legend("topright",legend=levels(Depth),bty="n",col=c("black","black","black"),
                pt.bg=colvec,cex=1.3))
arrows(0, 0, c50_spe_sc2[,1], c50_spe_sc2[,2], angle=20, col="black")
text(c_pca_50, display = "species", scaling=2, cex = 1, col = "black")

###########################################################################
# Now conduct RDA using C data and Env Data (DO, Temp, Flora) for station 50
rda_data <- data %>% 
  select(Date:Fmax1,Fmax3,BIX:'A/C',ch4_umolL,co2_umolL,DOC_mgL,temp,DO,Flora_ugL) %>% 
  filter(Station == "50")
rda_data <- rda_data[complete.cases(rda_data),]

rda_c_data <- rda_data %>% select(Fmax1,Fmax3,BIX:'A/C',ch4_umolL,co2_umolL)
rda_env_data <- rda_data %>% select(temp,DO,Flora_ugL,DOC_mgL)

# Double check correlations
chart.Correlation(rda_c_data,histogram=TRUE,method=c("spearman"))
chart.Correlation(rda_env_data,histogram=TRUE,method=c("spearman"))

# Correct data for right-skew
rda_c_skew <- (rda_c_data)^(1/3)

# Normalize data using scale
rda_c_scale <- as.data.frame(scale(rda_c_skew))
chart.Correlation(rda_c_scale,histogram=TRUE,method=c("spearman"))
rda_env_scale <- as.data.frame(scale(rda_env_data))

# Conduct RDA
c_rda <- rda(rda_c_scale~.,rda_env_scale,scale=FALSE)

# Check global model (R2 = 0.4)
(R2a_all <- RsquareAdj(c_rda)$adj.r.squared)

# Test of all canoical axes from full rda
anova(c_rda,by="axis",permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
c_rda_forsel <- forward.sel(rda_c_scale,rda_env_scale,Xscale=FALSE,Yscale=FALSE,
                            Ycenter=FALSE,adjR2thresh=R2a_all)

# All variables are important (DO, DOC, Flora)
# C RDA model: DO, DOC, Flora
c_rda_final <- rda(rda_c_scale ~ rda_env_scale$DO + rda_env_scale$DOC_mgL + rda_env_scale$Flora_ugL,scale=FALSE)

# Extract scores
c_rspe_sc2 <- scores(c_rda_final,display="sp",choices=c(1,2),scaling=2)
c_rbp.sc2 <- scores(c_rda_final,display="bp",choices=c(1,2),scaling=2)

# Plot
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,1))

rda_data$Depth <- factor(rda_data$Depth, levels=c("0.1", "5", "9"))
with(rda_data,levels(Depth))
colvec<-c("#7FC6A4","#7EBDC2","#393E41")
with(rda_data,levels(Depth))
sq<-c(21,22,23)

plot(c_rda_final,scaling=2,display="sites",xlab="RDA1 (56% fitted,23% total var.)",
     ylab="RDA2 (40% fitted, 16% total var.)",cex.axis=1.5,cex.lab=1.5)
with(rda_data,points(c_rda_final,display="sites",col=c("black","black","black"),scaling=2,pch=sq[Depth],
                  bg=colvec[Depth],cex=1.5))
with(rda_data,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black"),
                  pch=c(21,22,23),pt.bg=colvec,cex=1.5))
arrows(0,0,c_rspe_sc2[,1], c_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(c_rda_final,display = "species", scaling=2, cex = 1.5, 
     col = "black")
text(c_rda_final,display="bp",scaling=2,cex=1.5,col="black")
