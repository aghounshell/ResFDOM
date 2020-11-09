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
dom <- data %>% select(Date:'M/C') %>% filter(Station == 50)
dom <- dom[complete.cases(dom),]

# Double check correlations
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

c_corr <- c_data %>% select(-c(Date,Station,Depth,location,Fmax4,'M/C'))
chart.Correlation(c_corr,histogram=TRUE,method=c("spearman"))

# Normalize data using scale
c_scale <- scale(c_corr)

# Conduct PCA on DOM data - normalized but not transformed; no outliers removed
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
text(c_pca, display = "species", scaling=2, cex = 0.8, col = "black")

plot(c_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (33% var. explained)",ylab="PC3 (18% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(c_data,points(c_pca,choices=c(1,3),display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[location],
                     bg=colvec[location],cex=1.5))
with(c_data,legend("bottomleft",legend=levels(location),bty="n",col=c("black","black","black","black","black"),
                     pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, cspe_sc2[,1], cspe_sc2[,3], angle=20, col="black")
text(c_pca,choices=c(1,3), display = "species", scaling=2, cex = 0.8, col = "black")

plot(c_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (23% var. explained)",ylab="PC3 (18% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(c_data,points(c_pca,choices=c(2,3),display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[location],
                   bg=colvec[location],cex=1.5))
with(c_data,legend("bottomleft",legend=levels(location),bty="n",col=c("black","black","black","black","black"),
                   pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, cspe_sc2[,2], cspe_sc2[,3], angle=20, col="black")
text(c_pca,choices=c(2,3), display = "species", scaling=2, cex = 0.8, col = "black")