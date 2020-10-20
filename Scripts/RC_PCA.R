### Script to conduct super quick and dirty PCA with RC day nutrient+Chla data
### NOTE: NO OM data for this analysis!!!
### 20 Oct 2020, A Hounshell

# Load libraries
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc,ggplot2,tidyverse,vegan3d,
               scatterplot3d,rgl,ggpubr)

# Load data (compiled by WMW - nutrient, Chla data only)
data <- read_csv('C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/continuum_ww.csv')
data$DateTime <- as.POSIXct(strptime(data$DateTime, "%Y-%m-%d", tz = "EST"))

data <- data %>% select(Reservoir,Site,DateTime,TN_ugL,TP_ugL,NH4_ugL,NO3NO2_ugL,SRP_ugL,DOC_mgL,
                        Chla_ugL,connectivity,overflow)

data <- data[complete.cases(data),]

# Select data for PCA
data_all <- data %>% select(TN_ugL:Chla_ugL)

# Scale data (z scores)
data_all_scale <- scale(data_all)

# Check for linear correlations
chart.Correlation(data_all, histogram=TRUE, method=c("pearson"))
# Remove TN (correlated w/ Chla); Data is right skewed - try square root?

pca_all <- data %>% select(TP_ugL:Chla_ugL)
pca_all <- (pca_all)^(1/2)
pca_all <- scale(pca_all)

chart.Correlation(pca_all,histogram=TRUE,method=c("pearson"))

# Conduct pca: removed TN; transformed with square root method; z-score scaled
data_pca <- rda(pca_all)
summary(data_pca,axes=0)
plot(data_pca)
text(data_pca)
screeplot(data_pca, bstick = TRUE)

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
spe_sc2 <- scores(data_pca, choices=1:3, display="sp", scaling=2)

# PCA results plotted in 3D by Reservoir + site
data$location <- paste(data$Reservoir,data$Site)
data$location <- factor(data$location, levels=c("BVR 1","BVR 20","BVR 30","BVR 45","BVR 50", "BVR 100", 
                                                "BVR 200", "FCR 1", "FCR 20", "FCR 30", "FCR 45", 
                                                "FCR 50", "FCR 99", "FCR 100", "FCR 101", "FCR 102", "FCR 200"))


with(data,levels(location))
colvec<-c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#EEAA55","#7DAF4B","#BFACC8",
          "#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#E7804B","#DA2C38","#7DAF4B")
sq<-c(22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21)

graphs <- ordiplot3d(data_pca,display="sites",choices=1:3,scaling=2,xlab="PC1 (36%)",
                  zlab="PC3 (15%)",ylab="PC2 (28%)")
points(graphs,"points",pch=sq[data$location],col=c("black","black","black","black"),bg=colvec[data$location],
       cex=0.7)
text(graphs$xyz.convert(spe_sc2),rownames(spe_sc2),cex=0.8,xpd=TRUE)

# 2D plot
data$location <- factor(data$location, levels=c("BVR 1","BVR 20","BVR 30","BVR 45","BVR 50", "BVR 100", 
                                                "BVR 200", "FCR 1", "FCR 20", "FCR 30", "FCR 45", 
                                                "FCR 50", "FCR 99", "FCR 100", "FCR 101", "FCR 102", "FCR 200"))
with(data,levels(location))
colvec<-c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#EEAA55","#7DAF4B","#BFACC8",
          "#393E41","#5C7E82","#7EBDC2","#7FC6A4","#F4D35E","#EEAA55","#E7804B","#DA2C38","#7DAF4B")
sq<-c(22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21)

plot(data_pca,type="n",scaling=2,xlab="PC1 (36% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5,ylim=c(-1.5,1.5))
with(data,points(data_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[location],
                     bg=colvec[location],cex=1.5))
with(data,legend("topleft",legend=levels(location),bty="n",col=c("black","black","black","black"),
                     pch=c(22,22,22,22,22,22,22,21,21,21,21,21,21,21,21,21,21),pt.bg=colvec,cex=1.3))
arrows(0, 0, spe_sc2[,1], spe_sc2[,2], angle=20, col="black")
text(data_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")


# 2D plot by date
data$DateTime <- factor(data$DateTime, levels = c("2019-04-29", "2019-05-30", "2019-06-27", "2019-07-18", 
                                                  "2019-08-22", "2019-09-20", "2019-10-04"))
with(data,levels(DateTime))
colvec <- c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#EEAA55","#7DAF4B")
sq <- c(22,22,22,22,22,22,22)

plot(data_pca,type="n",scaling=2,xlab="PC1 (36% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5,ylim=c(-1.5,1.5))
with(data,points(data_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[DateTime],
                 bg=colvec[DateTime],cex=1.5))
with(data,legend("topleft",legend=levels(DateTime),bty="n",col=c("black","black","black","black"),
                 pch=c(22,22,22,22,22,22,22),pt.bg=colvec,cex=1.3))
arrows(0, 0, spe_sc2[,1], spe_sc2[,2], angle=20, col="black")
text(data_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")


########################## Conduct PCA on JUST reservoir samples ################################
res_data <- read_csv('C:/Users/ahoun/OneDrive/Desktop/ResFDOM/Data/continuum_ww_res.csv')
res_data$DateTime <- as.POSIXct(strptime(res_data$DateTime, "%Y-%m-%d", tz = "EST"))

res_data <- res_data %>% select(Reservoir,Site,DateTime,TN_ugL,TP_ugL,NH4_ugL,NO3NO2_ugL,SRP_ugL,DOC_mgL,
                        Chla_ugL,connectivity,overflow)

res_data <- res_data[complete.cases(res_data),]

# Same data procedures as above
res_all <- res_data %>% select(TN_ugL:Chla_ugL)
res_all <- (res_all)^(1/2)
res_all <- scale(res_all)

chart.Correlation(res_all,histogram=TRUE,method=c("pearson"))

# Remove TN,TP
res_all <- res_data %>% select(NH4_ugL:Chla_ugL)
res_all <- (res_all)^(1/2)
res_all <- scale(res_all)

chart.Correlation(res_all,histogram=TRUE,method=c("pearson"))

# Conduct pca: removed TN & TP; transformed with square root method; z-score scaled
res_pca <- rda(res_all)
summary(res_pca,axes=0)
plot(res_pca)
text(res_pca)
screeplot(res_pca, bstick = TRUE)

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
res_spe_sc2 <- scores(res_pca, choices=1:3, display="sp", scaling=2)

# PCA results plotted in 3D by Reservoir + site
res_data$location <- paste(res_data$Reservoir,res_data$Site)
res_data$location <- factor(res_data$location, levels=c("BVR 20","BVR 30","BVR 45","BVR 50", "FCR 20", 
                                                        "FCR 30", "FCR 45","FCR 50"))


with(res_data,levels(location))
colvec<-c("#393E41","#5C7E82","#7EBDC2","#7FC6A4",
          "#393E41","#5C7E82","#7EBDC2","#7FC6A4")
sq<-c(22,22,22,22,21,21,21,21,21)

plot(res_pca,type="n",scaling=2,xlab="PC1 (50% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(res_data,points(res_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[location],
                 bg=colvec[location],cex=1.5))
with(res_data,legend("topleft",legend=levels(location),bty="n",col=c("black","black","black","black"),
                 pch=c(22,22,22,22,21,21,21,21,21),pt.bg=colvec,cex=1.3))
arrows(0, 0, res_spe_sc2[,1], res_spe_sc2[,2], angle=20, col="black")
text(res_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")

# 2D plot by date
res_data$DateTime <- factor(res_data$DateTime, levels = c("2019-04-29", "2019-05-30", "2019-06-27", "2019-07-18", 
                                                  "2019-08-22", "2019-09-20", "2019-10-04"))
with(res_data,levels(DateTime))
colvec <- c("#BFACC8","#393E41","#5C7E82","#7EBDC2","#7FC6A4","#EEAA55","#7DAF4B")
sq <- c(22,22,22,22,22,22,22)

plot(res_pca,type="n",scaling=2,xlab="PC1 (50% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(res_data,points(res_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[DateTime],
                 bg=colvec[DateTime],cex=1.5))
with(res_data,legend("topleft",legend=levels(DateTime),bty="n",col=c("black","black","black","black"),
                 pch=c(22,22,22,22,22,22,22),pt.bg=colvec,cex=1.3))
arrows(0, 0, res_spe_sc2[,1], res_spe_sc2[,2], angle=20, col="black")
text(res_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")
