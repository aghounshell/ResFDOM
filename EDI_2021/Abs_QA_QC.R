##Script to QA/QC absorbance data
##R Thai, 11 May 2021

# Load abs data
abs <- read.csv("./Data/20210510_ResultsFiles_Abs2019.csv")

#Select data for a254 and a350
abs$Date <- as.POSIXct(strptime(abs$Date, "%m/%d/%Y", tz = "EST"))
a254 <- abs %>% select(Reservoir,Date,Station,Depth,a254)
a350 <- abs %>% select(Reservoir,Date,Station,Depth,a350)

#Filter for a254 and a350 in FCR and BVR for all samples
a254_FCR <- a254 %>% filter(Reservoir=="FCR")
a254_BVR <- a254 %>% filter(Reservoir=="BVR")
a350_FCR <- a350 %>% filter(Reservoir=="FCR")
a350_BVR <- a350 %>% filter(Reservoir=="BVR")


#Calculate mean and sd for a254 and a350 at FCR and BVR
m_a254_FCR <- a254_FCR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_FCR <- a254_FCR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a254_BVR <- a254_BVR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_BVR <- a254_BVR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_FCR <- a350_FCR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_FCR <- a350_FCR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_BVR <- a350_BVR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_BVR <- a350_BVR %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))

#Combine station and depth columns for plotting
m_a254_FCR$Location <- paste(m_a254_FCR$Station,m_a254_FCR$Depth)
sd_a254_FCR$Location <- paste(sd_a254_FCR$Station,sd_a254_FCR$Depth)
m_a254_BVR$Location <- paste(m_a254_BVR$Station,m_a254_BVR$Depth)
sd_a254_BVR$Location <- paste(sd_a254_BVR$Station,sd_a254_BVR$Depth)
m_a350_FCR$Location <- paste(m_a350_FCR$Station,m_a350_FCR$Depth)
sd_a350_FCR$Location <- paste(sd_a350_FCR$Station,sd_a350_FCR$Depth)
m_a350_BVR$Location <- paste(m_a350_BVR$Station,m_a350_BVR$Depth)
sd_a350_BVR$Location <- paste(sd_a350_BVR$Station,sd_a350_BVR$Depth)

#Plot by location for a254 and a350 at FCR and BVR
ggplot(m_a254_FCR,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_FCR,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_FCR$a254,ymax=a254+sd_a254_FCR$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a254_BVR,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_BVR,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_BVR$a254,ymax=a254+sd_a254_BVR$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_FCR,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_FCR,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_FCR$a350,ymax=a350+sd_a350_FCR$a350,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_BVR,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_BVR,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_BVR$a350,ymax=a350+sd_a350_BVR$a350,color=Location))+theme_classic(base_size=15)

#Filter for a254 and a350 in FCR and BVR for samples with a dilution of 1,2
a254_FCR_D1.2 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution<=2) %>% select(Reservoir,Date,Station,Depth,a254)
a254_BVR_D1.2 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution<=2) %>% select(Reservoir,Date,Station,Depth,a254)
a350_FCR_D1.2 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution<=2) %>% select(Reservoir,Date,Station,Depth,a350)
a350_BVR_D1.2 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution<=2) %>% select(Reservoir,Date,Station,Depth,a350)

#Calculate mean and sd for a254 and a350 at FCR and BVR with a dilution of 1,2
m_a254_FCR_D1.2 <- a254_FCR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_FCR_D1.2 <- a254_FCR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a254_BVR_D1.2 <- a254_BVR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_BVR_D1.2 <- a254_BVR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_FCR_D1.2 <- a350_FCR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_FCR_D1.2 <- a350_FCR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_BVR_D1.2 <- a350_BVR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_BVR_D1.2 <- a350_BVR_D1.2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))

#Combine station and depth columns for plotting
m_a254_FCR_D1.2$Location <- paste(m_a254_FCR_D1.2$Station,m_a254_FCR_D1.2$Depth)
sd_a254_FCR_D1.2$Location <- paste(sd_a254_FCR_D1.2$Station,sd_a254_FCR_D1.2$Depth)
m_a254_BVR_D1.2$Location <- paste(m_a254_BVR_D1.2$Station,m_a254_BVR_D1.2$Depth)
sd_a254_BVR_D1.2$Location <- paste(sd_a254_BVR_D1.2$Station,sd_a254_BVR_D1.2$Depth)
m_a350_FCR_D1.2$Location <- paste(m_a350_FCR_D1.2$Station,m_a350_FCR_D1.2$Depth)
sd_a350_FCR_D1.2$Location <- paste(sd_a350_FCR_D1.2$Station,sd_a350_FCR_D1.2$Depth)
m_a350_BVR_D1.2$Location <- paste(m_a350_BVR_D1.2$Station,m_a350_BVR_D1.2$Depth)
sd_a350_BVR_D1.2$Location <- paste(sd_a350_BVR_D1.2$Station,sd_a350_BVR_D1.2$Depth)

#Plot by location for a254 and a350 at FCR and BVR with a dilution of 1,2
ggplot(m_a254_FCR_D1.2,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_FCR_D1.2,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_FCR_D1.2$a254,ymax=a254+sd_a254_FCR_D1.2$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a254_BVR_D1.2,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_BVR_D1.2,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_BVR_D1.2$a254,ymax=a254+sd_a254_BVR_D1.2$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_FCR_D1.2,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_FCR_D1.2,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_FCR_D1.2$a350,ymax=a350+sd_a350_FCR_D1.2$a350,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_BVR_D1.2,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_BVR_D1.2,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_BVR_D1.2$a350,ymax=a350+sd_a350_BVR_D1.2$a350,color=Location))+theme_classic(base_size=15)

#Filter for a254 and a350 in FCR and BVR for samples wih a dilution of 1
a254_FCR_D1 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution==1) %>% select(Reservoir,Date,Station,Depth,a254)
a254_BVR_D1 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution==1) %>% select(Reservoir,Date,Station,Depth,a254)
a350_FCR_D1 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution==1) %>% select(Reservoir,Date,Station,Depth,a350)
a350_BVR_D1 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution==1) %>% select(Reservoir,Date,Station,Depth,a350)

#Calculate mean and sd for a254 and a350 at FCR and BVR with a dilution of 1
m_a254_FCR_D1 <- a254_FCR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_FCR_D1 <- a254_FCR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a254_BVR_D1 <- a254_BVR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_BVR_D1 <- a254_BVR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_FCR_D1 <- a350_FCR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_FCR_D1 <- a350_FCR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_BVR_D1 <- a350_BVR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_BVR_D1 <- a350_BVR_D1 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))

#Combine station and depth columns for plotting
m_a254_FCR_D1$Location <- paste(m_a254_FCR_D1$Station,m_a254_FCR_D1$Depth)
sd_a254_FCR_D1$Location <- paste(sd_a254_FCR_D1$Station,sd_a254_FCR_D1$Depth)
m_a254_BVR_D1$Location <- paste(m_a254_BVR_D1$Station,m_a254_BVR_D1$Depth)
sd_a254_BVR_D1$Location <- paste(sd_a254_BVR_D1$Station,sd_a254_BVR_D1$Depth)
m_a350_FCR_D1$Location <- paste(m_a350_FCR_D1$Station,m_a350_FCR_D1Depth)
sd_a350_FCR_D1$Location <- paste(sd_a350_FCR_D1$Station,sd_a350_FCR_D1$Depth)
m_a350_BVR_D1$Location <- paste(m_a350_BVR_D1$Station,m_a350_BVR_D1$Depth)
sd_a350_BVR_D1$Location <- paste(sd_a350_BVR_D1$Station,sd_a350_BVR_D1$Depth)

#Plot by location for a254 and a350 at FCR and BVR with a dilution of 1
ggplot(m_a254_FCR_D1,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_FCR_D1,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_FCR_D1$a254,ymax=a254+sd_a254_FCR_D1$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a254_BVR_D1,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_BVR_D1,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_BVR_D1$a254,ymax=a254+sd_a254_BVR_D1$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_FCR_D1,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_FCR_D1,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_FCR_D1$a350,ymax=a350+sd_a350_FCR_D1$a350,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_BVR_D1,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_BVR_D1,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_BVR_D1$a350,ymax=a350+sd_a350_BVR_D1$a350,color=Location))+theme_classic(base_size=15)

#Filter for a254 and a350 in FCR and BVR for samples wih a dilution of 2
a254_FCR_D2 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution==2) %>% select(Reservoir,Date,Station,Depth,a254)
a254_BVR_D2 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution==2) %>% select(Reservoir,Date,Station,Depth,a254)
a350_FCR_D2 <- abs %>% filter(Reservoir=="FCR") %>% filter(Dilution==2) %>% select(Reservoir,Date,Station,Depth,a350)
a350_BVR_D2 <- abs %>% filter(Reservoir=="BVR") %>% filter(Dilution==2) %>% select(Reservoir,Date,Station,Depth,a350)

#Calculate mean and sd for a254 and a350 at FCR and BVR with a dilution of 2
m_a254_FCR_D2 <- a254_FCR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_FCR_D2 <- a254_FCR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a254_BVR_D2 <- a254_BVR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a254_BVR_D2 <- a254_BVR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_FCR_D2 <- a350_FCR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_FCR_D2 <- a350_FCR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))
m_a350_BVR_D2 <- a350_BVR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(mean))
sd_a350_BVR_D2 <- a350_BVR_D2 %>% group_by(Date,Station,Depth,Reservoir) %>% summarize_all(funs(sd))

#Combine station and depth columns for plotting
m_a254_FCR_D2$Location <- paste(m_a254_FCR_D2$Station,m_a254_FCR_D2$Depth)
sd_a254_FCR_D2$Location <- paste(sd_a254_FCR_D2$Station,sd_a254_FCR_D2$Depth)
m_a254_BVR_D2$Location <- paste(m_a254_BVR_D2$Station,m_a254_BVR_D2$Depth)
sd_a254_BVR_D2$Location <- paste(sd_a254_BVR_D2$Station,sd_a254_BVR_D2$Depth)
m_a350_FCR_D2$Location <- paste(m_a350_FCR_D2$Station,m_a350_FCR_D2$Depth)
sd_a350_FCR_D2$Location <- paste(sd_a350_FCR_D2$Station,sd_a350_FCR_D2$Depth)
m_a350_BVR_D2$Location <- paste(m_a350_BVR_D2$Station,m_a350_BVR_D2$Depth)
sd_a350_BVR_D2$Location <- paste(sd_a350_BVR_D2$Station,sd_a350_BVR_D2$Depth)

#Plot by location for a254 and a350 at FCR and BVR with a dilution of 2
ggplot(m_a254_FCR_D2,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_FCR_D2,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_FCR_D2$a254,ymax=a254+sd_a254_FCR_D2$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a254_BVR_D2,mapping=aes(x=Date,y=a254,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a254_BVR_D2,mapping=aes(x=Date,y=a254,ymin=a254-sd_a254_BVR_D2$a254,ymax=a254+sd_a254_BVR_D2$a254,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_FCR_D2,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_FCR_D2,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_FCR_D2$a350,ymax=a350+sd_a350_FCR_D2$a350,color=Location))+theme_classic(base_size=15)
ggplot(m_a350_BVR_D2,mapping=aes(x=Date,y=a350,color=Location))+geom_line()+geom_point()+geom_errorbar(m_a350_BVR_D2,mapping=aes(x=Date,y=a350,ymin=a350-sd_a350_BVR_D2$a350,ymax=a350+sd_a350_BVR_D2$a350,color=Location))+theme_classic(base_size=15)

## Combine all abs,EEMs,and PARAFAC data by Reservoir, Date, Station, Depth, Rep, and Dilution ##

# Load abs data
abs <- read.csv("./Data/20210510_ResultsFiles_Abs2019.csv")
abs$Date <- as.POSIXct(strptime(abs$Date, "%m/%d/%Y", tz = "EST"))

#Load eems data
eems <- read.csv("./Data/20210210_ResultsFiles_ResEEMs2019_RAW.csv")
eems$Date <- as.POSIXct(strptime(eems$Date, "%m/%d/%Y", tz = "EST"))

#Load parafac data
parafac <- read.csv("./Data/20210511_ResultsFiles_PARAFAC.csv")
parafac$Date <- as.POSIXct(strptime(parafac$Date, "%m/%d/%Y", tz = "EST"))

#Filter abs data for dilution of 1,2 in FCR,BVR
abs_d <- abs %>% filter(Dilution<=2)
abs_r <- abs_d %>% filter(Reservoir%in%c("FCR","BVR"))
abs_join <- abs_r %>% select(-Sample.Name,-Date.analyzed)

#Filter eems data for dilution of 1,2 in FCR,BVR
eems_d <- eems %>% filter(Dilution<=2)
eems_r <- eems_d %>% filter(Reservoir%in%c("FCR","BVR"))
eems_join <- eems_r %>% select(-Sample.Name,-Date.analyzed,-a254,-a350,-Sr)

#Filter parafac data for dilution of 1,2 in FCR,BVR
parafac_d <- parafac %>% filter(Dilution<=2)
parafac_r <- parafac_d %>% filter(Reservoir%in%c("FCR","BVR"))
parafac_join <- parafac_r%>% select(Reservoir,Date,Station,Depth,Rep,Dilution,Fmax1,Fmax2,Fmax3,Fmax4)

#Join abs and eems data
abs_eems <- full_join(abs_join,eems_join,by=c("Reservoir","Date","Station","Depth","Rep","Dilution"))

#Join the previous dataset (abs_eems) with the parafac data
abs_eems_parafac <- full_join(abs_eems,parafac_join,by=c("Reservoir","Date","Station","Depth","Rep","Dilution"))

#Adjust full dataset column order and names
abs.eems.parafac <- abs_eems_parafac %>% arrange(Date,Reservoir,Station,Depth,Rep,Dilution)
abs.eems.parafac <- abs.eems.parafac %>% rename(Site=Station,DateTime=Date,Depth_m=Depth)
abs.eems.parafac <- abs.eems.parafac[,c(1,3,2,4:36)]
abs.eems.parafac <- abs.eems.parafac %>% rename(a254_m=a254,a350_m=a350,S275_295=S275.295)
abs.eems.parafac <- abs.eems.parafac %>% rename(S350_400=S350.400,Max_Fl_Ex=Max.Fl.Ex,Max_Fl_Em=Max.Fl.Em)
abs.eems.parafac <- abs.eems.parafac %>% rename(Max_Fl=Max.Fl)
abs.eems.parafac <- abs.eems.parafac %>% rename(T_B=T.B,T_M=T.M,T_N=T.N,T_C=T.C,A_T=A.T,A_C=A.C,A_M=A.M,M_C=M.C,C_N=C.N)

#Calculate mean and sd for a254, similar to QA/QC for absorbance data, based on full dataset
a254_alldata <- abs.eems.parafac %>% select(DateTime,Site,Depth_m,a254_m)
a254_mean <- a254_alldata %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(mean))
a254_sd <- a254_alldata %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(sd))

#Combine site and depth_m columns for plotting
a254_mean$Location <- paste(a254_mean$Site,a254_mean$Depth_m)
a254_sd$Location <- paste(a254_sd$Site,a254_sd$Depth_m)

#Plot by location for a254
ggplot(a254_mean,mapping=aes(x=DateTime,y=a254_m,color=Location))+geom_line()+geom_point()+geom_errorbar(a254_mean,mapping=aes(x=DateTime,y=a254_m,ymin=a254_m-a254_sd$a254_m,ymax=a254_m+a254_sd$a254_m,color=Location))+theme_classic(base_size=15)

#Calculate mean and sd for peak A
peakA <- abs.eems.parafac %>% select(DateTime,Site,Depth_m,A)
peakA_mean <- peakA %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(mean))
peakA_sd <- peakA %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(sd))

#Combine site and depth_m columns for plotting
peakA_mean$Location <- paste(peakA_mean$Site,peakA_mean$Depth_m)
peakA_sd$Location <- paste(peakA_sd$Site,peakA_sd$Depth_m)

#Plot by location for peak A
ggplot(peakA_mean,mapping=aes(x=DateTime,y=A,color=Location))+geom_line()+geom_point()+geom_errorbar(peakA_mean,mapping=aes(x=DateTime,y=A,ymin=A-peakA_sd$A,ymax=A+peakA_sd$A,color=Location))+theme_classic(base_size=15)

#Calculate mean and sd for peak T
peakT <- abs.eems.parafac %>% select(DateTime,Site,Depth_m,T)
peakT_mean <- peakT %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(mean))
peakT_sd <- peakT %>% group_by(DateTime,Site,Depth_m) %>% summarise_all(funs(sd))

#Combine site and depth_m columns for plotting
peakT_mean$Location <- paste(peakT_mean$Site,peakT_mean$Depth_m)
peakT_sd$Location <- paste(peakT_sd$Site,peakT_sd$Depth_m)

#Plot by location for peak T
ggplot(peakT_mean,mapping=aes(x=DateTime,y=T,color=Location))+geom_line()+geom_point()+geom_errorbar(peakT_mean,mapping=aes(x=DateTime,y=T,ymin=T-peakT_sd$T,ymax=T+peakT_sd$T,color=Location))+theme_classic(base_size=15)

#Plot sd for a254, peak A, and peak T
ggplot(a254_sd,mapping=aes(x=DateTime,y=a254_m,color=Location))+geom_line()+geom_point()+theme_classic(base_size=15)
ggplot(peakA_sd,mapping=aes(x=DateTime,y=A,color=Location))+geom_line()+geom_point()+theme_classic(base_size=15)
ggplot(peakT_sd,mapping=aes(x=DateTime,y=T,color=Location))+geom_line()+geom_point()+theme_classic(base_size=15)

#Export full dataset as .csv file
write.csv(abs.eems.parafac,"./Data/20210511_abs_eems_parafac.csv")

### Add in flags for absorbance and fluorescence data
#Absorbance data
# 0 = Data good!
# 1 = Absorbance not collected
# 2 = Low absorbance values; use with caution
abs_eems_parafac_flags <- abs.eems.parafac %>% 
  mutate(abs_Flag = ifelse(is.na(Sr) & is.na(a254_m), 1,
                           ifelse(is.na(Sr), 2, 
                                  0)))
#Fluorescence data
# 1 = Fluorescence not collected
abs_eems_parafac_flags <- abs_eems_parafac_flags %>% 
  mutate(fl_Flag = ifelse(is.na(Max_Fl_Ex),1,0))

#Export out file w/ flags for EDI!
write.csv(abs_eems_parafac_flags,"./EDI_2021/20210511_OpticalData.csv")
