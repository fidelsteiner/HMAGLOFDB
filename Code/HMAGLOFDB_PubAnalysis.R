################################################################################
# Analysis of all GLOF events from the HMAGLOFDB
# 
# HMAGLOFDB_Analysis.R
#
# ReadMe: 
# Read and analyze GLOF events
# 1) Read out database for basic statistics
# 2) Extract topographic variables for impact analysis
# 3) Pairing with RGI + dhdt datasets
# 4) Pairing with permafrost map
#
#
#
#
# Created:          2022/02/05
# Latest Revision:  2022/06/15
#
#
# Jakob F Steiner| ICIMOD | jakob.steiner@icimod.org | x-hydrolab.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

library('geosphere')
library(leastcostpath)
library(raster)
library(rgdal)
library(sf)
library(vioplot)
library('scales')
library(lwgeom)

################ 
# paths and raw data
################

projec_utm <- '+proj=utm +datum=WGS84'

data_path <- 'D:\\Work\\Repositories\\HMAGLOFDB\\Database\\GLOFs'       # location of the database on local station
db_file <- 'HMAGLOFDB_v0.1.csv'                                         # file name of database

figures_path <- 'D:\\Work\\ICIMODProjects\\GLOF_Database\\Figures'

DEM_path <- 'D:\\Work\\GeospatialData\\HMA\\SRTM\\VoidFilled'           # location of DEM for topographic analysis
DEM_file <- 'SRTM_Corrected_Extended_HMA.tif'                           # DEM filename

HKH_path <- 'D:\\Work\\GeospatialData\\HMA\\HKH_Boundary'               # location of HKH Outline (http://rds.icimod.org/Home/DataDetail?metadataId=3924&searchlist=True)
HKH_file <- 'HKH_Outline.shp'
HKH_outline<-readOGR(dsn=HKH_path&'\\'&HKH_file)                        # read in HKH outline
HKH_outline <- spTransform(HKH_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")                          
HKH_outlinesf <- read_sf(HKH_path&'\\'&HKH_file)
HKH_outlinesf <- st_transform(HKH_outlinesf,crs="+proj=longlat +datum=WGS84 +no_defs")

HKH2_path <- 'D:\\Work\\GeospatialData\\HMA\\Chinese_HKH_Boundary'      # location of the HKH outline used by Chinese scientists (REF!); the Hindukush outline is added from the RGI
HKH2_file <- 'Chinese_HKHOutline.shp'
HKH2_outlinesf <- read_sf(HKH2_path&'\\'&HKH2_file)
HKH2_outlinesf <- st_transform(HKH2_outlinesf,crs="+proj=longlat +datum=WGS84 +no_defs")

lakesHKH_path <- 'D:\\Work\\Repositories\\HMAGLOFDB\\Database\\Lakes\\HKHLAKEDB'   # location of the HKH lake database
lakes_file_2000 <- 'HKH_2000.shp'
lakes_file_2005 <- 'HKH_2005.shp'
lakes_file_2010 <- 'HKH_2010.shp'
lakes_file_2020 <- 'HKH_2020.shp'

lakesHiMAG_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacialLakes\\Chen2021\\Hi-MAG database'   # location of HiMAG (DOI: 10.5194/essd-13-741-2021)

lakesIGH_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacialLakes\\Wang2020\\igh_Asia_glacial_lake\\High_Asia_glacial_lake'   # location of IGH (DOI: 10.5194/essd-12-2169-2020)
lakesIGH_file_1990 <- 'High_Asia_glacial_lake_1990.shp'
lakesIGH_file_2018 <- 'High_Asia_glacial_lake_2018.shp'

#outline_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Lee2021\\SuppInfo_vectors\\SuppInfo_vectors'
#outline_file <- 'glaciers_LIA.shp'

cls <- c(Elev_lake="numeric", Elev_impact="numeric")            # Load database
db_data <-read.csv(data_path&'\\'&db_file,header = T,sep=',')

DEM_HMA <- raster(DEM_path&'\\'&DEM_file)                       # Load DEM

#ogrInfo(outline_path&'\\'&outline_file)
#LIA_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
#LIA_outline <- spTransform(LIA_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

# Load RGI
path_RGI <- 'D:\\Work\\GeospatialData\\RGI60'                                     # Folder for RGI glacier outlines
RGI15_filename <- '15_rgi60_SouthAsiaEast.shp'
RGI14_filename <- '14_rgi60_SouthAsiaWest.shp'                                    # RGI filename
RGI13_filename <- '13_rgi60_CentralAsia.shp' 

ogrInfo(path_RGI&'\\'&RGI15_filename)
RGI60_15<-read_sf(dsn=path_RGI&'\\'&RGI15_filename)


ogrInfo(path_RGI&'\\'&RGI14_filename)
RGI60_14<-read_sf(dsn=path_RGI&'\\'&RGI14_filename)
#RGI60_14 <- spTransform(RGI60_14, projec)

ogrInfo(path_RGI&'\\'&RGI13_filename)
RGI60_13<-read_sf(dsn=path_RGI&'\\'&RGI13_filename)
#RGI60_13 <- spTransform(RGI60_13, projec)

################ 
# Evaluate Lake databases
################

sf::sf_use_s2(FALSE) # turn off buffer for finding intersections

# ICIMOD database

ogrInfo(lakesHKH_path&'\\'&lakes_file_2005)
lakes2005_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2005)
lakes2005_outline <- spTransform(lakes2005_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2005_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2005)
lakes2005_sf$Area <- as.numeric(st_area(lakes2005_sf))

ogrInfo(lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_outline <- spTransform(lakes2020_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2020_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_sf <- st_transform(lakes2020_sf,crs=st_crs(lakes2005_sf))
lakes2020_sf$Area <- as.numeric(st_area(lakes2020_sf))

ogrInfo(lakesHKH_path&'\\'&lakes_file_2000)
lakes2000_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2000)
lakes2000_outline <- spTransform(lakes2000_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2000_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2000)
lakes2000_sf <- st_transform(lakes2000_sf,crs=st_crs(lakes2005_sf))
lakes2000_sf$Area <- as.numeric(st_area(lakes2000_sf))

ogrInfo(lakesHKH_path&'\\'&lakes_file_2010)
lakes2010_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2010)
lakes2010_outline <- spTransform(lakes2010_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2010_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2010)
lakes2010_sf <- st_transform(lakes2010_sf,crs=st_crs(lakes2005_sf))
lakes2010_sf$Area <- as.numeric(st_area(lakes2010_sf))

matchID <- which(!is.na(match(lakes2020_sf$GLIMS_ID,lakes2005_sf$GLIMS_ID)))
comb0520 <- rbind(lakes2020_sf[-matchID,],lakes2005_sf)
matchID2 <- which(!is.na(match(lakes2000_sf$GLIMS_ID,comb0520$GLIMS_ID)))
comb050020 <- rbind(lakes2000_sf[-matchID2,],comb0520)
matchID3 <- which(!is.na(match(lakes2010_sf$GLIMS_ID,comb050020$GLIMS_ID)))
combAll <- rbind(lakes2010_sf[-matchID3,],comb050020)

combAll$Altitude <- raster::extract(DEM_HMA, cbind(combAll$Longitude,combAll$Latitude), sp = T)[]
lakes2005_sf$Altitude <- raster::extract(DEM_HMA, cbind(lakes2005_sf$Longitude,lakes2005_sf$Latitude), sp = T)[]

# Wang database

ogrInfo(lakesIGH_path&'\\'&lakesIGH_file_1990)
lakes1990_outline <- readOGR(dsn=lakesIGH_path&'\\'&lakesIGH_file_1990)
lakes1990_outline <- spTransform(lakes1990_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes1990_sf <- read_sf(lakesIGH_path&'\\'&lakesIGH_file_1990)
lakes1990_sf <- st_transform(lakes1990_sf ,crs=st_crs(lakes2005_sf))

ogrInfo(lakesIGH_path&'\\'&lakesIGH_file_2018)
lakes2018_outline <- readOGR(dsn=lakesIGH_path&'\\'&lakesIGH_file_1990)
lakes2018_outline <- spTransform(lakes2018_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2018_sf <- read_sf(lakesIGH_path&'\\'&lakesIGH_file_2018)
lakes2018_sf <- st_transform(lakes2018_sf ,crs=st_crs(lakes2005_sf))

# combine all unique lakes from the two years for complete set
matchWang <- which(!is.na(match(lakes2018_sf$GLAKE_ID,lakes1990_sf$GLAKE_ID)))
wangComb <- rbind(lakes2018_sf[-matchWang,],lakes1990_sf)

# limit the lakes to same area that was used for the ICIMOD inventory (HKH = 2005 outline, HKH2 = 2000/2010/2020 outline)
wangComb_HKH2 <- st_intersection(wangComb,HKH2_outlinesf)
wangComb_HKH <- st_intersection(wangComb,HKH_outlinesf)

# Chen database

for(k in 2008:2017){
ogrInfo(lakesHiMAG_path&'\\Hi_MAG_database_'&k&'.shp')
lakesChen_outline <- readOGR(dsn=lakesHiMAG_path&'\\Hi_MAG_database_'&k&'.shp')
lakesChen_outline <- spTransform(lakesChen_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakesChen_sf <- read_sf(lakesHiMAG_path&'\\Hi_MAG_database_'&k&'.shp')
lakesChen_sf <- st_transform(lakesChen_sf ,crs=st_crs(lakes2005_sf))

# combine all unique lakes from the two years for complete set
if(k>2008){
matchChen <- which(!is.na(match(lakesChen_sf$GL_ID,chenComb$GL_ID)))
chenComb <- rbind(lakesChen_sf[-matchChen,],chenComb)}
if(k==2008){
  chenComb <- lakesChen_sf}
}
# limit the lakes to same area that was used for the ICIMOD inventory (HKH = 2005 outline, HKH2 = 2000/2010/2020 outline)
chenComb_HKH2 <- st_intersection(chenComb,HKH2_outlinesf)
chenComb_HKH <- st_intersection(chenComb,HKH_outlinesf)

# Extract GLOFs that occured within the HKH arch
db_data$Lon_lake[which(is.na(db_data$Lon_lake))] <- 0
db_data$Lat_lake[which(is.na(db_data$Lat_lake))] <- 0
GLOFLoc <- data.frame(Longitude = db_data$Lon_lake,
                  Latitude = db_data$Lat_lake)
coordinates(GLOFLoc) <- ~ Longitude + Latitude
proj4string(GLOFLoc) <- proj4string(HKH_outline)
GLOFwithin <- which(!is.na(over(GLOFLoc,HKH_outline)[,1]))

# Plot ELevation distributions in vertical histogram

histList <- c(combAll$Altitude,wangComb_HKH$GL_Elev,chenComb_HKH$GL_Elev)
Elevbreaks <- c(1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500)

png(file=figures_path&'\\ElevDistribution.png', res = 180,width=1600,height=900)
par(mar=c(2,4.1,2,0.2),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,3,4), nrow = 1, ncol = 4, byrow = TRUE))

xhist <- hist(lakes2005_sf$Altitude, breaks = Elevbreaks, plot = FALSE)
barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "", ylab="Elevation [m a.s.l.]",xpd = FALSE,xlim = c(0,9000),col='black')
axis(2,at=seq(1,10,by=1),labels=Elevbreaks[1:10])
grid(nx=NULL, ny=NULL)
xhist <- hist(wangComb_HKH$GL_Elev, breaks = Elevbreaks, plot = FALSE)
barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "", ylab="",xpd = FALSE,xlim = c(0,9000),col='red')
#axis(2,at=seq(1,9,by=1),labels=Elevbreaks[1:9])
grid(nx=NULL, ny=NULL)
xhist <- hist(chenComb_HKH$GL_Elev, breaks = Elevbreaks, plot = FALSE)
barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "", ylab="",xpd = FALSE,xlim = c(0,9000),col='blue')
grid(nx=NULL, ny=NULL)
xhist <- hist(db_data$Elev_lake[GLOFwithin], breaks = Elevbreaks, plot = FALSE)
barplot(xhist$counts, axes = TRUE, space = 0, horiz=TRUE, xlab= "", ylab="",xpd = FALSE,xlim = c(0,150),col='grey')
grid(nx=NULL, ny=NULL)
dev.off()

################ 
# Cumulative distribution of lake areas

Area_cdf_ICIMOD <- ecdf(lakes2005_sf$Area)
Area_cdf_Wang <- ecdf(wangComb_HKH$GL_Area)
Area_cdf_Chen <- ecdf(chenComb_HKH$GL_Area)

png(file=figures_path&'\\ECDF_lakedatabases.png', res = 160,width=1800,height=900)
par(mar=c(5,7,0,1),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(3,1))
layout(matrix(c(1,1,1,2), nrow = 4, ncol = 1, byrow = TRUE))
plot(sort(lakes2005_sf$Area),Area_cdf_ICIMOD(sort(lakes2005_sf$Area)),type='l',log="x",xlim=c(2000, 5000000),lwd=2, main='', panel.first=grid(equilogs=T),xaxt="n",xlab ='',ylab = 'Probability of exceedance [-]',do.points=F,verticals = T,col.01line = NULL,ylim=c(0,1))
points(sort(wangComb_HKH$GL_Area),Area_cdf_Wang(sort(wangComb_HKH$GL_Area)),type='l',add=T,col='red',xaxt="n",lwd=2,col.01line = NULL)
points(sort(chenComb_HKH$GL_Area),Area_cdf_Chen(sort(chenComb_HKH$GL_Area)),type='l',add=T,col='blue',xaxt="n",lwd=2,col.01line = NULL)
legend('topleft',bty='n',col=c('black','red','blue','grey'),lty = 1,lwd = 2,legend=c('ICIMOD inventory', 'Wang et al., 2020','Chen et al., 2021','recorded GLOFs'),cex=2)

uniqueGLOFs <- match(unique(db_data$GL_ID[GLOFwithin]),db_data$GL_ID[GLOFwithin])
uniqueGLOFs2 <- match(unique(db_data$GL_ID),db_data$GL_ID)

GLOFArea <- combAll$Area[match(db_data$GL_ID[GLOFwithin],combAll$GLIMS_ID)[uniqueGLOFs]]/ 10^6
GLOFArea[is.na(GLOFArea)] <- wangComb$GL_Area[match(db_data$GL_ID[GLOFwithin],wangComb$GLAKE_ID)][uniqueGLOFs][is.na(GLOFArea)] / 10^6
GLOFArea[is.na(GLOFArea)] <- chenComb$GL_Area[match(db_data$GL_ID[GLOFwithin],chenComb$GL_ID)][uniqueGLOFs][is.na(GLOFArea)] / 10^6

GLOFArea2 <- combAll$Area[match(db_data$GL_ID,combAll$GLIMS_ID)]/ 10^6
GLOFArea2[is.na(GLOFArea2)] <- wangComb$GL_Area[match(db_data$GL_ID,wangComb$GLAKE_ID)][is.na(GLOFArea2)] / 10^6
GLOFArea2[is.na(GLOFArea2)] <- chenComb$GL_Area[match(db_data$GL_ID,chenComb$GL_ID)][is.na(GLOFArea2)] / 10^6

GLOFArea3 <- combAll$Area[match(db_data$GL_ID,combAll$GLIMS_ID)[uniqueGLOFs2]]/ 10^6
GLOFArea3[is.na(GLOFArea3)] <- wangComb$GL_Area[match(db_data$GL_ID,wangComb$GLAKE_ID)][uniqueGLOFs2][is.na(GLOFArea3)] / 10^6
GLOFArea3[is.na(GLOFArea3)] <- chenComb$GL_Area[match(db_data$GL_ID,chenComb$GL_ID)][uniqueGLOFs2][is.na(GLOFArea3)] / 10^6

boxplot(GLOFArea*10^6,horizontal=TRUE,log="x",ylim=c(2000, 5e6),xaxt="n",xlab = expression(paste("Area [",m^2,"]", sep="")))
myjitter<-jitter(rep(1, length(GLOFArea*10^6)), amount=0.2)
points(GLOFArea*10^6,myjitter,  pch=20, col='grey',cex=1.5)

at.y <- outer(1:9, 10^(0:8))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(1, at=at.y, labels=lab.y, las=1)
dev.off()

################ 
# Extract the database
################

missingYear <- length(which(is.na(db_data$Year_exact))) / length(db_data$Year_exact) * 100   # Events where we do not know the year
availableDay <- length(which(!is.na(db_data$Day))) / length(db_data$Year_exact) * 100

# Distribution by country and RGI region

cntrDistrib <- c(length(which(db_data$Country=='Afghanistan')),length(which(db_data$Country=='Kyrgyzstan')),length(which(db_data$Country=='Tajikistan')),
                 length(which(db_data$Country=='Kazakhstan')),length(which(db_data$Country=='China')),length(which(db_data$Country=='India')),length(which(db_data$Country=='Nepal')),
                 length(which(db_data$Country=='Bhutan')),length(which(db_data$Country=='Pakistan'))) / length(db_data$Year_exact) * 100

regDistrib <- c(length(which(db_data$Region_RGI=='13-01')),length(which(db_data$Region_RGI=='13-02')),length(which(db_data$Region_RGI=='13-02')),
                 length(which(db_data$Region_RGI=='13-03')),length(which(db_data$Region_RGI=='13-04')),length(which(db_data$Region_RGI=='14-01')),length(which(db_data$Region_RGI=='14-02')),
                 length(which(db_data$Region_RGI=='14-03')),length(which(db_data$Region_RGI=='15-01')), length(which(db_data$Region_RGI=='15-02')), length(which(db_data$Region_RGI=='15-03'))) / length(db_data$Year_exact) * 100

regDistrib_HiMAP <- as.data.frame(table(db_data$Region_HiMAP))
regDistrib_HiMAP$Freq/length(db_data$Year_exact) * 100

# Elevation and drop height of GLOFs
meanElev_GLOF <- c(min(db_data$Elev_lake,na.rm=T),mean(db_data$Elev_lake,na.rm=T),max(db_data$Elev_lake,na.rm=T))
meandH_GLOF <- c(min(db_data$Elev_lake - db_data$Elev_impact,na.rm=T),mean(db_data$Elev_lake - db_data$Elev_impact,na.rm=T),max(db_data$Elev_lake - db_data$Elev_impact,na.rm=T))

# Data on lake formation/mechanisms
lakeFormatinRecorded <- length(which(db_data$Driver_lake!='unknown')) / length(db_data$Year_exact) * 100 
GLOFFormationRecorded <- length(which(db_data$Driver_GLOF!='unknown')) / length(db_data$Year_exact) * 100 
GLOFMechanismRecorded <- length(which(db_data$Mechanism!='unknown')) / length(db_data$Year_exact) * 100

# Data on discharge/volume
volRecorded <- length(which(!is.na(as.numeric(db_data$Volume)))) / length(db_data$Year_exact) * 100 
QRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_water)))) / length(db_data$Year_exact) * 100 
QsRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_solid)))) / length(db_data$Year_exact) * 100 

################ 
# Extract topographic data
################


png(file=figures_path&'\\GLOFLakesCharacteristics.png', res = 160,width=1800,height=900)
par(mar=c(5,7,0,1),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
# Area
Areadf <- data.frame(y=c(wangComb$GL_Area/10^6,GLOFArea3),x=c(rep(1,length(wangComb$GL_Area)),rep(2,length(GLOFArea3))))
with(Areadf, boxplot(y~x,log="y",ylim=c(0.001,10),border=c("blue", "red"),xaxt='n',xlab = '',yaxt ='n',ylab=expression(paste("Area [",km^2,"]", sep=""))))
at.y <- outer(1:8, 10^(-3:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1)

myjitter<-jitter(rep(1, length(wangComb$GL_Area/10^6)), amount=0.2)
points(myjitter,wangComb$GL_Area/10^6,  pch=1, col = alpha('blue', 0.006),cex=1.5)
myjitter<-jitter(rep(2, length(GLOFArea3)), amount=0.2)
points(myjitter,GLOFArea3,  pch=1, col= alpha('red', 0.2),cex=1.5)

# Area change of individual lakes - identify lakes in individual inventories that appear in more than one year

ICIMOD0500match <- which(!is.na(match(lakes2005_sf$GLIMS_ID,lakes2000_sf$GLIMS_ID)))
ICIMOD0510match <- which(!is.na(match(lakes2005_sf$GLIMS_ID,lakes2010_sf$GLIMS_ID)))
ICIMOD1020match <- which(!is.na(match(lakes2010_sf$GLIMS_ID,lakes2020_sf$GLIMS_ID)))

dA0500_ICIMOD <- (lakes2005_sf$Area[ICIMOD0500match] - lakes2000_sf$Area[match(lakes2005_sf$GLIMS_ID,lakes2000_sf$GLIMS_ID)[ICIMOD0500match]]) / 5
dA0510_ICIMOD <- (lakes2010_sf$Area[match(lakes2005_sf$GLIMS_ID,lakes2010_sf$GLIMS_ID)[ICIMOD0510match]] - lakes2005_sf$Area[ICIMOD0510match]) / 5
dA1020_ICIMOD <- (lakes2020_sf$Area[match(lakes2010_sf$GLIMS_ID,lakes2020_sf$GLIMS_ID)[ICIMOD1020match]] - lakes2010_sf$Area[ICIMOD1020match]) / 10

ID0500 <- cbind(lakes2005_sf$GLIMS_ID[ICIMOD0500match],dA0500_ICIMOD)
ID0510 <- cbind(lakes2005_sf$GLIMS_ID[ICIMOD0510match],dA0510_ICIMOD)
ID1020 <- cbind(lakes2010_sf$GLIMS_ID[ICIMOD1020match],dA1020_ICIMOD)

IDall <- cbind(ID0500,ID0510[match(ID0500[,1],ID0510[,1]),2],ID1020[match(ID0500[,1],ID1020[,1]),2])
IDall <- cbind(IDall, rowMeans(cbind(as.numeric(IDall[,2]),as.numeric(IDall[,3]),as.numeric(IDall[,4])),na.rm=T))
#IDall[which(as.numeric(IDall[,5])==0),5] <- NA

Wang19902018match <- which(!is.na(match(lakes1990_sf$GLAKE_ID,lakes2018_sf$GLAKE_ID)))
dA9018_Wang <- (lakes2018_sf$GL_Area[match(lakes1990_sf$GLAKE_ID,lakes2018_sf$GLAKE_ID)[Wang19902018match]] - lakes1990_sf$GL_Area[Wang19902018match]) / 28


lakesChen_sf1 <- read_sf(lakesHiMAG_path&'\\Hi_MAG_database_2008.shp')
lakesChen_sf1 <- st_transform(lakesChen_sf1 ,crs=st_crs(lakes2005_sf))
dAChen <- lakesChen_sf1$GL_ID

for(k in 2009:2017){
    lakesChen_sf_prev <- read_sf(lakesHiMAG_path&'\\Hi_MAG_database_'&k-1&'.shp')
    lakesChen_sf_prev <- st_transform(lakesChen_sf_prev ,crs=st_crs(lakes2005_sf))
    
    matchChen <- which(!is.na(match(lakesChen_sf$GL_ID,lakesChen_sf_prev$GL_ID)))
    matchInit <- which(!is.na(match(lakesChen_sf$GL_ID[matchChen],lakesChen_sf1$GL_ID)))
    dAChen <- cbind(dAChen,(lakesChen_sf$GL_Area[matchChen][matchInit] - lakesChen_sf_prev$GL_Area[match(lakesChen_sf$GL_ID,lakesChen_sf_prev$GL_ID)[matchChen][matchInit]]))
}

ChenMeans <- rowMeans(cbind(as.numeric(dAChen[,2]),as.numeric(dAChen[,3]),as.numeric(dAChen[,4]),as.numeric(dAChen[,5]),as.numeric(dAChen[,6]),as.numeric(dAChen[,7]),as.numeric(dAChen[,8]),as.numeric(dAChen[,9]),as.numeric(dAChen[,10])),na.rm=T)
ChenGLOFmatch <-  ChenMeans[which(!is.na(match(dAChen[,1],db_data$GL_ID)))]


# Get area change of GLOF lakes from the ICIMOD DB
ICIMOD0500Tab <- cbind(lakes2005_sf$GLIMS_ID[ICIMOD0500match],dA0500_ICIMOD)
ICIMOD0510Tab <- cbind(lakes2005_sf$GLIMS_ID[ICIMOD0510match],dA0510_ICIMOD)
ICIMOD1020Tab <- cbind(lakes2010_sf$GLIMS_ID[ICIMOD1020match],dA1020_ICIMOD)

db_data$dA0500 <- ICIMOD0500Tab[match(db_data$GL_ID,ICIMOD0500Tab[,1]),2]
db_data$dA0510 <- ICIMOD0510Tab[match(db_data$GL_ID,ICIMOD0510Tab[,1]),2]
db_data$dA1020 <- ICIMOD1020Tab[match(db_data$GL_ID,ICIMOD1020Tab[,1]),2]

GLOFICIMODdA <- rowMeans(cbind(as.numeric(db_data$dA0500),as.numeric(db_data$dA0510),as.numeric(db_data$dA1020)),na.rm=T)

# Get area change of GLOF lakes from the Wang DB
WangGLOFmatch <- which(!is.na(match(lakes1990_sf$GLAKE_ID[Wang19902018match],db_data$GL_ID)))
GLOFWangdA <- dA9018_Wang[WangGLOFmatch]

dAreadf <- data.frame(y=c(as.numeric(IDall[,5]),dA9018_Wang,ChenMeans,GLOFICIMODdA,GLOFWangdA,ChenGLOFmatch),x=c(rep(1,length(as.numeric(IDall[,5]))),rep(2,length(dA9018_Wang)),rep(3,length(ChenMeans)),rep(4,length(GLOFICIMODdA)),rep(5,length(GLOFWangdA)),rep(6,length(ChenGLOFmatch))))
with(dAreadf, boxplot(y~x,ylim=c(-2000,2000),border=c("blue","blue","blue", "red", "red","red"),xaxt='n',xlab = '',ylab=c(expression('Area change [m'^"2"*' yr'^"-1"*']'))))
grid(nx=NULL, ny=NULL)
#at.y <- outer(1:8, 10^(-3:4))
#lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
#axis(2, at=at.y, labels=lab.y, las=1)

dev.off()


# distance to glacier

RGI60_13_UTM <- st_transform(RGI60_13, crs=projec_utm)
RGI13BB <- st_bbox(RGI60_13_UTM)
RGI60_14_UTM <- st_transform(RGI60_14, crs=projec_utm)
RGI14BB <- st_bbox(RGI60_14_UTM)
RGI60_15_UTM <- st_transform(RGI60_15, crs=projec_utm)
RGI15BB <- st_bbox(RGI60_15_UTM)
wangComb_UTM <- st_transform(wangComb, crs=projec_utm)
combAll_UTM <- st_transform(combAll, crs=projec_utm)

lakeglacierDis <- vector()
for(i in 1: dim(wangComb_UTM)[1]){
  lakeBB <- st_bbox(wangComb_UTM$geometry[i])
  

  if(as.numeric(lakeBB[1])>as.numeric(RGI13BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI13BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI13BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI13BB[4])){
    
    cropRGI <- st_crop(RGI60_13_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis[i] <- NA
    }else{
      lakeglacierDis[i] <- min(st_distance(wangComb_UTM$geometry[i],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI14BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI14BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI14BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI14BB[4])){
    cropRGI <- st_crop(RGI60_14_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis[i] <- NA
    }else{
      lakeglacierDis[i] <- min(st_distance(wangComb_UTM$geometry[i],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI15BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI15BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI15BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI15BB[4])){
    cropRGI <- st_crop(RGI60_15_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis[i] <- NA
    }else{
      lakeglacierDis[i] <- min(st_distance(wangComb_UTM$geometry[i],cropRGI))
    }
  }else{
    lakeglacierDis[i] <- NA
  }
}

lakeglacierDis_ICIMOD <- vector()
subSetICIMOD <- which(is.na(match(combAll$GLIMS_ID,wangComb$GLAKE_ID)))

for(i in 1: length(subSetICIMOD)){
  lakeBB <- st_bbox(combAll_UTM$geometry[subSetICIMOD[i]])
  
  
  if(as.numeric(lakeBB[1])>as.numeric(RGI13BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI13BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI13BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI13BB[4])){
    
    cropRGI <- st_crop(RGI60_13_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_ICIMOD[i] <- NA
    }else{
      lakeglacierDis_ICIMOD[i] <- min(st_distance(combAll_UTM$geometry[subSetICIMOD[i]],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI14BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI14BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI14BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI14BB[4])){
    cropRGI <- st_crop(RGI60_14_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_ICIMOD[i] <- NA
    }else{
      lakeglacierDis_ICIMOD[i] <- min(st_distance(combAll_UTM$geometry[subSetICIMOD[i]],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI15BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI15BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI15BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI15BB[4])){
    cropRGI <- st_crop(RGI60_15_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_ICIMOD[i] <- NA
    }else{
      lakeglacierDis_ICIMOD[i] <- min(st_distance(combAll_UTM$geometry[subSetICIMOD[i]],cropRGI))
    }
  }else{
    lakeglacierDis_ICIMOD[i] <- NA
  }
}

lakeglacierDis_Chen <- vector()
subSetChen <- which(is.na(match(chenComb$GL_ID,c(wangComb$GLAKE_ID,combAll$GLIMS_ID))))

for(i in 1: length(subSetChen)){
  lakeBB <- st_bbox(chenComb$geometry[subSetChen[i]])
  
  
  if(as.numeric(lakeBB[1])>as.numeric(RGI13BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI13BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI13BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI13BB[4])){
    
    cropRGI <- st_crop(RGI60_13_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_Chen[i] <- NA
    }else{
      lakeglacierDis_Chen[i] <- min(st_distance(chenComb$geometry[subSetChen[i]],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI14BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI14BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI14BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI14BB[4])){
    cropRGI <- st_crop(RGI60_14_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_Chen[i] <- NA
    }else{
      lakeglacierDis_Chen[i] <- min(st_distance(chenComb$geometry[subSetChen[i]],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI15BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI15BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI15BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI15BB[4])){
    cropRGI <- st_crop(RGI60_15_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_Chen[i] <- NA
    }else{
      lakeglacierDis_Chen[i] <- min(st_distance(chenComb$geometry[subSetChen[i]],cropRGI))
    }
  }else{
    lakeglacierDis_Chen[i] <- NA
  }
}

dbLat <- db_data$Lat_lake
dbLat[duplicated(dbLat)] <- NA
uniquedb <- which(!is.na(dbLat))

for(i in 1: length(uniquedb)){
  
  cord.dec <- SpatialPoints(cbind(db_data$Lon_lake[uniquedb[i]],db_data$Lat_lake[uniquedb[i]]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  cord.UTM <- spTransform(cord.dec, CRS(projec_utm))
  
  lakeBB <- st_bbox(cord.UTM)
  
  if(as.numeric(lakeBB[1])>as.numeric(RGI13BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI13BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI13BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI13BB[4])){
    
    cropRGI <- st_crop(RGI60_13_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_DB[i] <- NA
    }else{
      lakeglacierDis_DB[i] <- geosphere::dist2Line(p = st_bbox(cord.dec), line = st_transform(cropRGI,crs="+proj=longlat +datum=WGS84 +no_defs"))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI14BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI14BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI14BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI14BB[4])){
    cropRGI <- st_crop(RGI60_14_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_Chen[i] <- NA
    }else{
      lakeglacierDis_Chen[i] <- min(st_distance(chenComb$geometry[subSetChen[i]],cropRGI))
    }
  }else if(as.numeric(lakeBB[1])>as.numeric(RGI15BB[1])&&as.numeric(lakeBB[2])>as.numeric(RGI15BB[2])&&as.numeric(lakeBB[3])<as.numeric(RGI15BB[3])&&as.numeric(lakeBB[4])<as.numeric(RGI15BB[4])){
    cropRGI <- st_crop(RGI60_15_UTM,lakeBB+c(-5000,-5000,5000,5000))
    if(dim(cropRGI)[1]==0){
      lakeglacierDis_Chen[i] <- NA
    }else{
      lakeglacierDis_Chen[i] <- min(st_distance(chenComb$geometry[subSetChen[i]],cropRGI))
    }
  }else{
    lakeglacierDis_Chen[i] <- NA
  } 
  
}

#par(mfrow=c(4,4),mar=c(1,2,1,3),mai = c(0.1, 0.1, 0.1, 0.1),cex.lab=1,cex.axis=1)  

flowpathL <- vector() # initialize empty vector for flow length of the GLOF
flowPathH <- vector() # initialize empty vector for drop height
flowPathalpha_mean <- vector()  # initialize empty vector for mean slope
flowPathalpha_max <- vector()   # initialize empty vector for max slope


for(k in 245:length(db_data$GF_ID)){
  if(!is.na(db_data$Lon_impact[k]) == T){
idGF <- k

extent_case <- c(min(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),max(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),
            min(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF]),max(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF])) +
  c(-0.03,0.03,-0.03,0.03)
  

DEM_case <- crop(DEM_HMA,extent_case)
LIA_case <- crop(LIA_outline,extent_case)

#r.fd <- terrain(DEM_case, opt='flowdir')
#r.p <- flowPath(r.fd, c(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]))

#p.xy <- xyFromCell(r.fd,r.p)
#plot(DEM_case)
#points(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF],col='blue',pch=4)
#points(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF],col='blue',pch=4)
#lines(p.xy,col='red', lwd=4)

stpnt <- sp::SpatialPoints(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]))
enpnt <- sp::SpatialPoints(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF]))

rast_slope <- create_slope_cs(DEM_case, cost_function = 'tobler')


#plot(raster(rast_slope_barrier), col = grey.colors(100))


lcp <- create_lcp(rast_slope, stpnt, enpnt, TRUE)

#plot(DEM_case)
#if(!is.null(LIA_case)){
#plot(LIA_case, add=T,col='red')}
#lines(lcp)
#points(stpnt)
#points(enpnt)



flowpathL[k] <- lengthLine(lcp)
flowPathH[k] <- as.numeric(db_data$Elev_lake[idGF]) - as.numeric(db_data$Elev_impact[idGF])
flowPathalpha_mean[k] <- mean(extract(terrain(DEM_case,'slope','degrees'),lcp)[[1]],na.rm=T)
flowPathalpha_max[k] <- max(extract(terrain(DEM_case,'slope','degrees'),lcp)[[1]],na.rm=T)
}
}