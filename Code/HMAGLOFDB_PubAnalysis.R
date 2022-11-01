################################################################################
# Analysis of all GLOF events from the HMAGLOFDB
# 
# HMAGLOFDB_Analysis.R
#
# ReadMe: 
# Read and analyze GLOF events
# 1) Evaluate the lake databases
# 2) Extract the GLOF database
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
db_file <- 'HMAGLOFDB_v1.0.csv'                                         # file name of database

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

outline_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Lee2021\\SuppInfo_vectors\\SuppInfo_vectors'
outline_file <- 'glaciers_LIA.shp'

cls <- c(Elev_lake="numeric", Elev_impact="numeric")            # Load database
db_data <-read.csv(data_path&'\\'&db_file,header = T,sep=',')

DEM_HMA <- raster(DEM_path&'\\'&DEM_file)                       # Load DEM

ogrInfo(outline_path&'\\'&outline_file)
LIA_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
LIA_outline <- spTransform(LIA_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

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

# LOad permafrost data
path_pf <- 'D:\\Work\\GeospatialData\\HMA\\Permafrost' 
pf_filename <- 'PZI_HMA.tif'                                          # PZI filename

regionalPF <- raster(path_pf&'\\'&pf_filename)

################ 
# 1 Evaluate lake databases
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
lakes1990_sf <- lakes1990_sf[-which(lakes1990_sf$GL_Lati>46),]

ogrInfo(lakesIGH_path&'\\'&lakesIGH_file_2018)
lakes2018_outline <- readOGR(dsn=lakesIGH_path&'\\'&lakesIGH_file_1990)
lakes2018_outline <- spTransform(lakes2018_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2018_sf <- read_sf(lakesIGH_path&'\\'&lakesIGH_file_2018)
lakes2018_sf <- st_transform(lakes2018_sf ,crs=st_crs(lakes2005_sf))
lakes2018_sf <- lakes2018_sf[-which(lakes2018_sf$GL_Lati>46),]

# combine all unique lakes from the two years for complete set
matchWang <- which(!is.na(match(lakes2018_sf$GLAKE_ID,lakes1990_sf$GLAKE_ID)))
wangComb <- rbind(lakes2018_sf[-matchWang,],lakes1990_sf)

# Limit the extent to HMA
wangComb<- wangComb[-which(wangComb$GL_Lati>46),]

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
legend('topleft',bty='n',col=c('black','red','blue','grey'),lty = 1,lwd = 2,legend=c('ICIMODDB', 'WangDB','ChenDB','recorded GLOFs'),cex=2)

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
# 2 Extract the GLOF database
################

missingYear <- length(which(is.na(db_data$Year_exact))) / length(db_data$Year_exact) * 100   # Events where we do not know the year
availableDay <- length(which(!is.na(db_data$Day))) / length(db_data$Year_exact) * 100
availableMonth <- length(which(!is.na(db_data$Month))) / length(db_data$Year_exact) * 100

summerGLOFs <- length(which(db_data$Month>=6&db_data$Month<9)) / length(which(!is.na(db_data$Month))) * 100
winterGLOFs <- length(which(db_data$Month<4|db_data$Month>=11)) / length(which(!is.na(db_data$Month))) * 100

# Time series of GLOFs

db_data$Region_RGI2 <- db_data$Region_RGI # Group Regions for better plotting (RGI)
db_data$Region_RGI2[db_data$Region_RGI == '13-01'] <- 'Pamir/Alay'
db_data$Region_RGI2[db_data$Region_RGI == '13-02'] <- 'Pamir/Alay'
db_data$Region_RGI2[db_data$Region_RGI == '13-03'] <- 'Tien Shan'
db_data$Region_RGI2[db_data$Region_RGI == '13-04'] <- 'Tien Shan'
db_data$Region_RGI2[db_data$Region_RGI == '13-08'] <- 'Tibet'
db_data$Region_RGI2[db_data$Region_RGI == '13-09'] <- 'Tibet'
db_data$Region_RGI2[db_data$Region_RGI == '14-01'] <- 'Hindukush'
db_data$Region_RGI2[db_data$Region_RGI == '14-02'] <- 'Karakoram'
db_data$Region_RGI2[db_data$Region_RGI == '14-03'] <- 'Himalaya West/Central'
db_data$Region_RGI2[db_data$Region_RGI == '15-01'] <- 'Himalaya West/Central'
db_data$Region_RGI2[db_data$Region_RGI == '15-02'] <- 'Himalaya East/Hengduan Shan'
db_data$Region_RGI2[db_data$Region_RGI == '15-03'] <- 'Himalaya East/Hengduan Shan'

db_data$Region_HiMAP2 <- db_data$Region_HiMAP # Group Regions for better plotting (HiMAP)
db_data$Region_HiMAP2[db_data$Region_HiMAP == '8'] <- 'Hindukush'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '10'] <- 'Himalaya West/Central'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '12'] <- 'Himalaya East/Hengduan Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '11'] <- 'Himalaya West/Central'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '9'] <- 'Karakoram'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '6'] <- 'Pamir'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '5'] <- 'Pamir'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '7'] <- 'Pamir'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '4'] <- 'Tien Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '1'] <- 'Tien Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '2'] <- 'Tien Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '3'] <- 'Tien Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '13'] <- 'Karakoram'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '14'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '18'] <- 'Himalaya East/Hengduan Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '17'] <- 'Himalaya East/Hengduan Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '20'] <- 'Himalaya East/Hengduan Shan'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '16'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '19'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '21'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '22'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '15'] <- 'Tibet'
db_data$Region_HiMAP2[db_data$Region_HiMAP == '10'] <- 'Himalaya East/Hengduan Shan'

db_data$decade <- (c(db_data$Year_exact) %/% 10) * 10

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

seasonalCol <- RColorBrewer::brewer.pal(8, "RdYlGn")

pdecadal <- ggplot(subset(db_data, Year_exact>=1830), aes(x=factor(decade),fill = factor(Region_HiMAP2)))+
  geom_bar(stat="count", width=1, position = position_dodge2(preserve = "single"))+
  scale_fill_manual(values = seasonalCol)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"),
        legend.position = c(0.15, 0.75),
        text = element_text(size = 15))+
  xlab("") + 
  ylab("")+
  labs(fill = "")

png(file=figures_path&'\\GLOFtemporal.png', res = 160,width=1800,height=900)
print(pdecadal)
dev.off()

# Seasonal distribution of GLOFs

pseasonal <- ggplot(subset(db_data, !is.na(Month) & Month>3 & Month <11), aes(x=factor(Month),fill = factor(Region_HiMAP2)))+
  geom_bar(stat="count", width=1, position = position_dodge2(preserve = "single"))+
  scale_fill_manual(values = seasonalCol)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"),
        legend.position =  "none",
        text = element_text(size = 15))+
  scale_x_discrete("Month") +
  xlab("") + 
  ylab("")+
  labs(fill = "")

png(file=figures_path&'\\GLOFseasonal.png', res = 160,width=900,height=900)
print(pseasonal)
dev.off()


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

# Lake type
MD <- length(which(db_data$Lake_type=='Moraine dammed')) / length(db_data$Year_exact) * 100
ID <- length(which(db_data$Lake_type=='Ice dammed')) / length(db_data$Year_exact) * 100
SGL <- length(which(db_data$Lake_type=='Supraglacial')) / length(db_data$Year_exact) * 100
EGL <- length(which(db_data$Lake_type=='Englacial')) / length(db_data$Year_exact) * 100
BGL <- length(which(db_data$Lake_type=='Bedrock')) / length(db_data$Year_exact) * 100
GOF <- length(which(db_data$Lake_type=='Glacier outburst flood')) / length(db_data$Year_exact) * 100

# Data on lake formation/mechanisms
lakeFormationRecorded <- length(which(db_data$Driver_lake!='Unknown')) / length(db_data$Year_exact) * 100 
GLOFFormationRecorded <- length(which(db_data$Driver_GLOF!='Unknown')) / length(db_data$Year_exact) * 100 
GLOFMechanismRecorded <- length(which(db_data$Mechanism!='Unknown')) / length(db_data$Year_exact) * 100

# Lake formation
dataLF <- data.frame(
  category = unique(subset(db_data, Driver_lake!='Unknown')$Driver_lake),
  count=c(length(which(subset(db_data, Driver_lake!='Unknown')$Driver_lake==unique(subset(db_data, Driver_lake!='Unknown')$Driver_lake)[1])), 
          length(which(subset(db_data, Driver_lake!='Unknown')$Driver_lake==unique(subset(db_data, Driver_lake!='Unknown')$Driver_lake)[2])), 
          length(which(subset(db_data, Driver_lake!='Unknown')$Driver_lake==unique(subset(db_data, Driver_lake!='Unknown')$Driver_lake)[3])))
)

dataLF$fraction <- dataLF$count / sum(dataLF$count)
dataLF$ymax <- cumsum(dataLF$fraction)
dataLF$ymin <- c(0, head(dataLF$ymax, n=-1))
dataLF$labelPosition <- (dataLF$ymax + dataLF$ymin) / 2
dataLF$label <- paste0(dataLF$category, ": ", round(dataLF$fraction*100),"%")

pLF <- ggplot(dataLF, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = 'A'), x=2,y=500,size=10)

# GLOF formation

db_data$Driver_GLOF2 <- db_data$Driver_GLOF
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Ice avalanche"] <- 'Mass movement'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "High temperatures; intense rainfall"] <- 'T/P+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Intense rainfall"] <- 'P+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Rockfall"] <- 'Mass movement'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Upstream flood"] <- 'Mass movement'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Landslide" ] <- 'Mass movement'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Intense rainfall; high temperatures"] <- 'T/P+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "High temperatures"] <- 'T+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Debris flow" ] <- 'Mass movement'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Water level rise" ] <- 'WT+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "Water level rise; water temperature rise" ] <- 'WT+'
db_data$Driver_GLOF2[db_data$Driver_GLOF == "High temperatures; intense rainfall; landslide" ] <- 'T/P+'

dataGF <- data.frame(
  category = unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2),
  count=c(length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[1])), 
          length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[2])), 
          length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[3])),
          length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[4])), 
          length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[5])), 
          length(which(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2==unique(subset(db_data, Driver_GLOF2!='Unknown')$Driver_GLOF2)[6])))
)

dataGF$fraction <- dataGF$count / sum(dataGF$count)
dataGF$ymax <- cumsum(dataGF$fraction)
dataGF$ymin <- c(0, head(dataGF$ymax, n=-1))
dataGF$labelPosition <- (dataGF$ymax + dataGF$ymin) / 2
dataGF$label <- paste0(dataGF$category, ": ", round(dataGF$fraction*100),"%")

pGF <- ggplot(dataGF, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = 'B'), x=2,y=500,size=10)

# GLOF mechanism

db_data$Mechanism2 <- db_data$Mechanism
db_data$Mechanism2[db_data$Mechanism == "Dam piping"] <- 'Piping/seepage'
db_data$Mechanism2[db_data$Mechanism == "Dam seepage"] <- 'Piping/seepage'
db_data$Mechanism2[db_data$Mechanism == "Supraglacial lake drainage"] <- 'Others'
db_data$Mechanism2[db_data$Mechanism == "Overtopping"] <- 'Others'
db_data$Mechanism2[db_data$Mechanism == "Ice core thawing"] <- 'Moraine collapse'

dataGM <- data.frame(
  category = unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2),
  count=c(length(which(subset(db_data, Mechanism2!='Unknown')$Mechanism2==unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2)[1])), 
          length(which(subset(db_data, Mechanism2!='Unknown')$Mechanism2==unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2)[2])), 
          length(which(subset(db_data, Mechanism2!='Unknown')$Mechanism2==unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2)[3])),
          length(which(subset(db_data, Mechanism2!='Unknown')$Mechanism2==unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2)[4])), 
          length(which(subset(db_data, Mechanism2!='Unknown')$Mechanism2==unique(subset(db_data, Mechanism2!='Unknown')$Mechanism2)[5])))
  )

dataGM$fraction <- dataGM$count / sum(dataGM$count)
dataGM$ymax <- cumsum(dataGM$fraction)
dataGM$ymin <- c(0, head(dataGM$ymax, n=-1))
dataGM$labelPosition <- (dataGM$ymax + dataGM$ymin) / 2
dataGM$label <- paste0(dataGM$category, ": ", round(dataGM$fraction*100),"%")

pGM <- ggplot(dataGM, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = 'C'), x=2,y=500,size=10)

figureProc <- ggarrange(pLF,pGF,pGM,
                    ncol = 3, nrow = 1)

png(file=figures_path&'\\GLOFprocesses.png', res = 160,width=1800,height=600)
print(figureProc)
dev.off()

# Data on discharge/volume
volRecorded <- length(which(!is.na(as.numeric(db_data$Volume)))) / length(db_data$Year_exact) * 100 
QRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_water)))) / length(db_data$Year_exact) * 100 
QsRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_solid)))) / length(db_data$Year_exact) * 100 


length(which(is.na(db_data$GL_ID)))
length(which(db_data$GL_ID=='Ephemeral'))
length(which(db_data$GL_ID=='Not mapped'))
length(which(db_data$GL_ID=='No lake'))

################ 
# Extract topographic data
################
              
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
              
              
lakeglacierDis[match(db_data$GL_ID,wangComb$GLAKE_ID)]
write.csv(lakeglacierDis,file=figures_path&'\\distancesGLOFWANG.csv',row.names=F)
              
              
dbLat <- db_data$Lat_lake
dbLat[duplicated(dbLat)] <- NA
uniquedb <- which(!is.na(dbLat))
              
minDistGLOF <- vector()
              
for(i in 1: length(uniquedb)){
                
                cord.dec <- SpatialPoints(cbind(db_data$Lon_lake[uniquedb[i]],db_data$Lat_lake[uniquedb[i]]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
                
                lakeBB <- st_bbox(cord.dec)
                
                r13Dist <- vector()
                r14Dist <- vector()
                r15Dist <- vector()
                
                if(as.numeric(lakeBB[1])>as.numeric(st_bbox(RGI60_13)[1])&&as.numeric(lakeBB[2])>as.numeric(st_bbox(RGI60_13)[2])&&as.numeric(lakeBB[3])<as.numeric(st_bbox(RGI60_13)[3])&&as.numeric(lakeBB[4])<as.numeric(st_bbox(RGI60_13)[4])){
                  
                  cropRGI <- st_crop(RGI60_13,lakeBB+c(-1,-1,1,1))
                  if(dim(cropRGI)[1]>0){  
                    for(k in 1:dim(cropRGI)[1]){
                      
                      r13Dist[k] <- distm(c(db_data$Lon_lake[uniquedb[i]], db_data$Lat_lake[uniquedb[i]]), c(cropRGI$CenLon[k], cropRGI$CenLat[k]), fun = distHaversine)
                      
                    }
                  }
                }
                
                if(as.numeric(lakeBB[1])>as.numeric(st_bbox(RGI60_14)[1])&&as.numeric(lakeBB[2])>as.numeric(st_bbox(RGI60_14)[2])&&as.numeric(lakeBB[3])<as.numeric(st_bbox(RGI60_14)[3])&&as.numeric(lakeBB[4])<as.numeric(st_bbox(RGI60_14)[4])){
                  
                  cropRGI <- st_crop(RGI60_14,lakeBB+c(-1,-1,1,1))
                  
                  if(dim(cropRGI)[1]>0){ 
                    for(k in 1:dim(cropRGI)[1]){
                      
                      r14Dist[k] <- distm(c(db_data$Lon_lake[uniquedb[i]], db_data$Lat_lake[uniquedb[i]]), c(cropRGI$CenLon[k], cropRGI$CenLat[k]), fun = distHaversine)
                      
                    }
                  }
                }
                
                if(as.numeric(lakeBB[1])>as.numeric(st_bbox(RGI60_15)[1])&&as.numeric(lakeBB[2])>as.numeric(st_bbox(RGI60_15)[2])&&as.numeric(lakeBB[3])<as.numeric(st_bbox(RGI60_15)[3])&&as.numeric(lakeBB[4])<as.numeric(st_bbox(RGI60_15)[4])){
                  
                  cropRGI <- st_crop(RGI60_15,lakeBB+c(-1,-1,1,1))
                  
                  if(dim(cropRGI)[1]>0){ 
                    for(k in 1:dim(cropRGI)[1]){
                      
                      r15Dist[k] <- distm(c(db_data$Lon_lake[uniquedb[i]], db_data$Lat_lake[uniquedb[i]]), c(cropRGI$CenLon[k], cropRGI$CenLat[k]), fun = distHaversine)
                      
                    }
                  }
                }
                
                minDistGLOF[i] <- min(min(r13Dist,na.rm=T),min(r14Dist,na.rm=T),min(r15Dist,na.rm=T))
              }

# slope in lake environment
slopVec <- vector()
slopVecGLOF <- vector()

for(k in 1:dim(wangComb)[1]){
slopEnv <- crop(DEM_HMA,st_bbox(wangComb$geometry[k])+c(-0.01,-0.01,0.01,0.01))
slopEnv <- terrain(slopEnv, 'slope', unit = 'degrees', neighbors=8)
slopVec[k] <- cellStats(slopEnv, 'mean')
}  

for(k in 1: length(uniquedb)){
  cord.dec <- SpatialPoints(cbind(db_data$Lon_lake[uniquedb[k]],db_data$Lat_lake[uniquedb[k]]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  lakeBB <- st_bbox(cord.dec)
slopEnv <- crop(DEM_HMA,lakeBB+c(-0.01,-0.01,0.01,0.01))
slopEnv <- terrain(slopEnv, 'slope', unit = 'degrees', neighbors=8)
slopVecGLOF[k] <- cellStats(slopEnv, 'mean')
}

png(file=figures_path&'\\GLOFLakesCharacteristics.png', res = 160,width=1800,height=900)
par(mar=c(5,7,0,1),cex.lab=1.5,cex.axis=1.5)
#par(mfrow=c(3,1))
layout(matrix(c(1,2,3,4), nrow = 1, ncol = 4, byrow = TRUE))
# Area
Areadf <- data.frame(y=c(wangComb$GL_Area/10^6,GLOFArea3),x=c(rep(1,length(wangComb$GL_Area)),rep(2,length(GLOFArea3))))
with(Areadf, boxplot(y~x,log="y",ylim=c(0.001,10),border=c("blue", "red"),xaxt='n',xlab = '',yaxt ='n',ylab=expression(paste("Area [",km^2,"]", sep=""))))
at.y <- outer(1:8, 10^(-3:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1)
grid(nx=NULL, ny=NULL)

myjitter<-jitter(rep(1, length(wangComb$GL_Area/10^6)), amount=0.2)
points(myjitter,wangComb$GL_Area/10^6,  pch=1, col = alpha('blue', 0.006),cex=1.5)
myjitter<-jitter(rep(2, length(GLOFArea3)), amount=0.2)
points(myjitter,GLOFArea3,  pch=1, col= alpha('red', 0.2),cex=1.5)

# Area change of individual lakes - identify lakes in individual inventories that appear in more than one year
Wang19902018match <- which(!is.na(match(lakes1990_sf$GLAKE_ID,lakes2018_sf$GLAKE_ID)))
dA9018_Wang <- (lakes2018_sf$GL_Area[match(lakes1990_sf$GLAKE_ID,lakes2018_sf$GLAKE_ID)[Wang19902018match]] - lakes1990_sf$GL_Area[Wang19902018match]) / 28

# Get area change of GLOF lakes from the Wang DB
WangGLOFmatch <- which(!is.na(match(lakes1990_sf$GLAKE_ID[Wang19902018match],db_data$GL_ID)))
GLOFWangdA <- dA9018_Wang[WangGLOFmatch]

dAreadf <- data.frame(y=c(dA9018_Wang,GLOFWangdA),x=c(rep(1,length(dA9018_Wang)),rep(2,length(GLOFWangdA))))
with(dAreadf, boxplot(y~x,ylim=c(-1500,1500),border=c("blue", "red"),xaxt='n',xlab = '',ylab=c(expression('Area change [m'^"2"*' yr'^"-1"*']'))))
grid(nx=NULL, ny=NULL)
myjitter<-jitter(rep(1, length(dA9018_Wang)), amount=0.2)
points(myjitter,dA9018_Wang,  pch=1, col = alpha('blue', 0.006),cex=1.5)
myjitter<-jitter(rep(2, length(GLOFWangdA)), amount=0.2)
points(myjitter,GLOFWangdA,  pch=1, col= alpha('red', 0.2),cex=1.5)

# Distance to glacier
glDistdf <- data.frame(y=c(lakeglacierDis,minDistGLOF),x=c(rep(1,length(lakeglacierDis)),rep(2,length(minDistGLOF))))
with(glDistdf, boxplot(y~x,ylim=c(0,4500),border=c("blue", "red"),xaxt='n',xlab = '',ylab=c(expression('Distance to glacier [m]'))))
grid(nx=NULL, ny=NULL)
myjitter<-jitter(rep(1, length(lakeglacierDis)), amount=0.2)
points(myjitter,lakeglacierDis,  pch=1, col = alpha('blue', 0.006),cex=1.5)
myjitter<-jitter(rep(2, length(minDistGLOF)), amount=0.2)
points(myjitter,minDistGLOF,  pch=1, col= alpha('red', 0.2),cex=1.5)

# Slopes in glacier environment
glSlodf <- data.frame(y=c(slopVec,slopVecGLOF),x=c(rep(1,length(slopVec)),rep(2,length(slopVecGLOF))))
with(glSlodf, boxplot(y~x,ylim=c(0,45),border=c("blue", "red"),xaxt='n',xlab = '',ylab=c(expression("Periglacial slope ["*~degree*"]"))))
grid(nx=NULL, ny=NULL)
myjitter<-jitter(rep(1, length(slopVec)), amount=0.2)
points(myjitter,slopVec,  pch=1, col = alpha('blue', 0.006),cex=1.5)
myjitter<-jitter(rep(2, length(slopVecGLOF)), amount=0.2)
points(myjitter,slopVecGLOF,  pch=1, col= alpha('red', 0.2),cex=1.5)

dev.off()


# Investigate flow paths of actual GLOFs

#par(mfrow=c(4,4),mar=c(1,2,1,3),mai = c(0.1, 0.1, 0.1, 0.1),cex.lab=1,cex.axis=1)  

flowpathL <- vector() # initialize empty vector for flow length of the GLOF
flowPathH <- vector() # initialize empty vector for drop height
flowPathalpha_mean <- vector()  # initialize empty vector for mean slope
flowPathalpha_max <- vector()   # initialize empty vector for max slope

#missing: 23 + 27 + 88 + 102 + 104

for(k in 1:length(db_data$GF_Id)){
  if(!is.na(db_data$Lon_impact[k]) == T){
idGF <- k

extent_case <- c(min(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),max(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),
            min(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF]),max(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF])) +
  c(-0.03,0.03,-0.03,0.03)
  

DEM_case <- crop(DEM_HMA,extent_case)
LIA_case <- crop(LIA_outline,extent_case)
RGI_case <- st_crop(RGI60_13,extent_case)

#r.fd <- terrain(DEM_case, opt='flowdir')
#r.p <- flowPath(r.fd, c(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]))

#p.xy <- xyFromCell(r.fd,r.p)
#plot(DEM_case)
#points(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF],col='blue',pch=4)
#points(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF],col='blue',pch=4)
#lines(p.xy,col='red', lwd=4)

stpnt <- sp::SpatialPoints(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]))
stpnt_SPDF <- SpatialPointsDataFrame(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]),as.data.frame(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF])))
enpnt <- sp::SpatialPoints(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF]))
enpnt_SPDF <- SpatialPointsDataFrame(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF]),as.data.frame(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF])))

rast_slope <- create_slope_cs(DEM_case, cost_function = 'tobler')

#plot(raster(rast_slope_barrier), col = grey.colors(100))

lcp <- create_lcp(rast_slope, stpnt, enpnt, TRUE)

writeOGR(obj=lcp, dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\FlowPaths\\"&db_data$GF_Id[k]&'_flowpath.shp', layer="lcp", driver="ESRI Shapefile",overwrite=T) 
writeOGR(obj=stpnt_SPDF, dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\FlowPaths\\"&db_data$GF_Id[k]&'_startP.shp', layer="lcp", driver="ESRI Shapefile",overwrite=T) 
writeOGR(obj=enpnt_SPDF, dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\FlowPaths\\"&db_data$GF_Id[k]&'_endP.shp', layer="lcp", driver="ESRI Shapefile",overwrite=T) 

flowpathL[k] <- lengthLine(lcp)
flowPathH[k] <- as.numeric(db_data$Elev_lake[idGF]) - as.numeric(db_data$Elev_impact[idGF])
flowPathalpha_mean[k] <- mean(extract(terrain(DEM_case,'slope','degrees'),lcp)[[1]],na.rm=T)
flowPathalpha_max[k] <- max(extract(terrain(DEM_case,'slope','degrees'),lcp)[[1]],na.rm=T)
}
}

# Example Plots

k <- 626
idGF <- k
extent_case <- c(min(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),max(db_data$Lon_impact[idGF],db_data$Lon_lake[idGF]),
                 min(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF]),max(db_data$Lat_impact[idGF],db_data$Lat_lake[idGF])) +
  c(-0.03,0.03,-0.03,0.03)

DEM_case <- crop(DEM_HMA,extent_case)
LIA_case <- crop(LIA_outline,extent_case)
PF_case <- crop(regionalPF,extent_case)
PF_case[PF_case<0.5] <- NA

ogrInfo("D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Agri.shp")
GLOF_Agri<-readOGR(dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Agri.shp")

ogrInfo("D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Infra.shp")
GLOF_Infra<-readOGR(dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Infra.shp")

ogrInfo("D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Roads.shp")
GLOF_Roads<-readOGR(dsn="D:\\Work\\ICIMODProjects\\GLOF_Database\\GLOF_Roads.shp")


stpnt <- sp::SpatialPoints(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]))
stpnt_SPDF <- SpatialPointsDataFrame(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF]),as.data.frame(cbind(db_data$Lon_lake[idGF],db_data$Lat_lake[idGF])))
enpnt <- sp::SpatialPoints(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF]))
enpnt_SPDF <- SpatialPointsDataFrame(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF]),as.data.frame(cbind(db_data$Lon_impact[idGF],db_data$Lat_impact[idGF])))

rast_slope <- create_slope_cs(DEM_case, cost_function = 'tobler')
lcp <- create_lcp(rast_slope, stpnt, enpnt, TRUE)

tiff(file=figures_path&'\\GLOF_path.tiff', res = 160,width=1200,height=1200)
DEMContour <- rasterToContour(DEM_case)
par(bg=NA,mar=c(2,2,1,1),oma=c(0,0,0,0))
image(DEM_case, axes=T,col = 'white')
plot(DEMContour,add=T,border='grey',lwd=0.5)
plot(wangComb$geometry,add=T,col='blue')
plot(GLOF_Agri,add=T,col='green',lwd=1)
plot(GLOF_Infra,add=T,col='grey',lwd=1)
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}
plot(PF_case,add=T,col=colorRampAlpha(brewer.pal(9,"Blues"), n = 10, alpha=.5), horizontal = TRUE, 
     legend.args = list(text='PZI [-]', side = 1, line = 2),legend.width=1, legend.shrink=0.75,
     smallplot=c(0.1,0.2, 0.6,0.65),legend.key = element_rect(fill = "white"))
plot(RGI60_15$geometry,add=T,border='blue',lwd=1)


plot(GLOF_Roads,add=T,border='grey',lwd=2,lty=3)

if(!is.null(LIA_case)){
  plot(LIA_case, add=T,border='red',lwd=1)}
legend('topleft',legend = c('GLOF path','RGI 6.0','LIA extent','Roads'),col=c('red','blue','red','black'),lty = c(2,1,1,3),lwd =c(2,1,1,2),box.lwd = 0,box.col = "white",bg = "white")
lines(lcp,lwd=3,col='red',lty = 2)
points(stpnt)
points(enpnt)
dev.off()

write.csv(flowpathL,file=figures_path&'\\flowpathL.csv',row.names=F)
write.csv(flowPathH,file=figures_path&'\\flowPathH.csv',row.names=F)
write.csv(flowPathalpha_mean,file=figures_path&'\\flowPathalpha_mean.csv',row.names=F)
write.csv(flowPathalpha_max,file=figures_path&'\\flowPathalpha_max.csv',row.names=F)



# visuzalize dL/dH relation
KaabData <- read.csv(file="D:\\Work\\ICIMODProjects\\GLOF_Database\\dHdL_Kaab.csv")


flowpathL <- c(flowpathL,NA,NA)
flowPathH <- c(flowPathH,NA,NA)

MSize <- as.numeric(db_data$Volume)

MSize <- combAll$Area[match(db_data$GL_ID,combAll$GLIMS_ID)]
MSize[which(is.na(MSize))] <- as.numeric(db_data$Area)[which(is.na(MSize))]                      
             
#MSize[which(is.na(MSize))] <- 5
dfx <- data.frame(dL = flowpathL[which(!is.na(MSize))]/1000,dH = flowPathH[which(!is.na(MSize))]/1000, MSize = MSize[which(!is.na(MSize))])
dfx2 <- data.frame(dL = flowpathL[which(is.na(MSize))]/1000,dH = flowPathH[which(is.na(MSize))]/1000)

dfx_kaab <- data.frame(dL = KaabData$dL/1000,dH = KaabData$dH/1000)



p <- ggplot(dfx, aes(dL, dH)) +
  geom_abline(intercept = 0, slope = 0.09,col='grey') +
  geom_abline(intercept = 0, slope = 0.18,col='grey') +
  geom_abline(intercept = 0, slope = 0.27,col='grey') +
  geom_point(data=dfx2,aes(dL, dH),colour='black',size=2.5) +
  geom_point(aes(colour = log(MSize)),size=2.5) +
  scale_colour_gradient(low = "yellow", high = "red", na.value = NA, name = expression(paste("log(Area) [",m^2,"]", sep=""))) +

  scale_x_continuous(expand = c(0, 0), limits = c(0, 30)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  geom_point(data=dfx_kaab,aes(dL, dH),colour='blue', alpha=0.25) +

  xlab("dL [km]") +
  ylab("dH [km]") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"),
        legend.position = c(0.1, 0.85),
        text = element_text(size = 15),
        legend.key.size = unit(0.5, 'cm'))


p2 <- ggplot(dfx, aes(dL, dH)) +
  geom_abline(intercept = 0, slope = 0.09,col='grey') +
  geom_abline(intercept = 0, slope = 0.18,col='grey') +
  geom_abline(intercept = 0, slope = 0.27,col='grey') +
  geom_point(aes(colour = log(MSize)),size=2.5) +
  scale_colour_gradient(low = "yellow", high = "red", na.value = NA, name = expression(paste("log(Area) [",m^2,"]", sep=""))) +
  geom_point(data=dfx2,aes(dL, dH),colour='black',size=2.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(30, 250)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  geom_point(data=dfx_kaab,aes(dL, dH),colour='blue', alpha=0.25) +
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"),
        legend.position = 'none',
        text = element_text(size = 20))


png(file=figures_path&'\\GLOF_dHdL.png', res = 160,width=1200,height=800)
par(mar=c(2,4.1,2,0.2),cex.lab=1.5,cex.axis=1.5)
print(p)
dev.off()

png(file=figures_path&'\\GLOF_dHdL2.png', res = 160,width=500,height=500)
print(p2)
dev.off()