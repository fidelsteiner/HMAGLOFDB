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

################ 
# paths and raw data
################

data_path <- 'D:\\Work\\Repositories\\HMAGLOFDB\\Database\\GLOFs'       # location of the database on local station
db_file <- 'HMAGLOFDB_v0.1.csv'                                         # file name of database

DEM_path <- 'D:\\Work\\GeospatialData\\HMA\\SRTM\\VoidFilled'           # location of DEM for topographic analysis
DEM_file <- 'SRTM_Corrected_Extended_HMA.tif'                           # DEM filename

HKH_path <- 'D:\\Work\\GeospatialData\\HMA\\HKH_Boundary'               # location of HKH Outline (http://rds.icimod.org/Home/DataDetail?metadataId=3924&searchlist=True)
HKH_file <- 'HKH_Outline.shp'
HKH_outline<-readOGR(dsn=HKH_path&'\\'&HKH_file)                       # read in HKH outline
HKH_outline <- spTransform(HKH_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")                          


lakesHKH_path <- 'D:\\Work\\Repositories\\HMAGLOFDB\\Database\\Lakes\\HKHLAKEDB'   # location of the HKH lake database
lakes_file_2000 <- 'HKH_2000.shp'
lakes_file_2005 <- 'HKH_2005.shp'
lakes_file_2010 <- 'HKH_2010.shp'
lakes_file_2020 <- 'HKH_2020.shp'

lakesHiMAG_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacialLakes\\Chen2021\\Hi-MAG database'   # location of HiMAG (DOI: 10.5194/essd-13-741-2021)
lakesHiMAG_file_2007 <- 'HiMAG_database_2007.shp'
lakesHiMAG_file_2008 <- 'HiMAG_database_2008.shp'
lakesHiMAG_file_2009 <- 'HiMAG_database_2009.shp'
lakesHiMAG_file_2010 <- 'HiMAG_database_2010.shp'
lakesHiMAG_file_2011 <- 'HiMAG_database_2011.shp'
lakesHiMAG_file_2012 <- 'HiMAG_database_2012.shp'
lakesHiMAG_file_2013 <- 'HiMAG_database_2013.shp'
lakesHiMAG_file_2014 <- 'HiMAG_database_2014.shp'
lakesHiMAG_file_2015 <- 'HiMAG_database_2015.shp'
lakesHiMAG_file_2016 <- 'HiMAG_database_2016.shp'
lakesHiMAG_file_2017 <- 'HiMAG_database_2017.shp'

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


################ 
# Evaluate Lake databases
################

sf::sf_use_s2(FALSE)

ogrInfo(lakesHKH_path&'\\'&lakes_file_2000)
lakes2000_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2000)
lakes2000_outline <- spTransform(lakes2000_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

ogrInfo(lakesHKH_path&'\\'&lakes_file_2005)
lakes2005_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2005)
lakes2005_outline <- spTransform(lakes2005_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2005_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2005)

ogrInfo(lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_outline <- readOGR(dsn=lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_outline <- spTransform(lakes2020_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")
lakes2020_sf <- read_sf(lakesHKH_path&'\\'&lakes_file_2020)
lakes2020_sf <- st_transform(lakes2020_sf,crs=st_crs(lakes2005_sf))

meow_eez_intersection <- st_intersection(lakes2005_sf,lakes2020_sf)

#lakes2000_outline$Altitude <-raster::extract(DEM_HMA, SpatialPoints(cbind(lakes2000_outline$Longitude,lakes2000_outline$Latitude)), sp = T)[]

#hist(lakes2000_outline$Altitude$SRTM_Corrected_Extended_HMA,breaks = c(2000,2500,3000,3500,4000,4500,5000,5500,6000,6500))


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
lakeFormationRecorded <- length(which(db_data$Driver_lake!='unknown')) / length(db_data$Year_exact) * 100 
GLOFFormationRecorded <- length(which(db_data$Driver_GLOF!='unknown')) / length(db_data$Year_exact) * 100 
GLOFMechanismRecorded <- length(which(db_data$Mechanism!='unknown')) / length(db_data$Year_exact) * 100

# Data on discharge/volume
volRecorded <- length(which(!is.na(as.numeric(db_data$Volume)))) / length(db_data$Year_exact) * 100 
QRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_water)))) / length(db_data$Year_exact) * 100 
QsRecorded <- length(which(!is.na(as.numeric(db_data$Discharge_solid)))) / length(db_data$Year_exact) * 100 

################ 
# Extract topographic data
################

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