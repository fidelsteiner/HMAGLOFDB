################################################################################
# Analysis of all GLOF events from the HMAGLOFDB
# 
# HMAGLOFDB_Analysis.R
#
# ReadMe: 
# Read and analyze GLOF events
# 1) Read out database
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

################ 
# paths and raw data
################

data_path <- 'D:\\Work\\ICIMODProjects\\GLOF_Database'
db_file <- 'HMAGLOFDB_20220205.csv'

DEM_path <- 'D:\\Work\\GeospatialData\\HMA\\SRTM\\VoidFilled'
DEM_file <- 'SRTM_Corrected_Extended_HMA.tif'

outline_path <- 'D:\\Work\\GeospatialData\\HMA\\GlacierChangeData\\Lee2021\\SuppInfo_vectors\\SuppInfo_vectors'
outline_file <- 'glaciers_LIA.shp'

lakes_path <- 'D:\\Work\\ICIMODProjects\\TBWG\\LakeStandardization\\LakeDatabase'
lakes_file <- 'Glacial_lake_EH_2020.shp'

cls <- c(Elev_lake="numeric", Elev_impact="numeric")
db_data <-read.csv(data_path&'\\'&db_file,header = T,sep=',')

DEM_HMA <- raster(DEM_path&'\\'&DEM_file)

ogrInfo(lakes_path&'\\'&lakes_file)
lakes_outline<-readOGR(dsn=lakes_path&'\\'&lakes_file)
lakes_outline <- spTransform(lakes_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

ogrInfo(outline_path&'\\'&outline_file)
LIA_outline<-readOGR(dsn=outline_path&'\\'&outline_file)
LIA_outline <- spTransform(LIA_outline,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

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