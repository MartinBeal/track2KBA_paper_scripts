####### Example script of methods used in Lascelles et al. 2016 to identify IBA areas ########
## Requires the source script of user-defined functions ('source_scripts(update_06.18).r) ##
## Downloaded: 25-09-2018 by Martin Beal 


############ data preparation ##############

source("Lascelles_16/source_scripts(update_06.18).r")

kpacks <- c('sp','maptools', 'rgdal', 'adehabitatHR','geosphere','fields','spatstat','maps','rgeos','mapdata','raster', 'igraph') #MB# added raster package as this is required by 'polyCount' fxn and 'igraph' package for 'thresholdRaster' fxn

new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]

if(length(new.packs)) install.packages(new.packs)
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)
DataGroup <- read.csv("Lascelles_16/lascelles_test.data.csv")
 # replace "table.csv" by the name of your table with tracking data (see above)

# to see the data
head(DataGroup)
plot(Latitude~Longitude, data=DataGroup, asp=1)
map("world", add=T)

# calculate TrackTime field
DataGroup$DateTime <- paste(DataGroup$DateGMT, DataGroup$TimeGMT, sep= " ")
DataGroup$DateTime <- as.POSIXct(strptime(DataGroup$DateTime, "%d/%m/%Y %H:%M:%S "), "GMT") ##see other options for date format in strptime help
DataGroup$TrackTime <- as.numeric(DataGroup$DateTime)
head(DataGroup)

# create colony data.frame
Colony <- data.frame(Longitude = -38.06, Latitude = -54) # replace by Longitude and Latitude values of the colony
#or
#Colony=read.csv("nests.csv") # a data.frame with a Latitude and Longitude value per individual ID nest
# see individual tracks and colony location
Tracks <- unique(DataGroup$ID)
DataGroup$TrackId <- as.factor(DataGroup$ID)
plot(Latitude~Longitude, data=DataGroup, col="white", asp=1)
for(i in Tracks)
{
  TempTrack <- subset(DataGroup, DataGroup$TrackId == i)
  lines(Latitude~Longitude, data=TempTrack, col=TrackId)
}
map("world", add=T)

with(Colony, points(Longitude, Latitude, pch=19, cex=3))


#### split the tracks into trips and removes locations on the colony ####
DataGroup.Df <- DataGroup
DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
DgProj <- CRS(paste("+proj=laea +lon_0=", Colony$Longitude, " +lat_0=", Colony$Latitude, sep=""))
DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
DataGroupTrips <- NULL
for(i in Tracks)
{
  TempTrack <- subset(DataGroup, DataGroup$TrackId == i)
  TempTrips <- tripSplit(Track = TempTrack, Colony = Colony, InnerBuff = 15, ReturnBuff = 45, plotit=T ,MidPoint = F, nests=F)
  if(which(Tracks == i) == 1)
  {
    DataGroupTrips <- TempTrips
  } else {DataGroupTrips <- spRbind(DataGroupTrips, TempTrips)}
}

names(DataGroupTrips@data)[names(DataGroupTrips@data) %in% c("DataGroup.Longitude" , "DataGroup.Latitude")] <- c("X", "Y")

DataGroupTrips <- DataGroupTrips[DataGroupTrips$Returns != "N",]
DataGroupTrips <- DataGroupTrips[DataGroupTrips$trip_id != "-1",]
# MB # remove trips with less than 5 points
DataGroupTrips <- DataGroupTrips[!DataGroupTrips$trip_id %in% names(which(table(DataGroupTrips$trip_id) < 5)),]
head(DataGroupTrips@data)
DataGroupTrips$originalID <- DataGroupTrips$ID #save the original IDs with a new name
DataGroupTrips$ID <- DataGroupTrips$trip_id #reset the ID field to individual trips rather than individual birds

# summary of the trip data
tsum <- tripSummary(DataGroupTrips, Colony=Colony, nests=F)
tsum[order(tsum$departure),]


######### identify the appropriate ARS Scale #########

Scales <- seq(10,250,5)
fpt.scales <- scaleARS(DataGroup=DataGroupTrips,Scales=Scales,Peak="Flexible")

####### calculate the kernels and export the resultant polygons as shapefiles ########

datagroupsUDd <- batchUD(DataGroupTrips, Scale=fpt.scales/2, UDLev=50)

td <- "output" # path to the folder where the shapefile will be saved
writeOGR(datagroupsUDd, td, "testUD", driver="ESRI Shapefile") #MB# switched 'td' and '"testUD"' to send file to folder
datagroupsUDd <- readOGR(td, "testUD", encoding="ESRI Shapefile")


if ("originalID" %in% names(DataGroupTrips@data))
{
  TableKernel <- datagroupsUDd@data
  linkids <- with(DataGroupTrips@data,aggregate(ID, list(ID=ID,originalID=originalID), length))
  TableKernel <- merge(TableKernel,linkids[,c("ID","originalID")])[,c("Name_0","Name_1","ID","originalID")] [order(TableKernel$Name_1),]
  TableKernel
  datagroupsUDd@data <- TableKernel
}
datagroupsUDd@data

###### variance test #######
bird_string <- as.character(datagroupsUDd$originalID)
vt <- varianceTest(datagroupsUDd, bird_string, Iteration=10)
vt
### to choose randomly just one trip per individual
if (vt < 0.25)
{
  bird_idtrip=datagroupsUDd@data
  birds=unique(bird_idtrip$originalID)
  trips=numeric()
  set.seed(1)
  for (x in 1:length(birds)) trips=c(trips, as.character(sample(bird_idtrip[bird_idtrip$originalID==(birds[x]),]$ID,1)))
  DataGroupTrips2=DataGroupTrips[DataGroupTrips$ID%in%trips,]
  datagroupsUDd2=batchUD(DataGroupTrips2, Scale=fpt.scales/2, UDLev=50)
  datagroupsUDd=datagroupsUDd2
  DataGroupTrips=DataGroupTrips2
}
# Bootstrap
# boot_out<-bootstrap(DataGroupTrips, Scale=fpt.scales/2, Iteration=50)
# boot_out

# polygons count
count_UD <- polyCount(datagroupsUDd, Res = 0.05)
writeRaster(count_UD, "output/count_UD_ras.tif")

# Threshold_raster
mIBA_site <- thresholdRaster(count_UD, Threshold = 15)
writeOGR(mIBA_site, layer="spcies_mIBA_candidatesite_thresh.10p", dsn="output", driver="ESRI Shapefile", verbose=TRUE)


