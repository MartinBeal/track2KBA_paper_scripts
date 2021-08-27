## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Analysis of green turtle movements post-laying ## 

pacman::p_load(track2KBA, ctmm, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter)

## Data input ~~~~~~~~~~~~~~~~~~
tracks_f <- readRDS("C:/Users/Martim Bill/Documents/mIBA_package/data/green_turtles/foraging_periods/foraging_unfiltered_GPS_PTT_n23.rds")

bbox <- raster::shapefile("C:/Users/Martim Bill/Documents/mIBA_package/data/green_turtles/west_africa_bbox.shp")
bbox@proj4string <- CRS(SRS_string = "EPSG:4326")

## View pre-filtered data 

# fSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_f$Longitude, tracks_f$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data=tracks_f) # foraging

# mapview::mapview(fSP)


#####

tracks_f <- tracks_f %>% rename(Argos.Location.Quality = Location.Quality, GPS.HDOP=HDOP) %>% mutate(
  sensor = ifelse(Argos.Location.Quality == "", "GPS", "Argos doppler shift")
) 
 
## Remove PTT points of lowest quality  ~~~~~~~~~~~~~~~~~~~~~
tracks_f <- tracks_f[!tracks_f$Argos.Location.Quality %in% c("Z","B","A"),]
# tracks <- tracks[!tracks$Argos.Location.Quality %in% c("Z","B","A","0"),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Speed and inner angle filter for GPS data with SDLfilter ~~~~~~~~~~~~~~~~~~~~~

gps <- tracks_f[tracks_f$sensor == "GPS", ]

colnames(gps)[c(1,5,6,10)] <- c("id","lat","lon","qi")

## SDLfilter steps ## 
# gps_nodups <- dupfilter(gps) # removes 2.3% of points (87)
gps_nodups <- dupfilter(gps, step.time=.1/60, step.dist=0.0001) # removes 2.3% of points (87)

vmax <- vmax(gps_nodups, qi=5, prob=0.98)    # maximum linear speed   # 9.0 km/h data-driven speed filtere
vmaxlp  <- vmaxlp(gps_nodups, qi=5, prob=0.98) # maximum 'loop-speed'

gps_fltrd <- ddfilter(gps_nodups, vmax=vmax, vmaxlp = vmaxlp , qi=4, method=1)
# gps_fltrd <- ddfilter.speed(gps_nodups, vmax=vmax, method=1)  # only speed filter
# gps_fltrd <- ddfilter.loop(gps_nodups, maxvlp = maxvlp, qi=4)

# gps_fltrdSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(gps_fltrd$lon, gps_fltrd$lat), proj4string=CRS("+proj=longlat + datum=wgs84")), data=gps_fltrd) # foraging
# gpsSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(gps$lon, gps$lat), proj4string=CRS("+proj=longlat + datum=wgs84")), data=gps) # foraging
# mapview::mapview(gps_fltrdSP) # compare filtered result to raw 
# mapview::mapview(gpsSP)

# table(gps$id) - table(gps_fltrd$id) # number of points rmvd per individual

## View filtered points
# anti_gps <- anti_join(gps, gps_fltrd)
# anti_gpsSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(anti_gps$lon, anti_gps$lat), proj4string=CRS("+proj=longlat + datum=wgs84")), data=anti_gps) # foraging
# mapview::mapview(anti_gpsSP)

gps_fltrd <- formatFields(gps_fltrd, fieldID = "id", fieldLat="lat", fieldLon="lon", fieldDateTime="DateTime") #%>% 
  #dplyr::select(-c(qi, pTime, sTime, pDist, sDist, pSpeed, sSpeed, inAng))

## Recombine filtered GPS data with PTT data
tracks_f <- tracks_f %>% filter(sensor == "Argos doppler shift") %>% full_join(gps_fltrd, by=) %>% arrange(ID, DateTime)

fSP2 <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_f$Longitude, tracks_f$Latitude), proj4string=CRS(SRS_string = "EPSG:4326")), data=tracks_f) # foraging


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## track2KBA analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TD <- projectTracks(fSP2, projType = "azim", custom=T)

HVALS <- findScale(TD)
HVALS

h <- HVALS$href
h

KDE <- estSpaceUse(TD, scale=h, polyOut=T, res = 1.5)
# saveRDS(KDE, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\KDE.rds")

n <- length(KDE$KDE.Surface)
ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/white_storks/indcores_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6)

# KDE <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\KDE.rds")

KDEmap <- mapKDE(KDE$UDPolygons)

n <- length(KDE$KDE)
  
# dev.new()

## Assess representativeness ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# before <- Sys.time()
# represent <- repAssess(TD, KDE$KDE.Surface, iteration=800, nCores=2, avgMethod = "mean", bootTable = T)
# Sys.time() - before
# rep_value <- represent[[1]]$out
rep_value <- 32.4

popsize <- 18573

aggs <- findSite(KDE=KDE$KDE.Surface, popSize = popsize, represent = rep_value, polyOut = T)
# aggs <- findSite(KDE=KDE$KDE.Surface, represent = rep_value)
# saveRDS(aggs, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\Site.rds")

# aggs <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\Site.rds")
mapSite(aggs)

## Make Maps ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722)
poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")


# custom map of individual ranges # 
coordsets <- st_bbox(KDE$UDPolygons)
coordsets[1] <- coordsets[1] - .9
# coordsets[3] <- coordsets[1] + 1

UDPLOT <- ggplot(KDE$UDPolygons) +
  geom_sf(data=KDE$UDPolygons, aes(col=id), size=1, fill=NA) +
  coord_sf(xlim = c(coordsets$xmin - .9, coordsets$xmax + .1), ylim = c(coordsets$ymin, coordsets$ymax), expand = TRUE) +
  borders("world",fill=scales::alpha("dark grey", 0.6), colour="grey20") +
  theme(
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=12, color="black"),
    axis.title=element_text(size=16),
    legend.position = "none") +
  ylab("Latitude") +
  xlab("Longitude")
UDPLOT

# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/sea_turtles/forage_indcores_", "h", round(h), "_", "n",n, ".png"), width = 6, height=12)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Map of overlap areas with base map # 
library(ggmap)

xmin <- min(TD$Longitude)
xmax <- max(TD$Longitude)
ymin <- min(TD$Latitude) 
ymax <- max(TD$Latitude) 

gt_gmap <- ggmap::get_map(location=c(left=coordsets$xmin[[1]]-.8, bottom=coordsets$ymin[[1]]-.3, right=coordsets$xmax[[1]]+.5, top=coordsets$ymax[[1]]+.4), zoom=8, maptype = "satellite")

# ggmap(gmap)

# Use the function:
source("C:/Users/Martim Bill/Documents/R/source_scripts/ggmap_bbox.r")

gt_aggs <- st_transform(aggs[aggs$N_animals > 0, ], crs = 3857)
# MPAs_3857 <- st_transform(MPAstack, crs = 3857)  #MPAstack from turtle_map_overlap_zones
poilao <- st_transform(poilao, crs = 3857)

gt_gmap <- ggmap_bbox(gt_gmap)

gt_map <- ggmap(gt_gmap) + 
  coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  # geom_sf(data = MPAs_3857, inherit.aes = FALSE, color="black", fill=NA) +
  geom_sf(data = gt_aggs, inherit.aes = FALSE, aes(fill=N_IND), color=NA) +     # number of tracked INDs
  # geom_sf(data = gt_aggs, inherit.aes = FALSE, aes(fill=N_animals), color=NA) + # number of animals
  geom_sf(data = poilao, inherit.aes = FALSE, color="black", fill="red", size = 3.5, shape=23) +
  scale_fill_gradientn(colours=sf.colors(n=3), breaks=c(1,2,3)) +
  # scale_fill_gradientn(colours=sf.colors(n=3), breaks=c(1,2,3)) +
  guides( 
    fill = guide_colorbar(
      barwidth  = 1.75,
      title = "N IND",
      # title = "N animals",
      barheight = 8)) +
  theme(
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=11, color="black"),
    axis.title=element_text(size=16),
    legend.position=c(.85, .5),
    legend.title=element_text(size=15),
    legend.text = element_text(size = 13)) +
  ylab("Latitude") +
  xlab("Longitude")
gt_map

# saveRDS(gt_map, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\GT_Site_ggmap.rds")
saveRDS(gt_map, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\GT_Site_ggmap_NINDs.rds")


## SAVE ##
ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/sea_turtles/forage_overlap_", "h", round(h), "_", "n",n, ".png"), width = 6, height=12 )
