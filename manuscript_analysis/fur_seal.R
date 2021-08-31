### Analysis of Antarctic Fur Seal foraging areas ### 

pacman::p_load(dplyr, lubridate, sf, track2KBA, ggplot2)

# load data 

# setwd("C:/Users/Martim Bill/Documents/track2iba")

load("C:/Users/Martim Bill/Documents/mIBA_package/all_orig_dev_files/example_data/SG_Seal_FinalData4mIBA_V3.Rdata")

# ses <- subset(Data3, Data3$common_name=="SES")
afs <- subset(Data3, Data3$common_name!="SES") # just fur seals

col1 <- subset(afs, afs$SpecificColony=="BirdIsl") # just Bird Island data
col2 <- subset(afs, afs$SpecificColony=="Husvik")

##########

tracks <- col1 # just Bird Island data
# tracks <- col2

## filter to breeding stage only 

tracks <- subset(tracks, tracks$breed_stage=="Breeding")

n_distinct(tracks$track_id)

if(tracks$SpecificColony[1] == "Husvik"){
  colony <- data.frame(Longitude = -36.7116, Latitude = -54.18)
} else if (tracks$SpecificColony[1] == "BirdIsl"){
  colony <- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
    rename(Longitude=lon_colony,Latitude=lat_colony)
}

colony$Longitude <- as.numeric(colony$Longitude)
colony$Latitude  <- as.numeric(colony$Latitude)

######
# Format for track2KBA
tracks_form <- formatFields(
  tracks, 
  fieldID = "track_id", 
  fieldLat="Latitude", 
  fieldLon="Longitude", 
  fieldDateTime="date", 
  cleanDF = F
  )
str(tracks_form)

######
# filter to years with at least n (10) tracked individuals (rmv 1995)
tracks_form$year <- year(tracks_form$DateTime)
yr_n <- tracks_form %>% 
  group_by(year) %>% 
  summarise(n = n_distinct(ID))
goodyrs <- yr_n %>% filter(n > 10)

tracks_form <- tracks_form %>% filter(year %in% goodyrs$year)

#####
## split into foraging trips

allTrips <- tripSplit(
  tracks_form, 
  colony=colony, 
  innerBuff = 5, 
  returnBuff = 50, 
  duration=12, 
  rmNonTrip = T
  )

## check yearly sample sizes 
allTrips@data %>% mutate(year = lubridate::year(DateTime)) %>% group_by(year) %>% summarise(N = n_distinct(ID))

# Tripmap <- mapTrips(allTrips, colony)
Tripmap <- mapTrips(allTrips[allTrips$ID %in% unique(allTrips$ID)[35:46], ], colony)

# ggsave( "C:/Users/Martim Bill/Documents/mIBA_package/figures/fur_seals/trips.png", width = 8, height=6)

######
# source("tripSummary.R")

# Trips <- subset(allTrips, allTrips$Returns == "Yes" ) # removes some 9 seals
Trips <- allTrips
TripSum <- tripSummary(Trips, colony)
TripSum

frange <- median(TripSum$max_dist)
frange
c(min(TripSum$max_dist), max(TripSum$max_dist))
fdur <- median(na.omit(TripSum$duration))/24
fdur
c(min(na.omit(TripSum$duration)), max(na.omit(TripSum$duration)))/24

###### Identify and filter out duplicate DT stamps #####

Trips_list <- vector(mode="list", dplyr::n_distinct(Trips$ID))
for(i in 1:dplyr::n_distinct(Trips$ID)){
  one <- Trips[Trips$ID==unique(Trips$ID)[i], ]
  one$duplicated <- duplicated(one$DateTime)
  Trips_list[[i]] <- one
}
Trips <- do.call(rbind, Trips_list)

Trips <- Trips[Trips$duplicated == FALSE, ]

## 3. ####

Trips <- projectTracks(dataGroup = Trips, projType = "azim", custom=T, reproject = T)

### findScale (get average foraging range, a list of H-value options, and test whether desired grid cell for kernel estimation makes sense given movement scale/tracking resolution) ~~~~~~~~~~~~~~~

HVALS <- findScale(Trips,
  scaleARS = T,
  sumTrips = TripSum
)
HVALS

## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~
h <- HVALS$mag
h

Trips <- Trips[Trips$ColDist > 5, ] # remove trip start and end points near colony

KDE.Surface <- estSpaceUse(tracks=Trips, scale = h, levelUD = 50, polyOut=T)
# KDE.Surface <- estSpaceUse(tracks=Trips, scale = h, levelUD = 50, polyOut=F)
# saveRDS(KDE.Surface, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\KDE.rds")
KDE.Surface <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\KDE.rds")

KDEmap <- mapKDE(KDE.Surface$UDPolygons, colony = colony)

n <- length(KDE.Surface$KDE.Surface)
# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/fur_seals/indcores_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6)



## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
repr <- repAssess(Trips, KDE=KDE.Surface$KDE.Surface, levelUD=50, iteration=200, nCores=2, bootTable = T)
Sys.time() - before
#
# rep_value <- repr[[1]]$out
rep_value <- 95.91
# rep_value <- 96.1 # no 1995 data

## 7. ####
### findSite (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

popsize <- 64545 # number of adult female fur seals

# sites <- findSite(KDE.Surface, Represent=repr$out, polyOut = F)
sites <- findSite(KDE.Surface$KDE.Surface, represent=rep_value, levelUD=50, popSize = popsize, polyOut = T) 
sites
# saveRDS(sites, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\site.rds")

sites <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\site.rds")

sitemap <- mapSite(Site=sites, colony = colony, show=F)

# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/fur_seals/potSite_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6 )

# how many seals use the 'site'?
sites %>% dplyr::filter(.data$potentialSite==TRUE) %>% 
  summarise(
    max_animals = max(na.omit(N_animals)),
    min_animals = min(na.omit(N_animals))
  )

# % of global pop 
sites$max_animals/700000*100  # min
sites$max_animals/1000000*100 # max

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add in some background maps for context 

library(ggmap)

# aggs <- sites[sites$N_animals > 0, ]
fs_aggs <- sites[sites$N_animals >= (0.1 * popsize), ]

xmin <- st_bbox(fs_aggs)[[1]] - 0.1
xmax <- st_bbox(fs_aggs)[[3]] + 0.15
ymin <- st_bbox(fs_aggs)[[2]] - 0.1
ymax <- st_bbox(fs_aggs)[[4]] + 0.1

gmap <- ggmap::get_map(location=c(xmin, ymin, xmax, ymax), zoom=10, maptype = "satellite")

fs_aggs <- st_transform(fs_aggs, crs = 3857)

# plot(st_transform(colony_sf, crs = 3857), add=T, col=2) # doesn't work for some reason!
colony_sf <- st_transform(st_as_sf(colony, coords = c("Longitude", "Latitude"), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'), 3857)
# 
# Use the function:
source("C:/Users/Martim Bill/Documents/R/source_scripts/ggmap_bbox.r")

gmap <- ggmap_bbox(gmap)

fs_map <- ggmap(gmap) + 
  coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  geom_sf(data = fs_aggs, inherit.aes = FALSE, aes(fill=N_animals/1000), color=NA) + 
  geom_sf(data=colony_sf, inherit.aes = FALSE,  color="black", fill="red", shape=23, size=3) +
  scale_fill_gradientn(colours=sf.colors(n=3)) +
  guides( 
    fill = guide_colorbar(
      title.position="top",
      title = "N animals (thousands)",
      barwidth  = 10,
      barheight = 1.75)) +
  theme(
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=11, color="black"),
    axis.title=element_text(size=16),
    legend.direction = "horizontal",
    legend.position=c(0.26, 0.18),
    legend.title=element_text(size=15),
    legend.text = element_text(size = 13)) +
  ylab("Latitude") +
  xlab("Longitude")
fs_map

# saveRDS(fs_map, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\FS_site_ggmap.rds")


## SAVE ##
ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/fur_seals/Xpotsite_N_animals", "h", round(h), "_", "n",n, ".png"), fs_map, width = 10, height=8 )


## get area of potentialSite ## 

site <- sites %>% filter( potentialKBA == T ) %>% group_by(potentialKBA) %>% summarise(N_animals = max(N_animals))
plot(site)
st_area(site) / 1000^2 # sq km
