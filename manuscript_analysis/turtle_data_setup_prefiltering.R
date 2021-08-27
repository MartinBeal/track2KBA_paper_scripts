## Turtle pre-pre-processing script --- selecting right species, individuals of interest, and separating foraging from migratory periods based on individuals' schedules ## 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pacman::p_load(track2KBA, dplyr, sp, stringr, SDLfilter)

## Data input ~~~~~~~~~~~~~~~~~~
raw <- read.csv("data/tracks/green_turtles_all_nofilter_2018-19.csv")
tracksum <- read.csv("C:/Users/Martim Bill/Documents/mIBA_package/data/green_turtles/summary_csv.csv", stringsAsFactors = F) # deployment summary 


## Identifying stationary vs. migratory periods ##

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tracks <- formatFields(raw, field_ID = "Tag_ID", field_Lat="Latitude", field_Lon="Longitude", field_Date="UTC_Date", field_Time="UTC_Time")

tracks$ID <- do.call(rbind,  stringr::str_split(tracks$ID, pattern = stringr::fixed(":")))[,2] # remove colon and number before it 

# remove points from NZ (calibration)
tracksSP <- sp::SpatialPointsDataFrame(sp::SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=sp::CRS("+proj=longlat + datum=wgs84")), data=tracks)
tracksSP <- tracksSP[bbox,] ## Keep only points in bbox (a polygon around West Africa)

# convert to move object
tracks <- data.frame(tracksSP)
tracks <- tracks[order(tracks$ID, tracks$DateTime), ] ## order date time stamps within each individual

# remove Leatherbacks and Hawksbills
tracks <- tracks[!tracks$ID %in% c("197219","197220","182460","60861", "60862","60895 ","197136","197139","197140","197141 "), ]
# remove data sets without a post-migration period
tracks <- tracks[!tracks$ID %in% c("182455","182456","60900","60866","60888", "182461", "197137", "205284", "205285"), ]
# remove individual(s) with fewer than 20 points
tracks <- tracks[!tracks$ID %in% c("60868"), ]


# Remove duplicate time stamps #
tracks_list <- vector(mode="list", dplyr::n_distinct(tracks$ID))
for(i in 1:dplyr::n_distinct(tracks$ID)){
  one <- tracks[tracks$ID==unique(tracks$ID)[i], ]
  one$duplicated <- duplicated(one$DateTime)
  tracks_list[[i]] <- one
}
tracks <- do.call(rbind, tracks_list)

tracks <- tracks[tracks$duplicated == FALSE, ] 

# optionally use only PTT (removes > 1/2 data for 8 indivs.)
# tracks <- dplyr::filter(tracks, sensor == "Argos Doppler shift") %>% dplyr::select(ID, DateTime, Argos.Location.Quality, Longitude, Latitude, sensor)


## Split data into migratory, and post migratory periods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ID_list <- list()
forage_list  <- list()
migrate_list <- list()

for(i in 1:n_distinct(tracks$ID)){
  
  one <- tracks[tracks$ID == unique(tracks$ID)[i], ]
  ID <- unique(tracks$ID)[i]
  
  start.mig <- as.Date(tracksum[tracksum$PTT==ID,]$date.leaving.Poilão)
  end.mig   <- as.Date(tracksum[tracksum$PTT==ID,]$first.date.at.FG)
  migration <- one[ (one$DateTime > start.mig & one$DateTime < end.mig), ]
  foraging  <- one[ (one$DateTime > end.mig), ]
  
  ID_list[[i]]      <- ID 
  forage_list[[i]]  <- foraging
  migrate_list[[i]] <- migration
}

# summarise number of points of post-migr. and migration for each individ.
summary <- data.frame(
  ID=do.call(rbind, lapply(ID_list, function(x) x[[1]])), 
  n_migrate=do.call(rbind, lapply(migrate_list, function(x) nrow(x))),
  n_forage=do.call(rbind, lapply(forage_list, function(x) nrow(x)))
)

tracks_m <- do.call(rbind, migrate_list)
tracks_f <- do.call(rbind, forage_list )


fSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_f$Longitude, tracks_f$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data=tracks_f) # foraging
mSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_m$Longitude, tracks_m$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data=tracks_m) # migration

dev.new()
mapview::mapview(fSP)
# mapview::mapview(mSP)


## SAVE ## 
saveRDS(tracks_f, "C:/Users/Martim Bill/Documents/mIBA_package/data/green_turtles/foraging_periods/foraging_unfiltered_GPS_PTT_n23.rds")
