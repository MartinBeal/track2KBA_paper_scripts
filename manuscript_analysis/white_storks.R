## White Storks - Portugal ##

pacman::p_load(dplyr, lubridate, sf, track2KBA, ggplot2, ggmap, raster)


# allTD <- read.csv("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\White Stork 2016_2019_9s_20200928\\White Stork 2016_2019_9s_20200928.csv")

# migdates <- read.csv("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\all_migration_dates.csv")

# filter to only migrant individuals
# migTD <- allTD[allTD$individual.local.identifier %in% migdates$ID, ] %>% 
#   dplyr::select(
#   individual.local.identifier, Age, burst, timestamp, location.long, location.lat, Mean_ground.speed, Record_type, Year
# ) %>% rename(ID = individual.local.identifier)

# saveRDS(migTD, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\migratory_storks.rds")


## Data from migrant birds only ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

migTD <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\migratory_storks.rds")

migdates <- read.csv("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\all_migration_dates.csv")

# some funny dates in second spring so need to isolate to convert to Date
migdates <- migdates %>%  mutate(
  start_autumn_1 = as.Date(start_autumn_1),
  end_autumn_1   = as.Date(end_autumn_1),
  start_spring_1 = as.Date(start_spring_1),
  end_spring_1   = as.Date(end_spring_1),
  start_autumn_2 = as.Date(start_autumn_2),
  end_autumn_2   = as.Date(end_autumn_2),
  start_spring_2 = as.Date(ifelse(start_spring_2 == "", NA, start_spring_2)),
  end_spring_2   = as.Date(ifelse(end_spring_2 == "", NA, end_spring_2))
  )

migTD <- migTD %>% left_join(migdates, by=c('ID'))

# one <- migTD %>% filter(ID == unique(ID)[1]) %>% st_as_sf(coords = c("location.long","location.lat"), crs = 4326, agr = "constant")
# mapview::mapview(one)

# filter to only birds with at least the first full autumn migration

tracks <- migTD %>% filter(!is.na(end_autumn_1))

# filter to only migratory periods
tracks <- tracks %>% filter( 
  (timestamp > start_autumn_1) & (timestamp < end_autumn_1)   # 1st autumn migration
  | (timestamp > start_spring_1) & (timestamp < end_spring_1) # 1st spring migration
  | (timestamp > start_autumn_2) & (timestamp < end_autumn_2) # 2nd autumn mig.
  | (timestamp > start_spring_2) & (timestamp < end_spring_2) # 2nd spring mig.
  )

# For juvenile birds, retain only first autumn and spring migrations
tracks <- tracks %>% mutate(
  filter = ifelse((age == 'juvenile') & (timestamp > start_autumn_2), TRUE, FALSE),
  filter = ifelse(is.na(filter), FALSE, filter)
) %>% filter(
  filter == FALSE
) %>% dplyr::select(-filter)

# summarise number of IDs with 1, 2, 3 migrations for autumn and spring
# tracks %>% filter(!is.na(end_autumn_1)) %>% summarise(n_distinct(ID))
# tracks %>% filter(!is.na(end_spring_1)) %>% summarise(n_distinct(ID))
# tracks %>% filter(!is.na(end_autumn_2) & age != "juvenile") %>% summarise(n_distinct(ID))
# tracks %>% filter(!is.na(end_spring_2) & age != "juvenile") %>% summarise(n_distinct(ID))


## rough downsample to speed things up ##
# tracks <- tracks %>%
#   filter(row_number() %% 2 == 1) # retain ever nth row

## track2KBA analysis ## 

rm(migTD) # save memory

# analyze all migrations, only fall, or only spring?
analysis <- "all"
# analysis <- "fall"
# analysis <- "spring"
# analysis <- "juvenile"
# analysis <- "adult"

if(analysis == "fall"){
  tracks <- tracks %>% filter( 
    ((timestamp > start_autumn_1) | (timestamp < end_autumn_1)) | 
      ((timestamp > start_autumn_2) | (timestamp < end_autumn_2))
    )
} else if(analysis == "spring"){
  tracks <- tracks %>% filter( 
    ((timestamp > start_spring_1) | (timestamp < end_spring_1)) | 
      ((timestamp > start_spring_2) | (timestamp < end_spring_2))
  )
} else if(analysis == "juvenile"){
  tracks <- tracks %>% filter( 
    age == "juvenile"
  )
} else if(analysis == "adult"){
  tracks <- tracks %>% filter( 
    age == "adult" | age == "immature"
  )
}

## Calculate average time lag btwn points ## 
tracks_move <- move::move(
  x      = tracks$Longitude, 
  y      = tracks$Latitude, 
  time   = tracks$DateTime, 
  proj   = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),  # common project. for all individuals, custom set by tripSplit
  animal = tracks$ID
)

lags <- lapply(seq_along(split(tracks_move)), function(x) {
  one <- tracks_move[[x]]
  ID  <- unique(tracks_move@trackId)[x]
  lag <- timeLag(one, units="mins")
  lag <- data.frame(
    ID = ID, mean=mean(lag), sd=sd(lag), min(lag), max(lag), median(lag),
    q1=quantile(lag, 1/4), q3 = quantile(lag, 3/4), IQR(lag))
  return(lag)
})

lags <- do.call(rbind, lags)
lags

mean(lags$median.lag.)
sd(lags$median.lag.)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Track2KBA analysis ## 

# Format data 
tracks <- formatFields(tracks, fieldID = "ID", fieldDateTime = "timestamp", fieldLat = "location.lat", fieldLon = "location.long")

yrsamps <- tracks %>% mutate(year=year(DateTime)) %>% group_by(year) %>% summarise(n_distinct(ID))

# tracks %>% st_as_sf(coords = c("Longitude","Latitude"), crs = 4326, agr = "constant") %>%
#   mapview::mapview()

# project to custom LAEA projection
TD <- projectTracks(tracks, projType = "azim", custom=T)

rm(tracks)

# get candidate smoothing parameters 
# HVALS <- findScale(TD)
HVALS <- findScale(TD, scalesFPT = c(1:100))

HVALS

# h <- HVALS$scaleARS
h
h <- 7.5

# estimate individual core areas
b4 <- Sys.time()
KDE <- estSpaceUse(TD, scale=h, res=2.5, polyOut=T)
Sys.time() - b4

# saveRDS(KDE, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\KDE_h7.5.rds")
KDE <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\KDE_h7.5.rds")

mapKDE(KDE$UDPolygons)

n <- length(KDE)
n

# save core area map 
# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/white_storks/", analysis, "_indcores_AS1AS2_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6)


b4 <- Sys.time()
# rep <- repAssess(TD, KDE$KDE.Surface, iteration=200, bootTable = T)
# rep <- repAssess(TD, 
#                  KDE=readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\KDE_h7.5.rds"),
#                  iteration=200, 
#                  bootTable = T,
#                  nCores = 2)

Sys.time() - b4
# represent <- rep[[1]]$out
represent <- 95.7

# saveRDS(rep, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\rep_200its_h7.5km.rds")

rm(KDE)
rm(TD)

# potSite <- findSite(KDE = KDE$KDE.Surface, represent = represent, popSize = 26200) # popSize based on calc. from 'estimate_pop_size' spreadsheet
potSite <- findSite(KDE = KDE, represent = represent, levelUD=50, popSize = 26200)

# potSite <- findSite(KDE = KDE$KDE.Surface, represent = 96)

# saveRDS(potSite, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\potSite_h7.5.rds")

## Crop grid and convert to SF polygons (too big to do within findSite) ##
# fullgrid <- potSite
fullgrid <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\potSite_h7.5.rds")

# mapSite(fullgrid)

potSite <- fullgrid[fullgrid$potentialSite == TRUE, ]

bbox_poly <- as(extent(as.vector(t(bbox(potSite)))), "SpatialPolygons")
proj4string(bbox_poly) <- proj4string(potSite)

# Split spatially disparate areas into different sites # 
potSite <- brick(potSite)
potSite <- crop(potSite, bbox_poly)
potSite <- brick(addLayer(
  potSite, 
  clump(potSite$potentialSite == T, directions=8, gaps=F) # give unique IDs for different potSite sites
  ))

OUTMAP <- raster::aggregate(
  as(potSite, "SpatialPolygonsDataFrame"), 
  c('N_animals','N_IND','potentialSite','clumps')
)

# rm(potSite)

### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
potSite <- sf::st_as_sf(OUTMAP) %>%
  sf::st_union(by_feature = TRUE) %>%
  # smoothr::smooth(method = "densify") %>%
  sf::st_transform(4326) %>%
  arrange(.data$N_IND)
# rm(OUTMAP)

plot(potSite)

potSite %>% group_by(clumps) %>% summarise(
  max_at_site = max(N_animals),
  perc_region = max_at_site / 450000 * 100, # number of WHST in Europe 
  perc_global = max_at_site / 700000 * 100,  # number of WHST globally 
  perc_port   = max_at_site / 46027  * 100  # % of total Portuguse pop.
) %>% st_drop_geometry()

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add in some background maps for context 

library(ggmap)

ws_aggs <- potSite[potSite$N_animals > 0, ]
# ws_aggs <- potSite[potSite$potentialSite == TRUE, ]
# ws_aggs <- potSite[potSite$N_animals >= (0.1 * popsize), ]

xmin <- st_bbox(ws_aggs)[[1]] - 1.1
xmax <- st_bbox(ws_aggs)[[3]] + 1.5
ymin <- st_bbox(ws_aggs)[[2]] - 0.2
ymax <- st_bbox(ws_aggs)[[4]] + 0.2

gmap <- ggmap::get_map(location=c(xmin, ymin, xmax, ymax), zoom=9, maptype = "satellite")

ws_aggs <- st_transform(ws_aggs, crs = 3857)

# Use the function:
source("C:/Users/Martim Bill/Documents/R/source_scripts/ggmap_bbox.r")

gmap <- ggmap_bbox(gmap)

ws_map <- ggmap(gmap) + 
  coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  geom_sf(data = ws_aggs, inherit.aes = FALSE, aes(fill=N_animals), color=NA) +
  # geom_sf(data = ws_aggs, inherit.aes = FALSE, aes(fill=N_IND), color=NA) + 
  scale_fill_gradientn(colours=sf.colors(n=3), breaks=c(3000,5000,7000,9000,11000), labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  guides( 
    fill = guide_colorbar(
      barwidth  = 1.75,
      title = "N animals",
      barheight = 10)) +
  theme(
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=11, color="black"),
    axis.title=element_text(size=16),
    legend.direction = "vertical",
    legend.position=c(.9, .25),
    legend.title=element_text(size=15),
    legend.text = element_text(size = 13)) +
  ylab("Latitude") +
  xlab("Longitude")
ws_map

saveRDS(ws_map, "C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\WS_Site_ggmap.rds")

## SAVE ##

## all overlap
# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/white_storks/", analysis, "_N_INDs_AS1AS2_", "h", round(h), "_", "n",n, ".png"), ws_map, width = 9, height=9 )

## only potential Sites
# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/white_storks/", analysis, "_potSite_N_INDs_AS1AS2_c2.5_", "h", round(h), "_", "n",n, ".png"), ws_map, width = 8, height=10 ) # NINDS
ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/white_storks/", analysis, "_potSite_N_animals_AS1AS2_c2.5_", "h", round(h), "_", "n",n, ".png"), ws_map, width = 8, height=10 ) # N animals 

## Calc. area of potentialSites ##
site <- potSite %>% filter( potentialSite == T ) %>% group_by(clumps) %>% summarise(N_animals = max(N_animals)) 
plot(site)

## combine small ones near large ones ## 
site <- site %>% ungroup() %>% mutate(
  clumps = ifelse(clumps == 3, 2, ifelse(
    clumps == 9, 8, clumps
  )),
  clumps = as.numeric(as.factor(clumps))
) %>% group_by(clumps) %>% summarise(N_animals = max(N_animals))

plot(site)

site$area <- st_area(site) / 1000^2 # sq km

plot(site %>% filter(N_animals>8000))

potSite

