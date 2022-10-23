# Code for basin delimitation  -------------------------------------------------------------------
# Developed by Rodrigo Aguayo (2020-2022)

rm(list=ls())
cat("\014")  

library("rgrass")
library("terra")
library("sf")
setwd("/home/rooda/Dropbox/Patagonia/")

# Data: Station locations and DEM (Initializate GRASS and files)--------------------------------
dem        <- rast("GIS South/dem_patagonia3f.tif")
dem        <- aggregate(dem, fact=2, fun="mean")
dem        <- crop(dem, ext(c(-75, -67, -56, -40)))
dem        <- project(dem, "EPSG:32719")
depre      <- is.na(dem) # depression DEM -> water/ocean -> use less RAM

initGRASS(gisBase = "/usr/lib/grass82/", home="/home/rooda/Dropbox/Patagonia/", 
          SG = dem, override = TRUE) 

write_RAST(dem,    vname = "dem_grass",   flags = c("overwrite", "o"))
write_RAST(depre,  vname = "depre_grass", flags = c("overwrite", "o"))

execGRASS("r.stream.extract", flags=c("overwrite"), 
          parameters = list(elevation="dem_grass", depression = "depre_grass", threshold=200, 
                            direction = "fdir", stream_vector="stream_v", stream_raster="stream_r"))

execGRASS("r.stream.basins", flags=c("overwrite", "l"),
          parameters=list(direction="fdir", stream_rast = "stream_r", basins="basins"))

all_basins      <- read_RAST("basins")
all_basins      <- as.polygons(all_basins)
all_basins      <- project(all_basins, "EPSG:4326")
all_basins$area <- round(expanse(all_basins, unit="km"), 2)
all_basins$lat  <- as.data.frame(geom(centroids(all_basins)))$y
all_basins$lon  <- as.data.frame(geom(centroids(all_basins)))$x
all_basins <- subset(all_basins, (all_basins$lat < -41) & (all_basins$lat > -55.5))
all_basins <- subset(all_basins, (all_basins$lon < -70) | (all_basins$lat < -54))
all_basins <- subset(all_basins, !((all_basins$lon > -71) & (all_basins$lat > -50) & (all_basins$lat < -48)))
all_basins <- subset(all_basins, (all_basins$lon < -68.5))
all_basins <- subset(all_basins, all_basins$area > 10)
plot(all_basins)

glaciers6 <- vect("GIS South/Glaciers/RGI6_v2.shp")
glaciers7 <- vect("GIS South/Glaciers/RGI7_v2.shp")
glaciers6 <- rasterize(glaciers6, project(dem, "EPSG:4326"), background = 0) * 100
glaciers7 <- rasterize(glaciers7, project(dem, "EPSG:4326"), background = 0) * 100
glaciers  <- max(glaciers6, glaciers7)

all_basins$glacier_cover <- extract(glaciers, all_basins, "mean")[,2] 
all_basins <- subset(all_basins, all_basins$glacier_cover > 0.1) # 0.1%

all_basins$file16212f2b9905 <- NULL
all_basins$file1cf0e63335d59 <- NULL
all_basins$ID <- seq(1, nrow(all_basins))

# Assing code to all basins
all_basins <- st_as_sf(all_basins)
all_basins$Zone <- 0 #Initilizate


all_basins$Zone <- ifelse(all_basins$lat > -43.4, 1, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -43.4 & all_basins$lat > -46,    2, all_basins$Zone)

all_basins$Zone <- ifelse(all_basins$lat < -46   & all_basins$lat > -47.8 & all_basins$lon > -73,   3, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -46   & all_basins$lat > -47.8 & all_basins$lon < -73,   4, all_basins$Zone)

all_basins$Zone <- ifelse(all_basins$lat < -47.8 & all_basins$lat > -49.4, 5, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -49.4 & all_basins$lat > -50.7, 6, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -50.7 & all_basins$lat > -52.1 , 7, all_basins$Zone)

all_basins$Zone <- ifelse(all_basins$lat < -52.1 & all_basins$lat > -54.1 , 8, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -54.1, 9, all_basins$Zone)

glaciers6 <- vect("GIS South/Glaciers/RGI6_v2.shp")
all_basins$glacier_n <-lengths(sf::st_intersects(all_basins, st_as_sf(centroids(glaciers6))   ))


writeVector(vect(all_basins), "GIS South/Basins_Patagonia_all.shp", overwrite=TRUE)
plot(all_basins)
