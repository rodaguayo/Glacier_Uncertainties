# Code for RGI and catchment preprocessing  ------------------------------------------------------
# Developed by Rodrigo Aguayo (2022-2023)

rm(list=ls())
cat("\014")  

setwd("/home/rooda/Dropbox/Patagonia/")

library("exactextractr")
library("rgrass")
library("terra")
library("stats")
library("dplyr")
library("sf")
sf_use_s2(FALSE)

setwd("/home/rooda/Dropbox/Patagonia/")

# 1. Delimitation of all basins -------------------------------------------------------------------
dem        <- rast("GIS South/dem_patagonia3f.tif") # DEM from NASADEM 90m
dem        <- aggregate(dem, fact=2, fun="mean")
dem        <- crop(dem, ext(c(-75, -67, -56, -40)))
dem        <- project(dem, "EPSG:32719")

initGRASS(gisBase = "/usr/lib/grass82/", home="/home/rooda/Dropbox/Patagonia/", 
          SG = dem, override = TRUE) 

write_RAST(dem,    vname = "dem_grass",   flags = c("overwrite", "o"))
write_RAST(is.na(dem),  vname = "depre_grass", flags = c("overwrite", "o")) # depression DEM -> water/ocean 

execGRASS("r.stream.extract", flags=c("overwrite"), 
          parameters = list(elevation="dem_grass", depression = "depre_grass", threshold = 200, 
                            direction = "fdir", stream_vector="stream_v", stream_raster="stream_r"))

execGRASS("r.stream.basins", flags=c("overwrite", "l"),
          parameters=list(direction="fdir", stream_rast = "stream_r", basins="basins"))

# 1.1 Attributes and subset by area and location --------------------------------------------------
all_basins      <- read_RAST("basins")
all_basins      <- as.polygons(all_basins)
all_basins      <- project(all_basins, "EPSG:4326")
all_basins$basin_area <- round(expanse(all_basins, unit="km"), 2)
all_basins$lat      <- as.data.frame(geom(centroids(all_basins)))$y
all_basins$lon      <- as.data.frame(geom(centroids(all_basins)))$x
all_basins <- subset(all_basins, (all_basins$lat < -41) & (all_basins$lat > -55.5))
all_basins <- subset(all_basins, (all_basins$lon < -70) | (all_basins$lat < -54))
all_basins <- subset(all_basins, !((all_basins$lon > -71) & (all_basins$lat > -50) & (all_basins$lat < -48)))
all_basins <- subset(all_basins, (all_basins$lon < -68.5))
all_basins <- subset(all_basins, all_basins$basin_area > 10)

# 1.2 Assign Zone and preliminary ID to all basins ------------------------------------------------
all_basins$Zone <- 0 # Initialization
all_basins$Zone <- ifelse(all_basins$lat > -43.4, 1, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -43.4 & all_basins$lat > -46,    2, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -46   & all_basins$lat > -47.8 & all_basins$lon > -73,   3, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -46   & all_basins$lat > -47.8 & all_basins$lon < -73,   4, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -47.8 & all_basins$lat > -49.4, 5, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -49.4 & all_basins$lat > -50.7, 6, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -50.7 & all_basins$lat > -52.1 , 7, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -52.1 & all_basins$lat > -54.1 , 8, all_basins$Zone)
all_basins$Zone <- ifelse(all_basins$lat < -54.1, 9, all_basins$Zone)
all_basins$ID   <- seq(0, nrow(all_basins)-1) # preliminary IDs
plot(all_basins, "Zone")

all_basins <- st_as_sf(all_basins)
all_basins <- st_transform(all_basins, 32718) # UTM 18S for the glaciers

# delete trash from GRASS (file12135term1k2 etc)
all_basins <- all_basins[,c("basin_area", "lat", "lon", "ID",  "Zone", "geometry")]

# 2. RGI glaciers: Selection and zone assignment --------------------------------------------------
RGI  <- list(st_read("GIS South/Glaciers/RGI6.shp"), st_read("GIS South/Glaciers/RGI7.shp"))

dem <- rast("GIS South/dem_patagonia3f.tif")
dem <- project(dem, "epsg:32718", method = "bilinear")
dem <- subst(dem, NA, 0)  0 # NAs to sea level (= 0)

for (i in 1:2) { 
  RGI_i  <- subset(RGI[[i]], RGI[[i]]$CenLat < -40.5)
  RGI_i  <- st_transform(RGI_i, 32718) # UTM 18S
  
  # terminus location
  RGI_ic <- extract(dem, as.lines(vect(RGI_i)), xy = TRUE)
  RGI_ic <- RGI_ic %>% group_by(ID) %>% slice_min(order_by = dem_patagonia3f)
  RGI_ic <- aggregate(RGI_ic, by = list(RGI_ic$ID), FUN = mean) # several options if terminus close to water
  RGI_ic <- st_as_sf(RGI_ic, coords = c("x","y"), crs = 32718)
  
  # intersection with basins
  RGI_ic <- st_join(RGI_ic["ID"], all_basins[,c("Zone", "ID", "geometry")], join = st_within) 
  RGI_i$Zone <- RGI_ic$Zone # Assign zone based on terminus location
  RGI_i$ID_basin <- RGI_ic$ID.y  # Assign catchment
  
  # centroid location if terminus is outside the basins
  RGI_ic <- st_centroid(RGI_i["geometry"]) 
  RGI_ic <- st_join(RGI_ic, all_basins[,c("Zone", "ID", "geometry")], join = st_within)
  RGI_i$Zone[is.na(RGI_i$Zone)] <- RGI_ic$Zone[is.na(RGI_i$Zone)]
  RGI_i$ID_basin[is.na(RGI_i$ID_basin)] <- RGI_ic$ID[is.na(RGI_i$ID_basin)]  # Assign catchment
  
  # select glaciers with a code for the Zone 
  RGI_i <- subset(RGI_i, RGI_i$Zone > 0) 
  RGI_i$O2Region <- 1 # useful to count
  RGI_i$area_km2 <- expanse(vect(RGI_i), unit="km")
  RGI[[i]] <- RGI_i
}    

# 3. Subset catchments if there are enough glaciers -----------------------------------------------

# if there at least one glacier in RGI6 or RGI7
RGI_groupby <- list(aggregate(st_drop_geometry(RGI[[1]][,c("O2Region", "area_km2")]), 
                          by = list(RGI[[1]]$ID_basin), FUN = "sum"),
                    aggregate(st_drop_geometry(RGI[[2]][,c("O2Region", "area_km2")]), 
                          by = list(RGI[[2]]$ID_basin), FUN = "sum"))
RGI_groupby <- merge(RGI_groupby[[1]], RGI_groupby[[2]], by = "Group.1", all=TRUE)
RGI_groupby[is.na(RGI_groupby)] <- 0
colnames(RGI_groupby) <- c("ID", "RGI6_ngla","RGI6_area","RGI7_ngla", "RGI7_area")

all_basins <- subset(all_basins, all_basins$ID %in% RGI_groupby$ID)

# if the glacier area in RGI6 or RGI7 is higher than 0.1%  
all_basins <- merge(all_basins, RGI_groupby, by = "ID", all=TRUE)
all_basins <- subset(all_basins, all_basins$RGI6_area/all_basins$basin_area > 0.001 | 
                                 all_basins$RGI7_area/all_basins$basin_area > 0.001)

# 4. Subset glaciers again ------------------------------------------------------------------------
RGI[[1]] <- subset(RGI[[1]], RGI[[1]]$ID_basin %in% all_basins$ID) # small effect
RGI[[2]] <- subset(RGI[[2]], RGI[[2]]$ID_basin %in% all_basins$ID)

# 5. Finish catchment file ------------------------------------------------------------------------
glacier_area <- vect("GIS South/Glaciers/RGI6.shp") # gridded area to compare just in case
glacier_area <- rasterize(glacier_area, project(dem, "EPSG:4326"), background = 0) * 100

all_basins <- st_transform(all_basins, 4326) # go back to WGS84
all_basins$RGI6_area_r <- extract(glacier_area, all_basins, "mean")[,2] 

writeVector(vect(all_basins), "GIS South/Basins_Patagonia_all.shp", overwrite=TRUE)
plot(all_basins) # visual check


# 6. Glacier attributes ---------------------------------------------------------------------------

## 6.1 Topography ----------------------------------------------------------------------------------
slope   <- terrain(dem, v="slope", neighbors=8, unit="degrees")  
aspect  <- terrain(dem,  v="aspect", neighbors=8, unit="degrees")  

for (i in 1:2) { 
  RGI[[i]]$slope_v2  <- exact_extract(slope,  RGI[[i]], "mean")
  RGI[[i]]$aspect_v2 <- exact_extract(aspect, RGI[[i]], "mean")
}

## 6.2 Volume data from Farinotti et al. and Millan et al. ----------------------------------------
vol_F19  <- rast("GIS South/Glaciers/Thickness_2019.tif")
vol_M22  <- rast("GIS South/Glaciers/Thickness_2022.tif")

### 6.2.1 Millan et al. 2022 (M22) ---------------------------------------------------------------

RGI6 <- RGI[[1]]
RGI7 <- RGI[[2]]

RGI6$vol_M22  <- exact_extract(vol_M22,  RGI6, "sum") * res(vol_M22)[1] * res(vol_M22)[2] * 1e-9
RGI7$vol_M22  <- exact_extract(vol_M22,  RGI7, "sum") * res(vol_M22)[1] * res(vol_M22)[2] * 1e-9
RGI6$vol_M22c <- round(exact_extract(not.na(vol_M22),  RGI6, "mean"), 3)
RGI7$vol_M22c <- round(exact_extract(not.na(vol_M22),  RGI7, "mean"), 3)
RGI6$vol_M22[RGI6$vol_M22c < 0.5] <- NA
RGI7$vol_M22[RGI7$vol_M22c < 0.5] <- NA

for (i in sort(unique(RGI6$Zone))) { # VAS based on M22 data 
  model6 <- lm(log(RGI6$vol_M22)[RGI6$Zone == i & RGI6$vol_M22c > 0.9] ~ log(RGI6$area_km2)[RGI6$Zone == i & RGI6$vol_M22c > 0.9] )
  model7 <- lm(log(RGI7$vol_M22)[RGI7$Zone == i & RGI7$vol_M22c > 0.9] ~ log(RGI7$area_km2)[RGI7$Zone == i & RGI7$vol_M22c > 0.9] )
  RGI6$vol_M22[RGI6$Zone == i & is.na(RGI6$vol_M22)] <- exp(coef(model6)[1]) * RGI6$area_km2[RGI6$Zone == i & is.na(RGI6$vol_M22)] ** coef(model6)[2]
  RGI7$vol_M22[RGI7$Zone == i & is.na(RGI7$vol_M22)] <- exp(coef(model7)[1]) * RGI7$area_km2[RGI7$Zone == i & is.na(RGI7$vol_M22)] ** coef(model7)[2]
}

### 6.2.2  Farinotti et al. 2019 (F19) ------------------------------------------------------------
RGI6$vol_F19  <- exact_extract(vol_F19,  RGI6, "sum") * res(vol_F19)[1] * res(vol_F19)[2] * 1e-9
RGI7$vol_F19  <- exact_extract(vol_F19,  RGI7, "sum") * res(vol_F19)[1] * res(vol_F19)[2] * 1e-9
RGI6$vol_F19c <- round(exact_extract(not.na(vol_F19),  RGI6, "mean"), 3)
RGI7$vol_F19c <- round(exact_extract(not.na(vol_F19),  RGI7, "mean"), 3)
RGI6$vol_F19[RGI6$vol_F19c < 0.5] <- NA
RGI7$vol_F19[RGI7$vol_F19c < 0.5] <- NA
RGI6$vol_F19[RGI6$vol_F19 == 0] <- 1e-6 # few glaciers with 0 volume (problem for log scale)
RGI7$vol_F19[RGI7$vol_F19 == 0] <- 1e-6

for (i in sort(unique(RGI6$Zone))) { # VAS based on F19 data
  model6 <- lm(log(RGI6$vol_F19)[RGI6$Zone == i & RGI6$vol_F19c > 0.9] ~ log(RGI6$area_km2)[RGI6$Zone == i & RGI6$vol_F19c > 0.9] )
  model7 <- lm(log(RGI7$vol_F19)[RGI7$Zone == i & RGI7$vol_F19c > 0.9] ~ log(RGI7$area_km2)[RGI7$Zone == i & RGI7$vol_F19c > 0.9] )
  RGI6$vol_F19[RGI6$Zone == i & is.na(RGI6$vol_F19)] <- exp(coef(model6)[1]) * RGI6$area_km2[RGI6$Zone == i & is.na(RGI6$vol_F19)] ** coef(model6)[2]
  RGI7$vol_F19[RGI7$Zone == i & is.na(RGI7$vol_F19)] <- exp(coef(model7)[1]) * RGI7$area_km2[RGI7$Zone == i & is.na(RGI7$vol_F19)] ** coef(model7)[2]
}

## 6.3 dmdtda: specific-mass change rate in meters water-equivalent per year ----------------------
dhdt_21   <- rast("GIS South/Glaciers/dhdt_2021.tif")
dhdt_e21  <- rast("GIS South/Glaciers/dhdt_error_2021.tif")
dhdt_oggm <- read.csv("MS2 Results/dhdt_origin_theia.csv")
dhdt_oggm <- subset(dhdt_oggm, dhdt_oggm$period == "2000-01-01_2020-01-01")
dhdt_oggm <- dhdt_oggm[dhdt_oggm$rgiid %in% RGI6$RGIId,] # same order

RGI6$dmdtda_21  <- exact_extract(dhdt_21,  RGI6, "mean") * 0.850 # dhdt to dmdtda
RGI7$dmdtda_21  <- exact_extract(dhdt_21,  RGI7, "mean") * 0.850 
RGI6$dmdtda_21c <- round(exact_extract(not.na(dhdt_21),  RGI6, "mean"), 3)
RGI7$dmdtda_21c <- round(exact_extract(not.na(dhdt_21),  RGI7, "mean"), 3)
RGI6$dmdtda_21[RGI6$dmdtda_21c < 0.9] <- NA
RGI7$dmdtda_21[RGI7$dmdtda_21c < 0.9] <- NA

RGI6$dmdtda_error_21  <- exact_extract(dhdt_e21,  RGI6, "mean") * 0.850 
RGI7$dmdtda_error_21  <- exact_extract(dhdt_e21,  RGI7, "mean") * 0.850 
RGI6$dmdtda_error_21[RGI6$dmdtda_21c < 0.9] <- NA
RGI7$dmdtda_error_21[RGI7$dmdtda_21c < 0.9] <- NA

## 6.3.1  Filling: Every glacier needs to have a dhdt (dmdadt) ------------------------------------
dmdtda_RGI6_mean <- sapply(split(RGI6, RGI6$Zone), function(d) weighted.mean(d$dmdtda_21, w = d$area_km2, na.rm = T)) # area-weighted is better :)
dmdtda_RGI7_mean <- sapply(split(RGI7, RGI7$Zone), function(d) weighted.mean(d$dmdtda_21, w = d$area_km2, na.rm = T))

for (i in sort(unique(RGI6$Zone))) {
  RGI6$dmdtda_21[RGI6$Zone == i & is.na(RGI6$dmdtda_21)] <- dmdtda_RGI6_mean[i]
  RGI7$dmdtda_21[RGI7$Zone == i & is.na(RGI7$dmdtda_21)] <- dmdtda_RGI7_mean[i]
}

# 7. Save glacier file ----------------------------------------------------------------------------

# go back to WGS84 (RGI format)
RGI6 <- st_transform(RGI6, 4326) 
RGI7 <- st_transform(RGI7, 4326)

# area is area_km2
RGI7$area <- NULL

# source year in RGI format
RGI7$src_date <- paste0(substr(RGI7$src_date, 0,4), substr(RGI7$src_date, 6,7), substr(RGI7$src_date, 9,10)) 

# save using terra (problems in sf)
writeVector(vect(RGI6), "GIS South/Glaciers/RGI6_v2.shp", overwrite=TRUE)
writeVector(vect(RGI7), "GIS South/Glaciers/RGI7_v2.shp", overwrite=TRUE)



