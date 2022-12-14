# Code for preprocess RGI inventories
# Developed by Rodrigo Aguayo (2022-2023)

rm(list=ls())
cat("\014")  

setwd("/home/rooda/Dropbox/Patagonia/")

library("exactextractr")
library("terra")
library("stats")
library("dplyr")
library("sf")
sf_use_s2(FALSE)

# 0. Data ---------------------------------------------------------------------

## DEM from NASADEM 90m
dem     <- project(rast("GIS South/dem_patagonia3f.tif"), "epsg:32718", method = "bilinear")
dem[is.na(dem)] <- 0 # NAs to sea level (= 0)

## basins in Patagonia
basins <- st_read("GIS South/Basins_Patagonia_all.shp")
basins <- st_transform(basins, 32718) # UTM 18S
basins <- basins[,c("Zone", "ID", "geometry")]

# 1. RGI glaciers: Selection and assignment -----------------------------------

## 1.1 RGI6 outlines ----------------------------------------------------------
RGI6  <- st_read("GIS South/Glaciers/RGI6.shp")
RGI6  <- subset(RGI6, RGI6$CenLat < -40.5)
RGI6  <- st_transform(RGI6, 32718) # UTM 18S

RGI6c <- extract(dem, as.lines(vect(RGI6)), xy = TRUE)
RGI6c <- RGI6c %>% group_by(ID) %>% slice_min(order_by = dem_patagonia3f)
RGI6c <- aggregate(RGI6c, by = list(RGI6c$ID), FUN = mean) # not ideal 
RGI6c <- st_as_sf(RGI6c, coords = c("x","y"), crs = 32718)
RGI6c <- st_join(RGI6c["ID"], basins, join = st_within)
RGI6$Zone <- RGI6c$Zone # Assign zone
RGI6$ID_basin <- RGI6c$ID.y # Assign catchment

RGI6c <- st_centroid(RGI6["geometry"]) # fill with the centroid 
RGI6c <- st_join(RGI6c, basins, join = st_within)

RGI6$Zone[is.na(RGI6$Zone)] <- RGI6c$Zone[is.na(RGI6$Zone)]
RGI6 <- subset(RGI6, RGI6$Zone > 0) # only glaciers with a code
RGI6$area_km2  <- expanse(vect(RGI6), unit="km")

## 1.2 RGI7 outlines ----------------------------------------------------------
RGI7 <- st_read("GIS South/Glaciers/RGI7.shp")
RGI7  <- subset(RGI7, RGI7$CenLat < -40.5)
RGI7 <- st_transform(RGI7, 32718) # UTM 18S

RGI7c <- extract(dem, as.lines(vect(RGI7)), xy = TRUE)
RGI7c <- RGI7c %>% group_by(ID) %>% slice_min(order_by = dem_patagonia3f)
RGI7c <- aggregate(RGI7c, by = list(RGI7c$ID), FUN = mean) # not ideal 
RGI7c <- st_as_sf(RGI7c, coords = c("x","y"), crs = 32718)
RGI7c <- st_join(RGI7c["ID"], basins, join = st_within)
RGI7$Zone <- RGI7c$Zone  # Assign zone
RGI7$ID_basin <- RGI7c$ID.y  # Assign catchment

RGI7c <- st_centroid(RGI7["geometry"]) # fill with the centroid 
RGI7c <- st_join(RGI7c, basins, join = st_within)

RGI7$Zone[is.na(RGI7$Zone)] <- RGI7c$Zone[is.na(RGI7$Zone)]
RGI7 <- subset(RGI7, RGI7$Zone > 0) # only glaciers with a code
RGI7$area_km2  <- expanse(vect(RGI7), unit="km")

# 2. Topographic attributes ---------------------------------------------------
slope   <- terrain(dem, v="slope", neighbors=8, unit="degrees")  
aspect  <- terrain(dem,  v="aspect", neighbors=8, unit="degrees")  
RGI6$slope_v2  <- exact_extract(slope,  RGI6, "mean")
RGI7$slope_v2  <- exact_extract(slope,  RGI7, "mean")
RGI6$aspect_v2 <- exact_extract(aspect,  RGI6, "mean")
RGI7$aspect_v2 <- exact_extract(aspect,  RGI7, "mean")

# 3. Volume data from Farinotti et al. and Millan et al. ----------------------
vol_F19  <- rast("GIS South/Glaciers/Thickness_2019.tif")
vol_M22  <- rast("GIS South/Glaciers/Thickness_2022.tif")

## 3.1 Millan et al. 2022 (M22) -----------------------------------------------
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


## 3.2 Farinotti et al. 2019 (F19) --------------------------------------------
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

# 4. dmdtda: specific-mass change rate in meters water-equivalent per year ----
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

# 4.1 Filling: Every glacier needs to have a dhdt (dmdadt) --------------------
dmdtda_RGI6_mean <- sapply(split(RGI6, RGI6$Zone), function(d) weighted.mean(d$dmdtda_21, w = d$area_km2, na.rm = T)) # area-weighted is better :)
dmdtda_RGI7_mean <- sapply(split(RGI7, RGI7$Zone), function(d) weighted.mean(d$dmdtda_21, w = d$area_km2, na.rm = T))

for (i in sort(unique(RGI6$Zone))) {
  RGI6$dmdtda_21[RGI6$Zone == i & is.na(RGI6$dmdtda_21)] <- dmdtda_RGI6_mean[i]
  RGI7$dmdtda_21[RGI7$Zone == i & is.na(RGI7$dmdtda_21)] <- dmdtda_RGI7_mean[i]
}

# 5. Save file ----------------------------------------------------------------

# go back to wgs84 (RGI format)
RGI6 <- st_transform(RGI6, 4326) 
RGI7 <- st_transform(RGI7, 4326)

# area is area_km2
RGI7$area <- NULL

# source year in RGI format
RGI7$src_date <- paste0(substr(RGI7$src_date, 0,4), substr(RGI7$src_date, 6,7), substr(RGI7$src_date, 9,10)) 

# save using terra (problems in sf)
writeVector(vect(RGI6), "GIS South/Glaciers/RGI6_v2.shp", overwrite=TRUE)
writeVector(vect(RGI7), "GIS South/Glaciers/RGI7_v2.shp", overwrite=TRUE)
