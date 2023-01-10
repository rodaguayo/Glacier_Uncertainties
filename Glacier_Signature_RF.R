# Code for XXXX ---------------------------------------------------------
# Developed by Rodrigo Aguayo (2022-2023)

rm(list=ls())
cat("\014")  

library("randomForest")
setwd("/home/rooda/Dropbox/Patagonia/MS2 Results/")

# data
list_data <- list(Timing      = read.csv("data_basin_peak_water_year.csv"),
                  Rate_change = read.csv("data_basin_rate_change.csv"),
                  Duration    = read.csv("data_basin_duration.csv"), 
                  Magnitude   = read.csv("data_basin_peak_water_value.csv"), 
                  Frequency   = read.csv("data_basin_freq.csv"))

for (i in 1:5){ # at least one glacier in RGI6 and RGI7
  dataset <- list_data[[i]]
  dataset <- dataset[ , colSums(is.na(dataset)) == 0]
  list_data[[i]] <- dataset
}

basins <- colnames(list_data$Rate_change)[-7:0]
variables <- c("Outline", "Climate", "Volume",  "GGM", "SSP", "BCM")
rf_importance <- data.frame(matrix(ncol = 8, nrow = 0))

for (i in 1:5){ # for each metric
  
  # select dataset of glacier signature
  dataset <- list_data[[i]]
  
  # building the dataset
  dataset$X <- NULL
  dataset$Outline <- as.factor(dataset$Outline)
  dataset$Volume  <- as.factor(dataset$Volume)
  dataset$Climate <- as.factor(dataset$Climate)
  dataset$GGM     <- as.factor(dataset$GGM)
  dataset$SSP     <- as.factor(dataset$SSP)
  dataset$BCM     <- as.factor(dataset$BCM)
  
  for (j in 1:length(basins)){ # for each basin
    
    # select zone
    basin       <- basins[j]
    dataset_i  <- dataset[,c(variables, basin)]
    
    # random forest regression model
    rf_model      <- randomForest(as.formula(paste0(basin, " ~.")), data = dataset_i, ntree = 500, importance = TRUE)
    rf_importance <- rbind(rf_importance, c(basin, names(list_data)[i], importance(rf_model)[,1]))
    print(paste0("Metric ", i," - ", "Basin ", j))
    }
}

colnames(rf_importance) <- c("Zone", "Glacier_Signature", variables)
write.csv(rf_importance, "RF_importance_basin.csv")

