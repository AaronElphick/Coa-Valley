#https://www.statisticshowto.com/mann-kendall-trend-test/
#The minimum number of recommended measurements is therefore at least 8 to 10.
install.packages("Kendall")

library(Kendall)
library(sp)
library(tidyverse)
library(raster) 
library(terra)
library(geodata)
library(terra)
library(sf)
library(ggplot2)
library(remotes)
library(ecbtools)
library(randomForest)
library(caret)
library(splitTools)
library(ranger)
library(rgdal)
remotes::install_github("fdetsch/MODIS")
library(MODIS)
?matplot

setwd("~/Ecology and Data Science/Dissertation")

MODISoptions(MODISserverOrder = "LPDAAC", quiet = FALSE, localArcPath = "OutputDestinationFolder/MODIS")
getHdf("MOD13Q1", collection = "006", tileH = 17, tileV = 4, 
       begin = "2001.01.01", end= "2022.12.31")

# input filenames
inf <- list.files("OutputDestinationFolder/MODIS/MOD13Q1.006", full.names=TRUE, all.files = TRUE, recursive=TRUE)

# GCV Shapefile
library(raster)
shape2 <- read_sf(file.path("R Stuff/Shapefiles/ELP-WIarea_PTM06.shp"))
shp <- st_transform(shape2, crs(ndvi_stack))
shp2<-st_zm(shp)


#Inserting the raster information (will take a long time)
allrasters<- stack(inf)

#Load the Mann Kendall Test Function

MKraster <- function(rasterstack, type=c("trend","pval","both")){
  
  # Values for (layers, ncell, ncol, nrow, method, crs, extent) come straight from the input raster stack
  # e.g. nlayers(rasterstack), ncell(rasterstack)... etc.
  print(paste("Start MKraster:",Sys.time()))
  print("Loading parameters")
  layers=nlayers(rasterstack);ncell=ncell(rasterstack);
  ncol=ncol(rasterstack);nrow=nrow(rasterstack);crs=crs(rasterstack);
  extent=extent(rasterstack);pb = txtProgressBar(min = 0, max = ncell, initial = 0)
  print("Done loading parameters")
  mtrx <- as.matrix(rasterstack,ncol=layers)
  empt <- matrix(nrow=ncell, ncol=2)
  
  print("Initiating loop operation")
  if (type == "trend"){
    
    for (i in 1:length(mtrx[,1])){
      if (all(is.na(mtrx[i,]))){ 
        empt[i,1] <- NA 
      } else 
        if (sum(!is.na(mtrx[i,])) < 4){
          empt[i,1] <- NA 
        } else 
          empt[i,1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
    }
    
    print("Creating empty raster")
    trend <- raster(nrows=nrow,ncols=ncol,crs=crs)
    extent(trend) <- extent
    print("Populating trend raster")
    values(trend) <- empt[,1]
    print(paste("Ending MKraster on",Sys.time()))
    trend
  } 
  else
    if (type == "pval"){
      
      for (i in 1:length(mtrx[,1])){
        if (all(is.na(mtrx[i,]))){ 
          empt[i,1] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,])) < 4){
            empt[i,1] <- NA 
          } else 
            empt[i,1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
      }
      
      pval <- raster(nrows=nrow,ncols=ncol,crs=crs)
      extent(pval) <- extent
      print("Populating significance raster")
      values(pval) <- empt[,1]
      print(paste("Ending MKraster on",Sys.time()))
      pval
    }
  else
    if (type == "both"){
      
      for (i in 1:length(mtrx[,1])){
        if (all(is.na(mtrx[i,]))){ 
          empt[i,1] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,])) < 4){
            empt[i,1] <- NA 
          } else 
            tryCatch({
              empt[i,1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
              empt[i,2] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
            })
      }
      
      tr <- raster(nrows=nrow,ncols=ncol,crs=crs)
      pv <- raster(nrows=nrow,ncols=ncol,crs=crs)
      print("Populating raster brick")
      values(tr) <- empt[,1]
      values(pv) <- empt[,2]
      brk <- brick(tr,pv)
      extent(brk) <- extent
      names(brk) <- c("trend","p.value")
      print(paste("Ending MKraster on",Sys.time()))
      brk
    }
}

#Selecting the NDVI bands and  masking them
seq_test<-seq(from=1, to=6348, by=12)
ndvi_stack<-allrasters[[c(seq_test)]]
ndvi_stack<-crop(ndvi_stack, shp)
ndvi_stack_mask<- mask(ndvi_stack, shp)

values(ndvi_stack_mask)[values(ndvi_stack_mask) < 0] = 0

#Number represents year
#Firstly splitting the modis files into years from 2001 and 2022 (1 and 22)

for (i in 1:22) {
  variable_name <- paste0("year_", i)
  assign(variable_name, (i - 1) * 23 + 21 : (i * 23))
}


for (i in 1:22) {
  variable = paste0("test", i, sep = "")
  stack_no = get(paste0("year_", i, sep = ""))
  stack_no = ndvi_stack_mask[[stack_no]]
  assign(variable, stack_no)

#Calculating yearly values i.e the max, min, sum and RREL.

for (i in 1:22) {
  variable = paste0("max_", i, sep = "")
  temp_data = get(paste0("test", i, sep = ""))
  temp_data = calc(temp_data, function(x){max(x)})
  assign(variable, temp_data)
}

for (i in 1:22) {
  variable = paste0("min_", i, sep = "")
  temp_data = get(paste0("test", i, sep = ""))
  temp_data = calc(temp_data, function(x){min(x)})
  assign(variable, temp_data)
}

for (i in 1:22) {
  variable = paste0("sum_", i, sep = "")
  temp_data = get(paste0("test", i, sep = ""))
  temp_data = calc(temp_data, function(x){sum(x)})
  assign(variable, temp_data)
}

for (i in 1:22) {
  variable = paste0("mean_", i, sep = "")
  temp_data = get(paste0("test", i, sep = ""))
  temp_data = calc(temp_data, function(x){mean(x)})
  assign(variable, temp_data)
}


for (i in 1:22) {
  variable <- paste0("RREL_", i)
  temp_data <- get(paste0("test", i))
  temp_data1 <- calc(temp_data, function(x) {min(x)})
  temp_data2 <- calc(temp_data, function(x) {max(x)})
  temp_data3 <- calc(temp_data, function(x) {mean(x)})
  temp_data <- (temp_data2 - temp_data1) / temp_data3
  assign(variable, temp_data)
}


min_stack<- stack(min_1, min_2, min_3, min_4, min_5, min_6, 
                  min_7, min_8, min_9, min_10, min_11, min_12,
                  min_13, min_14, min_15, min_16, min_17, min_18,
                  min_19, min_, min_21, min_22)

max_stack<-stack(max_1, max_2, max_3, max_4, max_5, max_6, 
                 max_7, max_8, max_9, max_10, max_11, max_12,
                 max_13, max_14, max_15, max_16, max_17, max_18,
                 max_19, max_, max_21, max_22)

sum_stack<-stack(sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, 
                 sum_7, sum_8, sum_9, sum_10, sum_11, sum_12,
                 sum_13, sum_14, sum_15, sum_16, sum_17, sum_18,
                 sum_19, sum_, sum_21, sum_22)

RREL_stack<-stack(RREL_1, RREL_2, RREL_3, RREL_4, RREL_5, RREL_6, 
                  RREL_7, RREL_8, RREL_9, RREL_10, RREL_11, RREL_12,
                  RREL_13, RREL_14, RREL_15, RREL_16, RREL_17, RREL_18,
                  RREL_19, RREL_, RREL_21, RREL_22)

mean_stack<- stack(mean_1, mean_2, mean_3, mean_4, mean_5, mean_6, 
                   mean_7, mean_8, mean_9, mean_10, mean_11, mean_12,
                   mean_13, mean_14, mean_15, mean_16, mean_17, mean_18,
                   mean_19, mean_, mean_21, mean_22)

#Applying the man kendall test

min_Kendall = MKraster(rasterstack = min_stack, type = "both")

max_Kendall = MKraster(rasterstack = max_stack, type = "both")

sum_Kendall = MKraster(rasterstack = sum_stack, type = "both")

RREL_Kendall = MKraster(rasterstack = RREL_stack, type = "both")

#Write out the results
writeRaster(RREL_Kendall$trend, filename = "OutputDestinationFolder/Final/RREL_NDVI_Trend.TIF", options=c('TFW=YES'))
writeRaster(RREL_Kendall$p.value, filename = "OutputDestinationFolder/Final/RREL_NDVI_PValue.TIF", options=c('TFW=YES'))
writeRaster(allrasters_c[[1]], filename = "OutputDestinationFolder/Final/test#3.TIF", options=c('TFW=YES'))
writeRaster(RREL_mask[[1]], filename = "OutputDestinationFolder/Final/RREL_Test.TIF", options=c('TFW=YES'))

