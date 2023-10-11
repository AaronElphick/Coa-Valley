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


setwd("~/Ecology and Data Science/Dissertation/R Stuff")

#2022 raster information
blue = raster("Landsat/Landsat 2022 July/LC08_L2SP_203032_20220708_20220721_02_T1_SR_B2.TIF")
green = raster("Landsat/Landsat 2022 July/LC08_L2SP_203032_20220708_20220721_02_T1_SR_B3.TIF")
red = raster("Landsat/Landsat 2022 July/LC08_L2SP_203032_20220708_20220721_02_T1_SR_B4.TIF")
near.infrared= raster("Landsat/Landsat 2022 July/LC08_L2SP_203032_20220708_20220721_02_T1_SR_B5.TIF")
SWIR = raster("Landsat/Landsat 2022 July/LC08_L2SP_203032_20220708_20220721_02_T1_SR_B6.TIF")

shape <- read_sf(file.path("Shapefiles/ELP-WIarea_PTM06.shp"))

shp <- st_transform(shape, crs(blue))

shp<-st_zm(shp)

setwd("~/Ecology and Data Science/Dissertation/QGIS")

#Additional Indices calculated in QGIS
NDVI = raster("Outputs/2022_Landsat_NDVI.TIF")
NBR = raster("Outputs/Landsat_NBR.TIF")
MSAVI2 = raster("Outputs/MSAVI2.TIF")

#Load in your training data
shapetest<- read_sf(file.path("Outputs/Land Cover/Testing class/Test_2006#43.shp"))
shapetest <- st_transform(shapetest, crs(blue))
shapetest<-as(shapetest, 'Spatial')
shapetest <- setNames(shapetest, c("Class", "id"))

#Stacking the relevant band information.

Landsat <- stack(blue, green, red, near.infrared, SWIR, NDVI, NBR, MSAVI2)

Landsat<- crop(Landsat, shp)

#Producing class points out of the shapefiles
GetClassPoints<-function(data, polygons, MaxPointsPerPolygon=10, StratifyingCellSize=2){
  
  #get stratiefied poitns from polygons
  pts_list<-list()
  for (i in 1:length(polygons)){
    pts <- try(spsample(polygons[i,], type = "stratified", cellsize=StratifyingCellSize))
    if(class(pts)=="try-error") {return(print("StratifyingCellSize is too large for smallest polygon - there is some randomness in this process"))}
    pts<-sample(pts,MaxPointsPerPolygon, replace=TRUE)
    pts$class <- rep(polygons[i,]$Class, length(pts))
    pts_list<-c(pts_list, pts)
  }
  allpts <- do.call("rbind", pts_list)
  
  #extract values from raster based on raster cell each point falls within
  trainingvals <- raster::extract(data, y=allpts, cellnumbers=TRUE, method="simple")
  trainingvals <- data.frame(response = allpts$class, trainingvals)
  
  # remove raster cells that are selected multiple times
  if (any(duplicated(trainingvals$cells))) {
    print(paste0(sum(duplicated(trainingvals$cells)), " duplicated cells removed"))
    allpts<-allpts[!duplicated(trainingvals$cells),]
    trainingvals <- trainingvals[!duplicated(trainingvals$cells), -2]
  }

  #ensure class is a factor
  trainingvals$response<-as.factor(trainingvals$response)
  
  return( list( "pointVals"=trainingvals,
                "NumberCellsPerCategory" = table(trainingvals$response),
                "points"=allpts))
}

#Producing Training test split
set.seed(13)
inds <- partition(shapetest$Class, p = c(train = 0.7, test = 0.3))
train<- shapetest[inds$train, ]
test<- shapetest[inds$test, ]

train_pts<- GetClassPoints(Landsat, train, MaxPointsPerPolygon=10, StratifyingCellSize=2)
test_pts<- GetClassPoints(Landsat, test, MaxPointsPerPolygon=10, StratifyingCellSize=2)

#Random Forest Model selection
mod <- randomForest(response~., data=train_pts[["pointVals"]], na.action=na.omit, ntree=500, confusion=TRUE)

modsc <- predict(Landsat, mod)
modsc_unmerged<- predict(Landsat, mod)

Visual<-mask(crop(modsc,shp),shp)

plot(Visual)

obs<-test_pts[["points"]]$class %>% as.factor()

pred <- raster::extract(modsc, test_pts[["points"]], cellnumbers = TRUE)
preds<-pred[,"layer"] %>% as.factor()
levels(preds)<- c("Agriculture", "Conifer", "Grass", "Mixed", "Rocky", "Scrub", "Urban", "Water")

confMat<-confusionMatrix(obs, reference = preds)

setwd("~/Ecology and Data Science/Dissertation/R Stuff")

shape3 <- read_sf(file.path("Shapefiles/Ermo-aguias-limits_02.2022.shp"))
shape4 <-read_sf(file.path("Shapefiles/Paul-ToirΣes.shp"))
shape5 <-read_sf(file.path("Shapefiles/Vale-Carapito.shp"))

shp5<- st_transform(shape5, crs(blue))

shp5<-st_zm(shp5)

#Repeating the classification 

setwd("~/Ecology and Data Science/Dissertation/QGIS/Class Outputs_2022")

# where we will save our outputs

output_loc = "./Class Outputs_2022/"




# run bootstrap prediction 40 times, i is used as the seed selection

for(i in 41:60){
  
  
  
  # set a different seed for each
  
  
  random_seed = i
  
  
  set.seed(random_seed)
  
  
  
  # run_name
  
  run_name = paste("model_2022", random_seed, i, sep="_")
  test_name = paste("Test_Pts_2022", random_seed, i, sep="_")
  train_name = paste("Train_Pts_2022", random_seed, i, sep="_")
  table_train = paste("Train_Table_2022", random_seed, i, sep="_")
  table_test = paste("Test_Table_2022", random_seed, i, sep="_")
  conf_table = paste("Conf_Table_2022", random_seed, i, sep="_")
  conf_stats = paste("Conf_Stats_2022", random_seed, i, sep="_")
  conf_overall = paste("Conf_Overall_2022", random_seed, i, sep="_")
  # reporting
  
  print(paste("Running for", run_name, sep=" "))
  

  
  # subset the data to a random subset of 70% for training
  
  inds_i <- partition(shapetest$Class, p = c(train = 0.7, test = 0.3))
  train_i<- shapetest[inds_i$train, ]
  test_i<- shapetest[inds_i$test, ]
  
  train_pts_i<- GetClassPoints(Landsat, train_i, MaxPointsPerPolygon=10, StratifyingCellSize=2)
  test_pts_i<- GetClassPoints(Landsat, test_i, MaxPointsPerPolygon=10, StratifyingCellSize=2)
  
  
  # RandomForest
  
  mod_i <- randomForest(response~., data=train_pts_i[["pointVals"]], na.action=na.omit, ntree=500, confusion=TRUE)
  
  
  # predict on the new data
  
  modsc_i <- predict(Landsat, mod_i)
  
  obs_i<-test_pts_i[["points"]]$class %>% as.factor()
  
  pred_i <- raster::extract(modsc, test_pts_i[["points"]], cellnumbers = TRUE)
  preds_i<-pred_i[,"layer"] %>% as.factor()
  levels(preds_i)<- c("Agriculture", "Conifer", "Grass", "Mixed", "Rocky", "Scrub", "Urban", "Water")
  
  #Confusion Matrix
  
  confMat_i<-confusionMatrix(obs_i, reference = preds_i)
  
  # save preds to output location
  
  capture.output(confMat_i, append = TRUE)
  
  file_name = paste(run_name, ".TIF", sep="")
  
  table_train_name = paste(table_train, ".csv", sep="")
  
  table_test_name = paste(table_test, ".csv", sep="")
  
  conf_table_name = paste(conf_table, ".csv", sep="")
  
  conf_stats_name = paste(conf_stats, ".csv", sep="")
  
  conf_overall_name = paste(conf_overall, ".csv", sep="")
  
  
  writeRaster(modsc_i, filename = file_name, options=c('TFW=YES'))
  
  writeOGR(train_pts_i$points, dsn = '.' , layer = train_name, driver = "ESRI Shapefile")
  
  writeOGR(test_pts_i$points, dsn = '.' , layer = test_name, driver = "ESRI Shapefile")
  
  write.csv(train_pts_i$NumberCellsPerCategory, table_train_name)
  
  write.csv(test_pts_i$NumberCellsPerCategory, table_test_name)
  
  write.csv(confMat_i$table, conf_table_name)
  
  write.table(confMat_i$byClass, conf_stats_name)
  
  write.table(confMat_i$overall, conf_overall_name)
  
} # end loop



