# Coa-Valley
Code for my MSc Ecology and Data Science Dissertation Project

The Landsat_Code file includes the code I used to produce my RandomForest classifications of the Greater Coa Valley, Portugal. The shapefiles imported are labelled land cover classifications imported from QGIS and ground truthed using google earth imagery data. This classification is repeated 20 times in a for loop and automatically saved to an output file so that the classification is repeatable and robust.

The Modis_Code file uses every 16 day MODIS NDVI dataset between 2001 and 2022 (using complete years) in order to produce a Man Kendall test on four indices including the Maximum, Minimum, Mean and RREL NDVI values. Which produces a significance test where significance is assumed at <0.05, and therefore indicates whether primary productivity (indicative from NDVI values) has changed significantly or insignificantly, and in a positive or negative direction.
