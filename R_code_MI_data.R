library(ggplot2)
library(spNNGP) #Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes

data("MI_TSCA")
head(MI_TSCA)
summary(MI_TSCA)
# https://rdrr.io/cran/spNNGP/man/MI_TSCA.rda.html
#Covariates included minimum winter temperature (MIN), maximum summer temperature (MAX),
#total precipitation in the coldest quarter of the year (WIP), 
# total precipitation in the warmest quarter of the year (SUP), 
# annual actual evapotranspiration (AET) and annual climatic water deficit (DEF).
# Spatial coordinates are recorded in Longitude (long) and latitude (lat) which are Albers Equal Area 
# A data frame containing 17,743 rows and 9 columns.

ggplot(data=MI_TSCA,aes(x=lat,y=long,color = WIP))+geom_point() +
  scale_color_viridis_c(option="magma")
# colorspace examples
ggplot(data=MI_TSCA,aes(x=lat,y=long,color =SUP))+geom_point() +
  scale_color_viridis_c(option="magma")

cor(MI_TSCA[,4:9])
