library(ggplot2)
library(spNNGP)
library(schnellerGP)
library(Matrix)
data("MI_TSCA")
head(MI_TSCA)


library(splines)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(pracma)

set.seed(12345)
ssize = 1000  #Number of points to subsample from Y for inference

NCOMIX=5
M=5
x = cbind(MI_TSCA$lat,MI_TSCA$long)
x[,1] = (x[,1]-min(x[,1]))/max(x[,1])
x[,2] = (x[,2]-min(x[,2]))/max(x[,2])

indxsamp = sample(1:dim(x)[1],ssize) #Randomly sample indices of training sample

y1 = MI_TSCA$SUP[indxsamp]
y2 = MI_TSCA$WIP[indxsamp]
y  = cbind(y1,y2)
x = x[indxsamp,]

x[,1] = (x[,1]-min(x[,1]))/max(x[,1])
x[,2] = (x[,2]-min(x[,2]))/max(x[,2])
####################################
#
#
# S is the lat and lon
# D is the index of either Y1 or Y2

##################################
#Determine the size of the hold out data 
#set for predictive process
trainSize = 400
lengthSim = length(y1) - trainSize # number of objects to remove from 
                                   # the datasets below. Note 1000 is hardcoded
                                   # as the number of unique values of s \in S
                                   # change accordingly given the dataset
                                   #################################

######################################################################
# Get unique rows to train data on
# Get unique rows for hold out
######################################################################
set.seed(8675309)
temp_idx = 1:(trainSize + lengthSim)
train_idx = sample(temp_idx,trainSize,replace=FALSE) #Array index of train Samples
hout_idx = temp_idx[!(temp_idx %in% train_idx )]     #Array Index of Hold Out Samples

X <-      x[train_idx,]
Y <-      y[train_idx,]
X_hout <- x[hout_idx,]
Y_hout <- y[hout_idx,]

#source('help_files.R')
#source('TIES_SAMPLING_ALG2.R')


#out = tensorProductEstimate(Y,X,X_hout,NMCMC=1500,BURNIN=200,M,NCOMIX)
#save(out,file=paste0("TIES_",NCOMIX,"_",M,".Rdata"))

#data = data.frame(y = out$function.estimate[,1],lat = hTXVEC[,1],lon = hTXVEC[,2])
#data = data.frame(y = out$m2[unique(hID)] ,lat = hTXVEC[,1],lon = hTXVEC[,2])

#data_real = data.frame(y = y1, lat = x[,1],lon = x[,2])

#ggplot(data=data,aes(x=lat,y=lon,color = y))+geom_point() +
#       scale_color_viridis_c(option="magma")
