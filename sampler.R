library(ggplot2)
library(spNNGP)
data("MI_TSCA")
head(MI_TSCA)
library(schnellerGP)
library(Matrix)

#ggplot(data=MI_TSCA,aes(x=lat,y=long,color = WIP))+geom_point() +
#  scale_color_viridis_c(option="magma")
# colorspace examples
#ggplot(data=MI_TSCA,aes(x=lat,y=long,color =SUP))+geom_point() +
#  scale_color_viridis_c(option="magma") 

x = cbind(MI_TSCA$lat,MI_TSCA$long)
x[,1] = (x[,1]-min(x[,1]))/max(x[,1])
x[,2] = (x[,2]-min(x[,2]))/max(x[,2])
  
idx_lat_mi = order(x[,1])
idx_lon_mi = order(x[,2])
TPs <- array(1,dim= c(nrow(x),2,2,2))
idxs   <- matrix(NA,nrow=nrow(x),2)
idxs[,1]    <- idx_lat_mi 
idxs[,2]    <- idx_lon_mi

lsMat <- array(1,c(2,2,2))            #Matrix of length Scales
TPscale <- matrix(sort(rnorm(2)),2,2) #Matrix of Variances
cur_tau_1 = 0.5
cur_tau_2 = 0.5

y1 = MI_TSCA$SUP
y2 = MI_TSCA$WIP

# For Factor analysis
lambda = matrix(rnorm(4),2,2)
kronecker(Diagonal(10),solve(t(lambda)%*%lambda + diag(2)))

temp_tau = 1000
etas = rep(0,2*length(y1))
t_etas = matrix(etas,nrow=2)
#toDo : integrate the factor analysis into the 
#       model. I.E. build the residual into the factor model

for (ii in 1:25){
  #find the factor 'mean'
  t_etas = matrix(etas,nrow=2)
  fact_mean = lambda%*%t_etas #latent factor mean
  ##############################
  # Sample from the additive TP y1
  for (jj in 1:(dim(TPs)[3])){
    ids = 1:(dim(TPs)[3])
    ids = ids[-jj]
    rY = matrix(y1-fact_mean[1,])
    for (kk in ids){
      rY = rY - as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod))
    }
    
    #sample the two Tensor Products
    TPs[,1,jj,1] = HODLR_TP_sample(rY,TPscale[jj,1]*TPs[,2,jj,1],x[,1],temp_tau,
                                 lsMat[jj,1,1],idxs[,1],"matern")
    TPs[,2,jj,1] = HODLR_TP_sample(rY,TPscale[jj,1]*TPs[,1,jj,1],x[,2],temp_tau,
                                 lsMat[jj,2,1],idxs[,2],"matern")
    
    #sample the length scale
    lsMat[jj,1,1] = HODLR_TP_sample_ls(lsMat[jj,1,1],rY,TPscale[jj,1]*TPs[,2,jj,1],x[,1],
                                     temp_tau,c(0.1,10),idxs[,1],kernel="matern")
    lsMat[jj,2,1] = HODLR_TP_sample_ls(lsMat[jj,2,1],rY,TPscale[jj,1]*TPs[,1,jj,1],x[,2],
                                     temp_tau,c(0.1,10),idxs[,2],kernel="matern")
    
    TPscale[jj,1] = HODLR_TP_sample_scale(rY,TPs[,,jj,1],temp_tau)
  }
 
  m_mean1 = rY*0; 
  for (kk in 1:(dim(TPs)[3])){
    m_mean1 = m_mean1 + as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod))
  }
  
    
  ##############################
  # Sample from the additive TP y2
  for (jj in 1:(dim(TPs)[3])){
    ids = 1:(dim(TPs)[3])
    ids = ids[-jj]
    rY = matrix(y2-fact_mean[2,])
    for (kk in ids){
      rY = rY - as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))
    }
    TPs[,1,jj,2] = HODLR_TP_sample(rY,TPscale[jj,2]*TPs[,2,jj,2],x[,1],temp_tau,
                                   lsMat[jj,1,2],idxs[,1],"matern")
    TPs[,2,jj,2] = HODLR_TP_sample(rY,TPscale[jj,2]*TPs[,1,jj,2],x[,2],temp_tau,
                                   lsMat[jj,2,2],idxs[,2],"matern")
    
    
    lsMat[jj,1,2] = HODLR_TP_sample_ls(lsMat[jj,2,2],rY,TPscale[jj,2]*TPs[,2,jj,2],x[,1],
                                     temp_tau,c(0.1,10),idxs[,1],kernel="matern")
    lsMat[jj,2,2] = HODLR_TP_sample_ls(lsMat[jj,2,2],rY,TPscale[jj,2]*TPs[,1,jj,2],x[,2],
                                     temp_tau,c(0.1,10),idxs[,2],kernel="matern")
    
    TPscale[jj,2] = HODLR_TP_sample_scale(rY,TPs[,,jj,2],temp_tau)
  }
  
  m_mean2 = rY*0; 
  for (kk in 1:(dim(TPs)[3])){
    m_mean2 = m_mean2 + as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))
  }
  
  ## Location Covariance
  #latent factor model Y_i = \lambda_i * \eta + N(0,diag(\tau1,\tau2))
  # sample all of the \eta_i temp_tau = 1000
 
  W  = matrix(c(temp_tau,0,0,temp_tau),2,2)
  tX = kronecker(Diagonal(length(y1)),lambda%*%W) 
  tV = kronecker(Diagonal(length(y1)),solve(t(lambda)%*%W%*%lambda + diag(2)))
  tChol = kronecker(Diagonal(length(y1)),t(chol(solve(t(lambda)%*%W%*%lambda + diag(2)))))
  tY = matrix(t(cbind(y1-m_mean1,y2-m_mean2)))
  tM = tV%*%t(tX)%*%tY 
  etas = tChol%*%matrix(rnorm(2*length(y1))) + tM
  #sample lambda
  # 1) Reorder the etas to build regression X matrix where
  # the lambdas are our unknown quantities 
  # lambda = (lambda_{1,1},lambda_{1,2},lambda_{2,1},lambda_{2,2})'
  temp_etas = matrix(etas,ncol=2,byrow =T)
  temp_etas = t(cbind(temp_etas,temp_etas))
  t_e       = matrix(temp_etas,nrow = 2*length(y1),ncol=2,byrow=T)
  TEMP1  <- Matrix(c(1,1,0,0),nrow=2*length(y1),ncol=2,byrow=T,sparse = T)
  TEMP2  <- Matrix(c(0,0,1,1),nrow=2*length(y1),ncol=2,byrow=T,sparse = T)
  temp_X = cbind(TEMP1*t_e,TEMP2*t_e)
  tV = solve(t(temp_X)%*%temp_X*temp_tau+diag(4))
  tM = tV%*%t(temp_X*temp_tau)%*%tY
  lambdas = t(chol(as.matrix(tV)))%*%matrix(rnorm(4)) + tM
  lambda = matrix(lambdas,2,2,byrow=T)
  ##################################################################
}

