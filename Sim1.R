##############################################################################
# R code to simulate data Y=(Y1,Y2) from the bivariate spatial model:
#
# Y1 = Z1 + W1 + e1,    e1 ~ N(0,tausq1*I)
# Y2 = Z2 + W2 + e2,    e2 ~ N(0,tausq2*I)
#
# Z1 ~ N(F*beta1,sigmasq1*R1)    (Spatial Gaussian Process 1 (R1=Matern))
# Z2 ~ N(F*beta2,sigmasq2*R2)    (Spatial Gaussian Process 2 (R2=Matern))
#
# [W1,W2] = [eta1,eta2] Lambda^T   (Factor Model for Cross Correlations)
#=============================================================================
#
# Created: Friday, 4 February, 2022 by Jonathan Stroud
# Updated: Sunday, 6 February, 2022 by Jonathan Stroud
#
# Purpose: Simulate bivariate spatial data to test FIFA code for
#          Environmetrics paper with Monica Pirani, Ed Boone, Matt Wheeler,
#          Wesley Burr, and Ben Swallow.
#
##############################################################################

# Set the working directory
# setwd("/Users/ed/Documents/GitHub/Spatial/SpatialProjectCode/")
setwd("~/Documents/GitHub/SpatialProjectCode/")
#install.packages("fields")
#install.packages("geoR")
#library(geoR)
#library(RandomFieldsUtils)
library(fields)
library(Matrix)
library(ggplot2)
library(spNNGP)
library(spBayes)
# data("MI_TSCA")
# head(MI_TSCA)
library(schnellerGP)
source( "var_pred_code.R" )
source( "help_files.R" )
source( "set_up_data_TIES.R" )

# set.seed(1234)
# Simulation parameters...
nSim1 <- 10               # Number of simulations
n.samples <-2000          # Number of MCMC samples for spNNGP and spBayes

# Lambda... the value we are varying
lambda1 <- 0.5

#n=number of observations  (n1,n2 are grid dimensions if spatial.design="grid.2d")
n=100;  n1=10; n2=10
#n=1024; n1=32; n2=32
#n=400;  n1=20; n2=20

########################################
# spatial design (grid or random)
########################################

#spatial.design = "grid.1d"
#spatial.design = "random.2d"
spatial.design = "grid.2d"

if (spatial.design=="grid.1d"){    # Grid on unit interval (1d)
  x1 = 1:n/n    
  x2 = rep(0,n)
}
if (spatial.design=="grid.2d"){    # Grid on unit interval (2d)
  lon = seq(0,1,len=n1)
  lat = seq(0,1,len=n2)
  xy = expand.grid(lon,lat)
  x1 = xy[,2]
  x2 = xy[,1]  
  # For predictive analysis
  lat_p <- (lat + min( diff( lat ) )/2 )[1:(n1-1)]
  lon_p <- (lon + min( diff( lon) )/2 )[1:(n2-1)]
  xy_p <- expand.grid( lon_p, lat_p )
  x1_p <- xy_p[,2]
  x2_p <- xy_p[,1]
}
if (spatial.design=="random.2d"){  # Random on unit square (2d)
  x1 = sort(runif(n))  
  x2 = runif(n)
}
locs = cbind(x1,x2)
locs_p = cbind( x1_p, x2_p )
n_p <- nrow( locs_p )

#########################################################################
#  True parameters 
p = 1
F = matrix(1,n,p)                             # Design Matrix 
F_f = matrix( 1, n+n_p, p )
#Lambda = matrix(0,2,2)                       # Zero Out Factor Loadings
I = diag( n )
I_p = diag( (n1-1)*(n2-1) )

# Parameters (1=Variable 1;  2=Variable 2)
beta1=rep(0,p); sigmasq1=1; phi1=0.3; kappa1=0.5; tausq1=.001
beta2=rep(0,p); sigmasq2=2; phi2=1.0; kappa2=0.5; tausq2=.001
Lambda = matrix(c(1,0,lambda1 ,1),byrow=T,ncol=2)   # Factor Loading Matrix

# Distance Matrix
D = as.matrix(dist(locs))
D_f = as.matrix( dist( rbind( locs, locs_p ) ) )

# Containers to hold the results.
HodlR.time1 <- rep( 0, nSim1 )
HodlR.mse1 <- rep( 0, nSim1 )
sgNNP.mse1 <- rep( 0, nSim1 )
sgNNP.time1 <- rep( 0, nSim1 )
spB.time1 <- rep( 0, nSim1 )
spB.mse1 <- rep( 0, nSim1 )


#########################################################################
# Simulation Run
#########################################################################
for( iter1 in 1:nSim1 ){
  #########################################################################
  #  Data Generation
  # Spatial Correlation Matrices
  #R1 = Matern(D,range=phi1,smoothness=kappa1)
  #R2 = Matern(D,range=phi2,smoothness=kappa2)
  
  # Full Spatial Correlation matrices
  R1_f = Matern(D_f, range = phi1, smoothness = kappa1 )
  R2_f = Matern(D_f, range = phi2, smoothness = kappa2 )
  
  #R1 = R1_f[ 1:n, 1:n ]
  #R2 = R2_f[ 1:n, 1:n ]
  
  # For predictive 
  #R1_p = R1_f[ (n+1):(n+n_p), (n+1):(n+n_p) ]
  #R2_p = R2_f[ (n+1):(n+n_p), (n+1):(n+n_p) ]
  
  # Cross Spatial Correlation between sites
  #R1_cp = R1_f[ 1:n, (n+1):(n+n_p)]
  #R2_cp = R2_f[ 1:n, (n+1):(n+n_p)]
  
  # Gaussian Processes
  # Z1 = F %*% beta1 + sqrt(sigmasq1) * t(chol(R1)) %*% rnorm(n)
  # Z2 = F %*% beta2 + sqrt(sigmasq2) * t(chol(R2)) %*% rnorm(n)
  
  Z1_f = F_f %*% beta1 + sqrt(sigmasq1) * t(chol(R1_f)) %*% rnorm(n+n_p)
  Z2_f = F_f %*% beta2 + sqrt(sigmasq2) * t(chol(R2_f)) %*% rnorm(n+n_p)
  
  
  # Factor Models
  eta_f = matrix(rnorm( (n+n_p)*2),n+n_p,2)
  W1_f = eta_f %*% t(Lambda)[1,]
  W2_f = eta_f %*% t(Lambda)[2,]
  
  # Observed Data
  Y1_f = Z1_f + W1_f + sqrt(tausq1)*rnorm(n+n_p)
  Y2_f = Z2_f + W2_f + sqrt(tausq2)*rnorm(n+n_p)
  
  y1 = Y1_f[ 1:n ]
  y2 = Y2_f[ 1:n ]
  
  y1_p = Y1_f[ (n+1):(n+n_p) ]
  y2_p = Y2_f[ (n+1):(n+n_p) ]
  
  HodlR.tic1 <- Sys.time()
############################################################################
#
#  Hodlr Fit the model
#
############################################################################
  x = cbind( x1, x2 )
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
  
  # For Factor analysis
  lambda = matrix(rnorm(4),2,2)
  #kronecker(Diagonal(10),solve(t(lambda)%*%lambda + diag(2)))
  
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
      lsMat[jj,1,2] = HODLR_TP_sample_ls(lsMat[jj,1,2],rY,TPscale[jj,2]*TPs[,2,jj,2],x[,1],
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
    #W1_p = eta_p %*% t(lambda)[1,]
    #W2_p = eta_p %*% t(lambda)[2,]
    
    # Data to predict
    # Y1_p = Z1_p + W1_p + sqrt(tausq1)*rnorm((n1-1)*(n2-1))
    # Y2_p = Z2_p + W2_p + sqrt(tausq2)*rnorm((n1-1)*(n2-1))
    # Data to predict
    # Y1_p = m_mean1 + W1_p + sqrt(tausq1)*rnorm((n1-1)*(n2-1))
    # Y2_p = m_mean2 + W2_p + sqrt(tausq2)*rnorm((n1-1)*(n2-1))
    m_X <- matrix( 0, nrow = length( x1 ), ncol = 1 )
    m_new <- matrix(0, nrow = length( x1_p ) , ncol = 1 )
    HodlR.pred.y <- pred_Mean_Matern( xy, xy_p, m_X, m_new, 
                                      obsv_X = m_mean1, 
                                      var = 1,  lsMat)
  }
  HodlR.time1[ iter1 ] <- Sys.time() - HodlR.tic1


tic1.sgnnp <- Sys.time()
###############################################################################
#
#  Fit standard spNNGP model
#
###############################################################################
  #n.samples <-2000
  starting <-list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=1) 
  priors <-list("phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 1), "tau.sq.IG"=c(2, 1)) 
  cov.model <-"exponential" 
  tuning <-list("phi"=0.2)
  ord <-order(xy[,1]+xy[,2]) 
  # For y1
  sim.s_y1 <-spNNGP(formula=y1~1, coords = xy , 
                 starting=starting, 
                 tuning=tuning, priors=priors, cov.model=cov.model, 
                 n.samples=n.samples, 
                 n.neighbors=10, method="latent", 
                 ord=ord, n.omp.threads=1,
                 n.report=1000, fit.rep=TRUE, 
                 sub.sample=list(start=1), 
                 return.neighbor.info = TRUE)
  sim.s_y1_pred <- predict( sim.s_y1, X.0 = matrix(1, nrow = length(Y1_p), ncol = 1 ), 
                            coords.0 = as.matrix(xy_p),
                            sub.sample=list(start=1000,thin=10),
                            n.omp.threads=1,n.report=10000) 
  sim.s_y1.mse <- sum( (Y1_p - apply( sim.s_y1_pred$p.y.0, 1, median))^2)/length(Y1_p) 
  # For y2
  sim.s_y2 <-spNNGP(formula=y2~1, coords = xy , 
                    starting=starting, 
                    tuning=tuning, priors=priors, cov.model=cov.model, 
                    n.samples=n.samples, 
                    n.neighbors=10, method="latent", 
                    ord=ord, n.omp.threads=1,
                    n.report=1000, fit.rep=TRUE, 
                    sub.sample=list(start=1), 
                    return.neighbor.info = TRUE)
  sim.s_y2_pred <- predict( sim.s_y2, X.0 = matrix(1, nrow = length(Y2_p), ncol = 1 ), 
                            coords.0 = as.matrix(xy_p),
                            sub.sample=list(start=1000,thin=10),
                            n.omp.threads=2,n.report=10000) 
  sim.s_y2.mse <- sum( (Y2_p - apply( sim.s_y2_pred$p.y.0, 1, median))^2)/length(Y1_p) 
  sgNNP.mse1[ iter1 ] <- (sim.s_y1.mse + sim.s_y2.mse)/2 
  sgNNP.time1[ iter1 ] <- Sys.time() - tic1.sgnnp


#################################################################################################
#
# spBayes
#
#################################################################################################
  
  tic.spB <- Sys.time()
  # Fit the model
  q <- 2
  # y.1 <- y[seq(1,length(y),q)]
  # y.2 <- y[seq(2,length(y),q)]
  
  A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
  #n.samples <- 1000
  burn.in <- 1
  thin <- 1
  
  starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
  tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
  priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                 "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
  
  spB.m.1 <- spMvLM(list(y1 ~ 1, y2 ~ 1), 
                coords = as.matrix(xy), starting=starting, tuning=tuning, priors=priors,
                n.samples=n.samples, cov.model="exponential", n.report=100)
  spB.m.1.pred <- spPredict(spB.m.1, start = burn.in, thin = thin, joint.TRUE,
                         pred.covars = matrix(1, nrow = length(Y2_p), ncol = 2), 
                         pred.coords = as.matrix(xy_p), verbose = FALSE)
  spB.y.hat <- apply(spB.m.1.pred$p.y.predictive.samples, 1, median )
  sqB.y1.hat <- spB.y.hat[ 1:81 ]
  sqB.y2.hat <- spB.y.hat[ 82:162 ]
  sim.sB_y1.mse <- sum( (Y1_p - sqB.y1.hat)^2)/length(Y1_p) 
  sim.sB_y2.mse <- sum( (Y2_p - sqB.y2.hat)^2)/length(Y2_p) 
  spB.mse1[ iter1 ] <- ( sim.sB_y1.mse + sim.sB_y2.mse )/2
  spB.time1[ iter1 ] <- Sys.time() - tic.spB 
  
}



sgNNP.mse1
sgNNP.time1
HodlR.time1
