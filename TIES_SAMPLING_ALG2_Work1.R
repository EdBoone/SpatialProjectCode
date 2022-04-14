#'
#' Function: tensorProductEstimate()
#' Purpose: MCMC estimation of 
#' @param Y a (n x 2) matrix of observations  
#' @param X a (n x 2) matrix of spatial locations where each row
#'                    corresponds to the observation in the matrix Y
#' @param X_hout a (n x 2) matrix of spatial locations that we need
#'                         to predict a new location. 
#' @param NMCMC the number of MCMC iterations
#' @param BURNIN the number of burnin iterations
#' @param M the number of tensor producs to estimate the latent surface.  By
#'          default M = 5
#' @param K the number of splines for the predictive process: Note K^2 is the 
#'          number of unique locations, by default K = 5
#
tensorProductEstimate<-function(Y,X,X_hout,NMCMC=1500,BURNIN=200,M=5,K=5){
  #M = no TP bases
  ###################################################
  # Determine the knot set for the predictive process
  ###################################################
  tD <- K #no more than 100 knots
  X1 <- X[,1]
  tmin <- min(X1); tmax <-max(X1)
  xknots1 <- seq(tmin,tmax,(tmax - tmin)/tD)
  X2 <- X[,2]
  tmin <- min(X1); tmax <-max(X1)
  xknots2 <- seq(tmin,tmax,(tmax - tmin)/tD)
  
  
  KNOT = as.matrix(expand.grid(xknots1,xknots2))
  # Distance between knots, locs and new locs
  DTOT <- as.matrix(dist(rbind(KNOT,X,X_hout),diag=T,upper=T))
  # Distance between the knots and the locs
  DXs  <- DTOT[(nrow(KNOT)+1):(nrow(KNOT)+nrow(X)),1:nrow(KNOT)]
  # Distance between the knots and the new locs and knots
  DNs  <- DTOT[(nrow(KNOT)+nrow(X)+1):(nrow(KNOT)+nrow(X)+nrow(X_hout)),1:nrow(KNOT)]
  # Distance between the locs
  DXX  <- DTOT[(nrow(KNOT)+1):(nrow(KNOT)+nrow(X)),(nrow(KNOT)+1):(nrow(KNOT)+nrow(X))]
  diagDXX <- diag(DTOT[(nrow(KNOT)+1):(nrow(KNOT)+nrow(X)),(nrow(KNOT)+1):(nrow(KNOT)+nrow(X))])
  # Distance between the knots
  Dss  = as.matrix(dist(KNOT,diag=T,upper=T))
  rm(DTOT)
  pp1_sig2 <- 1
  pp2_sig2 <- 1
  pp1_rho <- 2
  pp2_rho <- 2
  NTOTITER = NMCMC+BURNIN
  ################################################### 
  pb <- txtProgressBar(min = 0, max = NTOTITER, style = 3)
   
  y1 <- Y[,1]
  y2 <- Y[,2]
  
  idx_lat_mi = order(X[,1])
  idx_lon_mi = order(X[,2])
  TPs <- array(1,dim= c(nrow(X),2,M,2))
  idxs   <- matrix(NA,nrow=nrow(X),2)
  idxs[,1]    <- idx_lat_mi 
  idxs[,2]    <- idx_lon_mi
    
  lsMat <- array(1,c(M,2,2))            #Matrix of length Scales
  TPscale  <- matrix(sort(rnorm(2)),M,2) #Matrix of Variances
  MUL_GAM1 <- cumprod(rgamma(M,3,1))     # Multiplicative Gamma Process 1
  MUL_GAM2 <- cumprod(rgamma(M,3,1))     # Multiplicative Gamma Process 2
  cur_tau_1 = 0.5
  cur_tau_2 = 0.5
    
  temp_tau = 100000
  cur_tau_1 = temp_tau
  cur_tau_2 = temp_tau
  etas     = rep(0,2*length(y1))
  t_etas   = matrix(etas,nrow=2)
  lambda = matrix(rnorm(4),2,2)  
  #############
  #MGP
  global_tau <- rgamma(1,1,1)
  mpg        <- matrix(rgamma(2*M,3,1),2,M)
  tp_taus   <- t(apply(mpg,1,cumprod))*global_tau
  #############
  reff_tau <- 1

  ####MCMC RUNS 
  pb <- txtProgressBar(min = 0, max = NTOTITER, style = 3)
  for(ijk in 1:NTOTITER){
      setTxtProgressBar(pb, ijk)
        
      #######################
      # Predictive Process 1
      t_etas = matrix(etas,nrow=2)
      fact_mean = lambda%*%t_etas#t(matrix(etas,ncol=2)) #latent factor mean
      m_mean1 = as.matrix(TPscale[1,1]*apply(TPs[,,1,1],1,prod))
      for (kk in 2:(dim(TPs)[3])){
        m_mean1 = m_mean1 + as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod))
      }
      rY    <- y1 - fact_mean[1,] - m_mean1
      PP1   <- sample_GIBBS_pp(rY,Dss,DXs,diagDXX,cur_tau_1,pp1_sig2,pp1_rho)
   
      temp <- sample_pp_ls_var(rY,Dss,DXs,diagDXX,pp1_sig2,pp1_rho,cur_tau_1)
      pp1_sig2 <- temp[1]
      pp1_rho  <- temp[2]
      ###
      
      #######################
      # Predictive Process 2
      t_etas = matrix(etas,nrow=2)
      fact_mean = lambda%*%t_etas#t(matrix(etas,ncol=2)) #latent factor mean
      m_mean1 = as.matrix(TPscale[1,2]*apply(TPs[,,1,2],1,prod))
      for (kk in 2:(dim(TPs)[3])){
        m_mean1 = m_mean1 + as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))
      }
     
      rY    <- y2 - fact_mean[2,] - m_mean1 
      PP2   <- sample_GIBBS_pp(rY,Dss,DXs,diagDXX,cur_tau_2,pp2_sig2,pp2_rho)
    
      temp <- sample_pp_ls_var(rY,Dss,DXs,diagDXX,pp2_sig2,pp2_rho,cur_tau_2)
      pp2_sig2 <- temp[1]
      pp2_rho  <- temp[2]
      
      #NOW UPDATE THE TP HODLR BASED ON THE RESIDUAL
      t_etas = matrix(etas,nrow=2)
      fact_mean = lambda%*%t_etas #latent factor mean
  
      ##############################
      # Sample from the additive TP y1
      for (jj in 1:(dim(TPs)[3])){
          ids = 1:(dim(TPs)[3])
          ids = ids[-jj]
          rY = matrix(y1-fact_mean[1,]-PP1) #Residual of y, mean and splines
          for (kk in ids){
            rY = rY - as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod)) #Subtract product of other tensor products
          }
          
          #sample the two Tensor Products
          TPs[,1,jj,1] = HODLR_TP_sample(rY,TPscale[jj,1]*TPs[,2,jj,1],X[,1],cur_tau_1,
                                         lsMat[jj,1,1],idxs[,1],"matern")
          TPs[,2,jj,1] = HODLR_TP_sample(rY,TPscale[jj,1]*TPs[,1,jj,1],X[,2],cur_tau_1,
                                             lsMat[jj,2,1],idxs[,2],"matern")
              
          #sample the length scale
          lsMat[jj,1,1] = HODLR_TP_sample_ls(lsMat[jj,1,1],rY,TPscale[jj,1]*TPs[,2,jj,1],X[,1],
                                             cur_tau_1,c(0.1,10),idxs[,1],kernel="matern")
          lsMat[jj,2,1] = HODLR_TP_sample_ls(lsMat[jj,2,1],rY,TPscale[jj,1]*TPs[,1,jj,1],X[,2],
                                             cur_tau_1,c(0.1,10),idxs[,2],kernel="matern")
              
          TPscale[jj,1] = HODLR_TP_sample_scale(rY,TPs[,,jj,1],cur_tau_1,mean=0,prec=tp_taus[1,jj])
      }
        
      #sample spline from residual and length scale
      rY = matrix(y1-fact_mean[1,]-PP1)
      ids = 1:(dim(TPs)[3])
      for (kk in ids){
          rY = rY - as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod)) 
      }
        
     
        
      ##############################
      # Sample from the additive TP y2
      for (jj in 1:(dim(TPs)[3])){
          ids = 1:(dim(TPs)[3])
          ids = ids[-jj]
          rY = matrix(y2-fact_mean[2,]-PP2)
          for (kk in ids){
            rY = rY - as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))
          }
          
          TPs[,1,jj,2] = HODLR_TP_sample(rY,TPscale[jj,2]*TPs[,2,jj,2],X[,1],cur_tau_2,
                                         lsMat[jj,1,2],idxs[,1],"matern")
          TPs[,2,jj,2] = HODLR_TP_sample(rY,TPscale[jj,2]*TPs[,1,jj,2],X[,2],cur_tau_2,
                                         lsMat[jj,2,2],idxs[,2],"matern")
          
          lsMat[jj,1,2] = HODLR_TP_sample_ls(lsMat[jj,1,2],rY,TPscale[jj,2]*TPs[,2,jj,2],X[,1],
                                             cur_tau_2,c(0.1,10),idxs[,1],kernel="matern")
          lsMat[jj,2,2] = HODLR_TP_sample_ls(lsMat[jj,2,2],rY,TPscale[jj,2]*TPs[,1,jj,2],X[,2],
                                             cur_tau_2,c(0.1,10),idxs[,2],kernel="matern")
          
          TPscale[jj,2] = HODLR_TP_sample_scale(rY,TPs[,,jj,2],cur_tau_2,mean=0,prec=tp_taus[2,jj])
      }
      m_mean1 =  PP1 + fact_mean[1,]; 
      for (kk in 1:(dim(TPs)[3])){
        m_mean1 = m_mean1 + as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod))
      }
      
      m_mean2 = PP2 + fact_mean[2,]; 
      for (kk in 1:(dim(TPs)[3])){
          m_mean2 = m_mean2 + as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))#Add spline TP
      }
      ##############################################
      ########## Sample Multiplicative Gamma Process
      tp_taus = tp_taus/global_tau
      temp_b = sum(t(TPscale*TPscale)*tp_taus*0.5) + 1
      temp_a = length(TPscale)*0.5 + 3
      global_tau = rgamma(1,temp_a,temp_b)
      tp_taus = tp_taus * global_tau
      for (nn in 1:2){
        for (mm in M:1){
          temp    = tp_taus[nn,mm:M]/mpg[nn,mm]
          t_scale = t(TPscale[mm:M,nn,drop=F])
          tb      = sum(t_scale^2*temp)*0.5 + 1
          ta      = length(t_scale)*0.5 + 4
          mpg[nn,mm] = rgamma(1,ta,tb)
          tp_taus[nn,mm:M] = as.numeric(temp*mpg[nn,mm])
        }
      }
      ##############################################
      # Variance
      tb = 0.5*sum((m_mean1 - y1)^2) + 1
      ta = 0.5* length(y1) + 1
      cur_tau_1 = rgamma(1,ta,tb)

      tb = 0.5*sum((m_mean2 - y2)^2) + 1
      ta = 0.5* length(y2) + 1
      cur_tau_2 = rgamma(1,ta,tb)
      ############################################### 
      m_mean1 =  PP1; 
      for (kk in 1:(dim(TPs)[3])){
        m_mean1 = m_mean1 + as.matrix(TPscale[kk,1]*apply(TPs[,,kk,1],1,prod))
      }
      
      m_mean2 = PP2; 
      for (kk in 1:(dim(TPs)[3])){
        m_mean2 = m_mean2 + as.matrix(TPscale[kk,2]*apply(TPs[,,kk,2],1,prod))#Add spline TP
      }
      # ## Random Effect
      # temp = Diagonal(n=length(y1))*(cur_tau_1 + cur_tau_2)
      # TV = solve(temp + Diagonal(n=length(y1))*reff_tau)
      # rY <- Y
      # rY[,1] <- (rY[,1] - m_mean1)*cur_tau_1
      # rY[,2] <- (rY[,2] - m_mean2)*cur_tau_2
      # rY = rowSums(rY)
      # TM = TV%*%rY
      # etas = as.numeric(TM) + as.numeric(sqrt(diag(TV)))*rnorm(length(TM))
      # etas = c(etas,etas)
      # # Random Effect Hyper-Prior Precision
      # fact_mean <- t(matrix(etas,ncol=2))
      # tB        <- 0.5*sum((fact_mean[,1]^2))+1
      # tA        <- 0.5*length(y1) + 1
      # reff_tau  <- rgamma(1,tA,tB) 
      
      #########################################
      # latent factor model Y_i = \lambda_i * \eta + N(0,diag(\tau1,\tau2))
      # sample all of the \eta_i temp_tau = 1000
      W  = matrix(c(cur_tau_1,0,0,cur_tau_2),2,2)
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
      W = c(rep(cur_tau_1,length(y1)),rep(cur_tau_2,length(y2)))
      temp_etas = matrix(etas,ncol=2,byrow =T)
      temp_etas = t(cbind(temp_etas,temp_etas))
      t_e       = matrix(temp_etas,nrow = 2*length(y1),ncol=2,byrow=T)
      TEMP1  <- Matrix(c(1,1,0,0),nrow=2*length(y1),ncol=2,byrow=T,sparse = T)
      TEMP2  <- Matrix(c(0,0,1,1),nrow=2*length(y1),ncol=2,byrow=T,sparse = T)
      temp_X = cbind(TEMP1*t_e,TEMP2*t_e)
      tV = solve(t(temp_X)%*%(temp_X*W)+diag(4))
      tM = tV%*%t(temp_X)%*%(tY*W)
      lambdas = t(chol(as.matrix(tV)))%*%matrix(rnorm(4)) + tM
      lambda = matrix(lambdas,2,2,byrow=T)

    }
  
  cov2cor(lambda%*%t(lambda) + diag(c(1/cur_tau_1,1/cur_tau_2)))
  
  TEMP_R = rowMeans(pred.func)
  PRED.FUN = matrix(TEMP_R,ncol=length(V.HASH),nrow=lengthSim,byrow=TRUE)
    
  #ALSO NEED TO OUTPUT TP OBJECTS
  out = list(function.estimate = PRED.FUN,m1=m_mean1,m2=m_mean2) #PRED.FUN is spline estimate of OOS
  return(out)
    
}
 




