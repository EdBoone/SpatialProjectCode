matern <- function(d,sig2,rho){
  rV <- sig2*(1+sqrt(3)*d/rho)*exp(-sqrt(3)*d/rho)
  return(rV)
}


#Compute the log-determinant given a Matrix
#in the Sherman Morrison Woodberry "nice" form. 
log_det_smw <- function(A,W,U){
  d_a = log(det(A + t(U)%*%solve(W)%*%U))
  d_b = sum(log(diag(W)))
  d_c = log(det(solve(A)))
  return(d_a+d_b+d_c)
}

#Gibbs sample of Predictive Process 
sample_GIBBS_pp <- function(Y,Dss,DXs,diagDXX,cur_tau,sig2,rho)
{
  Css_1 <- matern(Dss,sig2,rho)
  CXs_1 <- matern(DXs,sig2,rho)
  Css_tI <- solve(Css_1)
  Cpp_1  <- CXs_1%*%Css_tI%*%t(CXs_1)
  offset <- (matern(diagDXX,sig2,rho) - diag(Cpp_1)) 
  Cpp_1  <-  (Cpp_1 + diag(offset))
  A_inv  <- diag(1/(cur_tau*offset+1))
  
  U      <- CXs_1
  V      <- t(U)
  TMP    <- (A_inv - A_inv%*%U%*%solve((1/cur_tau)*Css_1 + V%*%A_inv%*%U)%*%V%*%A_inv)
  TV     <- Cpp_1%*%TMP
  TM   <- TV%*%(cur_tau*Y) 
  
  # Simulate Random Variable for Predictive Process
  t1 <- U%*%t(chol(Css_tI))%*%rnorm(nrow(Css_tI)) + sqrt(offset)*rnorm(length(offset))# t1~ N(0,Cpp)
  t2 <- Cpp_1%*%rnorm(length(offset),0,sqrt(cur_tau)) + t1 # t2 ~ N(0,Cpp(tau*Cpp + I))
  PP <- as.numeric(TMP%*%t2 + TM) # PP1 ~ N(TM,Cpp*[tau*Cpp + I]^-
  return(PP)
}

#Metropolis sample of length scale and variance for
#predictive process. 
sample_pp_ls_var <- function(Y,Dss,DXs,diagDXX,c_sig2,c_rho,cur_tau_1)  {
  newSig2  <- c_sig2 + runif(1,-0.25,0.25)
  newrho   <- c_rho  + runif(1,-0.25,0.25)
  if ((newrho > 0.5 & newrho < 10) & (newSig2 > 0.1 & newSig2 <5)){
    Css_1_A  <- matern(Dss,c_sig2,c_rho)
    CXs_1_A  <- matern(DXs,c_sig2,c_rho)
    Css_tI   <- solve(Css_1_A)
    Cpp_1_A  <- CXs_1_A%*%Css_tI%*%t(CXs_1_A)
    offset_A <- (matern(diagDXX,c_sig2,c_rho) - diag(Cpp_1_A)) 
    A_inv_A  <- Diagonal(n=length(offset_A),x =1/(offset_A+1/cur_tau_1))
    U_A      <- CXs_1_A
    V_A      <- t(U_A)
    CINV_A      <- (A_inv_A - A_inv_A%*%U_A%*%solve(Css_1_A + V_A%*%A_inv_A%*%U_A)%*%V_A%*%A_inv_A)
    
    log_post_den <- -0.5*(t(Y)%*%CINV_A%*%Y)
    
    Css_1_B  <- matern(Dss,newSig2,newrho)
    CXs_1_B  <- matern(DXs,newSig2,newrho)
    Css_tI   <- solve(Css_1_B)
    Cpp_1_B  <- CXs_1_B%*%Css_tI%*%t(CXs_1_B)
    offset_B <- (matern(diagDXX,newSig2,newSig2) - diag(Cpp_1_B)) 
    A_inv_B  <- Diagonal(n=length(offset_B),x =1/(offset_B+1/cur_tau_1))
    U_B      <- CXs_1_B
    V_B      <- t(U_A)
    CINV_B   <- (A_inv_B - A_inv_B%*%U_B%*%solve(Css_1_B + V_B%*%A_inv_B%*%U_B)%*%V_B%*%A_inv_B)
    log_post_num <- -0.5*(t(Y)%*%CINV_B%*%Y)
    
    log_det_den  <- log_det_smw(Css_1_A,solve(A_inv_A),U_A)
    log_det_num  <- log_det_smw(Css_1_B,solve(A_inv_B),U_B)
    
    log_num     <- -0.5*log_det_num + log_post_num
    log_den     <- -0.5*log_det_den + log_post_den
    if( runif(1) < as.numeric(exp(log_num-log_den))){
      c_sig2 = newSig2
      c_rho  = newrho
    }
  }
  return(c(c_sig2,c_rho))
}