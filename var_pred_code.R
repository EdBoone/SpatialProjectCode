#################################################
#' pred_Mean_Matern
#'
#' Find the conditional mean values
#'
#' @param X The locations of the current field. 
#' @param X_new New locations we wish to predict 
#' @param m_X   Prior mean of the observed field [Usually matrix(0,nrow=length(X)]
#' @param m_new Prior mean of the new points on the observed field 
#'              [Usually matrix(0,nrow=length(X_new))]
#' @param obsv_X Observed mean of the current GP
#' @param var 
#' @param ls    Length Scale 
#' @export
#'
#' @examples
#'
pred_Mean_Matern  <- function(X,X_new,m_X,m_new,obsv_X,var,ls){
  ###########################
  X = sort(X) #make sure it is sorted
  X12 = matrix(X_new,nrow=length(X_new),ncol=length(X))
  X12 = t(X12)- X
  X12 = X12^2
  
  dist_m <- as.matrix(dist(as.matrix(X_new),diag = T,upper = T))
  dist_m = dist_m*dist_m
  
  eX12 = exp(-sqrt(5)*sqrt(X12)/ls)*(1+sqrt(5)*sqrt(X12)/ls + 5/3*(X12)/(ls^2))
  h2m <- setup_compressedMatrixGP_Matern(as.matrix(X),var,ls,1e-14,50)
  
  tX12 = eX12
  #########################
  return(m_new+t(eX12)%*%solve_HODLR(h2m,as.matrix(obsv_X - m_X)))
  #########################
}

##############################################
#' pred_COV_Matern
#'
#' Find the conditional covaraince values
#'
#' @param X The locations of the current field. 
#' @param X_new New locations we wish to predict 
#' @param var 
#' @param ls    Length Scale 
#' @export
#'
#' @examples
#'
pred_COV_matern <- function(X,X_new,var,ls){
  X = sort(X) #make sure it is sorted
  X12 = matrix(X_new,nrow=length(X_new),ncol=length(X))
  X12 = t(X12)-X
  X12 = X12^2
  
  dist_m <- as.matrix(dist(as.matrix(X_new),diag = T,upper = T))
  dist_m = dist_m*dist_m
  
  COV <- exp(-sqrt(5)*sqrt(dist_m)/ls)*(1+sqrt(5)*sqrt(dist_m)/ls + 5*dist_m/(3*ls^2)) + diag(length(X_new))*1e-8
  eX12 = exp(-sqrt(5)*sqrt(X12)/ls)*(1+sqrt(5)*sqrt(X12)/ls + 5/3*(X12)/(ls^2))
  h2m <- setup_compressedMatrixGP_Matern(as.matrix(X),var,ls,1e-14,50)
  
  tX12 = eX12
  for( ii in 1:length(X_new)){
    tX12[,ii] = solve_HODLR(h2m,eX12[,ii,drop=F])
  }
  R = t(eX12)%*%solve(RRR)%*%eX12
  V = t(eX12)%*%tX12
  return(COV - V)
}
