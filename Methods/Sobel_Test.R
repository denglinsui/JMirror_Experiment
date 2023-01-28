#=================================================================================
# To implement Sobel test in Sobel, 1982
# We consider it but do not present it because its power is not competitive
#===============================================================================
P_sobel <- function(coef_XM.est,coef_XM.sd,
                    coef_MY.est,coef_MY.sd,H_sign=NULL){
  if(is.null(H_sign)){
    alternative <- "two.sided"
  }else{
    alternative <- "one.sided"
  }
  
  T.sobel <- coef_XM.est*coef_MY.est/
    sqrt(coef_MY.est^2* coef_XM.sd^2+
           coef_XM.est^2*coef_MY.sd^2)
  
  # Obtian p-values
  # For consistency, we use normal distribution to calculate pvalues
  if(alternative == "two.sided"){
    # p_XM <- (1-pt(abs(z_XM),n-1))*2
    # p_MY <- (1-pt(abs(z_MY),n-1))*2
    p.sobel <- (1-pnorm(abs(T.sobel)))*2
  }else{
    # p_XM <- (1-pt(H_sign*z_XM,n-1))
    # p_MY <- (1-pt(H_sign*z_MY,n-1))
    p.sobel <- (1-pnorm(H_sign*T.sobel))
  }
  
  return(p.sobel)
}

Oracle_Sobel <- function(z_XM,z_MY,H0,H_sign=NULL){
  if(!is.null(H_sign)){
    z_test <- z_XM*z_MY*H_sign
    fdr_est.tmp <- sapply(z_test,function(t){
      (sum(z_test[H0==1]>=t))/max(1,sum(z_test>=t))
    })
    
    # Take Inf
    fdr_est <- sapply(z_test,function(t){
      min(fdr_est.tmp[which(z_test<=t)])
    })
  }else{
    z_test <- z_XM*z_MY
    fdr_est.tmp <- sapply(abs(z_test),function(t){
      (sum(abs(z_test)[H0==1]>=t))/max(1,sum(abs(z_test)>=t))
    })
    
    # Take Inf
    fdr_est <- sapply(abs(z_test),function(t){
      min(fdr_est.tmp[which(abs(z_test)<=t)])
    })
  }
  
  return(list(fdr_est=fdr_est,
              z_test=z_test))
}
