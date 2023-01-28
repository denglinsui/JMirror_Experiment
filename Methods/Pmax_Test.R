
Oracle_Pmax <- function(p_XM,p_MY,H0){
    p_test <- pmax(p_XM,p_MY)
    fdr_est.tmp <- sapply(p_test,function(t){
      (sum(p_test[H0==1]<=t))/max(1,sum(p_test<=t))
    })
    
    # Take Inf
    fdr_est <- sapply(p_test,function(t){
      min(fdr_est.tmp[which(p_test>=t)])
    })
  
  return(list(fdr_est=fdr_est,
              p_test=p_test))
}
