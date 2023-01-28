#=================================================================================
# To implement conditional partial conjunction test in Dickhaus et al., 2021
#===============================================================================
conditionalPC <- function(Pval,gamma,s,
                          sign.level = 0.05,
                          ErrType="FDR",
                          partial.adj = "Fisher",
                          tau=0.1){
  
  if(partial.adj=="Fisher"){
    adj.fun <- function(x){
      test.val <- sum(-2*log(sort(x)[gamma:s]))
      1-pchisq(test.val,2*(s-gamma+1))
    }
  }
  if(partial.adj=="Stouffer"){
    adj.fun <- function(x){
      test.val <- sum(qnorm(1-x[gamma:s]))
      1-pnorm(test.val/sqrt(s-gamma+1))
    }
  }
  partial.p <- apply(Pval,1,adj.fun)
  
  cond.sel <- partial.p<tau
  cond.ind <- which(cond.sel)
  partial.cond.p <- pmin(1,partial.p[cond.sel]/tau )
  
  if(sum(cond.sel)!=0){
    if(ErrType=="FDR"){
      cond.qval <- qvalue(partial.cond.p,pi0=1)$qvalue
      cond.ind.sel <- cond.qval<=sign.level
      final.sel <- cond.ind[cond.ind.sel]
    }
  }else{
    final.sel <- integer(0)
  }
  
  return(list(partial.p = partial.p,
              cond.qval = cond.qval,
              final.sel=final.sel))
}
