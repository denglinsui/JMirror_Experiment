# Fit Mediation model
# Obtain p-values for simulations of mediation analysis
Fit_Med <- function(X=X,M=M,Y=Y,H_sign=NULL,H_XM_sign=NULL,
                    H_MY_sign=NULL){
  ## It seems H_sign is not proper here
  if(is.null(H_XM_sign)&is.null(H_MY_sign)){
    alternative <- "two.sided"
  }else{
    alternative <- "one.sided"
  }
  
  simu_n=dim(M)[1]
  simu_p=dim(M)[2]
  
  #p_XM <- numeric(simu_p); p_MY <- numeric(simu_p)
  #coef_XM.est <- coef_XM.sd <- numeric(simu_p)
  #coef_MY.est <- coef_MY.sd <- numeric(simu_p)
  
  coef_XM.est <- as.vector((t(X)%*%M)/(t(X) %*% X)[1])
  coef_XM.sd<- sqrt(sapply(1:simu_p,
                           function(i){
                             sum((M[,i]-coef_XM.est[i]*X)^2)
                           })*solve(crossprod(X))[1]/(simu_n-1))
  
  P <- diag(simu_n)-X%*%solve((t(X)%*%X),t(X))
  
  # When Y is a vector
  if(is.null(Y)){
    y <- Y
    coef_MY.est0 <- sapply(1:simu_p,
                           function(i){m <- M[,i];mX <- cbind(m,X)
                           solve(t(mX)%*%mX,t(mX)%*%y)
                           #(t(m)%*%P%*%y)/(t(m)%*%P%*%m)
                           })
    coef_MY.est <- coef_MY.est0[1,]
    coef_XY.est <- coef_MY.est0[2,]
    coef_MY.sd<- sqrt(sapply(1:simu_p,
                             function(i){m <- M[,i]; mX <- cbind(m,X)
                             sum((y-m*coef_MY.est[i]-X*coef_XY.est[i])^2)*(solve(crossprod(mX))[1,1])
                             })/(simu_n-1))
  }else if(all(dim(Y)==dim(M))){
    coef_MY.est0 <- sapply(1:simu_p,
                           function(i){m <- M[,i]; y <- Y[,i];mX <- cbind(m,X)
                           solve(t(mX)%*%mX,t(mX)%*%y)
                           #(t(m)%*%P%*%y)/(t(m)%*%P%*%m)
                           })
    coef_MY.est <- coef_MY.est0[1,]
    coef_XY.est <- coef_MY.est0[2,]
    coef_MY.sd<- sqrt(sapply(1:simu_p,
                             function(i){m <- M[,i]; y <- Y[,i];mX <- cbind(m,X)
                             sum((y-m*coef_MY.est[i]-X*coef_XY.est[i])^2)*(solve(crossprod(mX))[1,1])
                             })/(simu_n-1))
    
  }
  
  # Obtain statistics
  z_XM <- coef_XM.est/coef_XM.sd
  z_MY <- coef_MY.est/coef_MY.sd
  
  # Obtian p-values
  # For consistency, we use normal distribution to calculate pvalues
  if(alternative == "two.sided"){
    # p_XM <- (1-pt(abs(z_XM),n-1))*2
    # p_MY <- (1-pt(abs(z_MY),n-1))*2
    p_XM <- (1-pnorm(abs(z_XM)))*2
    p_MY <- (1-pnorm(abs(z_MY)))*2
  }else{
    # p_XM <- (1-pt(H_sign*z_XM,n-1))
    # p_MY <- (1-pt(H_sign*z_MY,n-1))
    p_XM <- (1-pnorm(H_XM_sign*z_XM))
    p_MY <- (1-pnorm(H_MY_sign*z_MY))
  }
  
  return(list(p_XM=p_XM,
              p_MY=p_MY,
              z_XM=z_XM,
              z_MY=z_MY,
              coef_XM.est=coef_XM.est,
              coef_XM.sd=coef_XM.sd,
              coef_MY.est=coef_MY.est,
              coef_MY.sd=coef_MY.sd))
}

# Fit Mediation model: we just need sign, so we remove the estimation of sd to accelarate
Fit_Med_Sign <- function(X,M,Y){
  simu_n=dim(M)[1]
  simu_p=dim(M)[2]
  
  #p_XM <- numeric(simu_p); p_MY <- numeric(simu_p)
  #coef_XM.est <- coef_XM.sd <- numeric(simu_p)
  #coef_MY.est <- coef_MY.sd <- numeric(simu_p)
  
  coef_XM.est <- as.vector((t(X)%*%M)/(t(X) %*% X)[1])
  
  # When Y is a vector
  
  P <- diag(simu_n)-X%*%solve((t(X)%*%X),t(X))
  if(is.null(Y)){
    y <- Y
    coef_MY.est <- sapply(1:simu_p,
                           function(i){m <- M[,i]; 
                           (t(m)%*%P%*%y)/(t(m)%*%P%*%m)
                           })
  }else if(all(dim(Y)==dim(M))){
    # coef_MY.est0 <- sapply(1:simu_p,
    #                        function(i){m <- M[,i]; y <- Y[,i];mX <- cbind(m,X)
    #                        solve(t(mX)%*%mX,t(mX)%*%y)
    #                        #(t(m)%*%P%*%y)/(t(m)%*%P%*%m)
    #                        })
    
    coef_MY.est <- sapply(1:simu_p,
                           function(i){m <- M[,i]; y <- Y[,i]
                           (t(m)%*%P%*%y)/(t(m)%*%P%*%m)
                           })
  }
  
  Sign.est <- sign((coef_MY.est*coef_XM.est>0)+0.5)
  Sign.XM.est <- sign((coef_XM.est>0)+0.5)
  Sign.MY.est <- sign((coef_MY.est>0)+0.5)
  return(list(coef_XM.est=coef_XM.est,
              coef_MY.est=coef_MY.est,
              Sign.XM.est=Sign.XM.est,
              Sign.MY.est=Sign.MY.est,
              Sign.est=Sign.est))
}
