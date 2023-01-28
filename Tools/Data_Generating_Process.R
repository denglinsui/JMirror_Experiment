## Data Generating Function for Mediation Studies
DGP2 <- function(simu_n,simu_p,
                 pi_alpha1=NULL,pi_beta1=NULL,
                 simu_pi=NULL,
                 alpha1,beta0,beta1,
                 #sigma.M=1,
                 #sigma.Y=1,
                 Sigma.MY=diag(2),
                 X_type='norm',sigma_X=5,X_prob=0.2,...){
  if(X_type=='norm'){
    X <- rnorm(simu_n,0,sigma_X)
  }
  
  if(X_type=="binom"){
    X <- rbinom(simu_n,1,X_prob)
  }
  
  if((!is.null(pi_alpha1))&(!is.null(pi_beta1))){
    alpha_vec1 <- rbinom(simu_p,1,pi_alpha1)*alpha1
    beta_vec1 <- rbinom(simu_p,1,pi_beta1)*beta1
    H <- as.numeric(alpha_vec1*beta_vec1!=0)
  }else if(!is.null(simu_pi)){
    H.mat <- rmultinom(simu_p,1,simu_pi)
    alpha_vec1 <- pmax(H.mat[2,],H.mat[4,])*alpha1
    beta_vec1 <- pmax(H.mat[3,],H.mat[4,])*beta1
    H <- H.mat[4,]
  }
  
  H0 <- 1-H
  H_sign <- sign((alpha_vec1*beta_vec1>0)+0.5)
  H_XM_sign <- sign((alpha_vec1>0)+0.5)
  H_MY_sign <- sign((beta_vec1>0)+0.5)
  
  Err_vec.MY <- mvrnorm(simu_p*simu_n,c(0,0),Sigma.MY)
  Err_vec.M <- matrix(Err_vec.MY[,1],simu_n,simu_p)
  Err_vec.Y <- matrix(Err_vec.MY[,2],simu_n,simu_p)
  
  M <- t(sapply(X,function(x){x*alpha_vec1})) + Err_vec.M
  Y <- M*matrix(beta_vec1,simu_n,simu_p,byrow=T)+ matrix(X*beta0,simu_n,simu_p) + Err_vec.Y
  #Y <- M%*%beta_vec1 + X*beta0 + Err_vec.Y
  
  return(list(X=X,M=M,Y=Y,
              H=H,H0=H0,H_sign=H_sign,
              H_XM_sign=H_XM_sign,H_MY_sign=H_MY_sign,
              alpha_vec=alpha_vec1,beta_vec=beta_vec1))
}


## Data Generating Function for Replicability Studies, modified from Wang2021.
## Add decay last.
DGP_GenPMat_K <- function (M = 1000, n = 2, r = 2, all.zero.frac = 0.8, 
                       alternative.frac = 0.1, mu = c(3, 4, 5), rho = 0, 
                       decay.last=1,
                       block.size = 10, one.sided = F) {
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 
                                                         1), n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - 
                                                           n.cases - 1), nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, 
                              prob = weights), ]
  truth.pc <- rowSums(mean.matrix) >= r
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)*(seq(2-decay.last,decay.last,length.out=n))[v == 1]
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  noise.matrix <- matrix(rnorm(n * M), ncol = n)
  if (rho != 0) {
    Sigma <- matrix(rep(rho, block.size * block.size), nrow = block.size)
    diag(Sigma) <- 1
    svd.Sigma <- svd(Sigma)
    rot <- t(t(svd.Sigma$u) * sqrt(svd.Sigma$d))
    for (i in 1:(M/block.size)) {
      noise.matrix[(i - 1) * block.size + 1:block.size, 
      ] <- rot %*% noise.matrix[(i - 1) * block.size + 
                                  1:block.size, ]
    }
  }
  zmat <- noise.matrix + mean.matrix
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  }
  else {
    pvalue.matrix <- 2 * raw.pvalues
  }
  return(list(truth.pc = truth.pc, pvalue.mat = pvalue.matrix, 
              zmat = zmat,
              mean.matrix=mean.matrix))
}
