## The main function to run the simulations

main.DGP.SepY.parallel <- function(seeds,nsim,
                                   simu_n,simu_p,
                                   pi_alpha1=NULL,pi_beta1=NULL,
                                   simu_pi=NULL,
                                   alpha1,beta0,beta1,
                                   sigma.M=1,
                                   sigma.Y=1,
                                   Sigma.MY=diag(2),
                                   X_type='norm',sigma_X=5,X_prob=0.2,
                                   is.oneside = T,
                                   prop_Split = 1/3,
                                   prop_Split2 = 1/2,
                                   num_Split = 50,
                                   trgt.fdr.level = 0.2,
                                   offset =1,
                                   DGP=DGP2,
                                   #Ker.base="ProjPval",
                                   Ker.base="InProjPval",
                                   num_cores=12){
  one_step <- function(jj=1){
    seed1 <- seeds[jj]
    set.seed(seed1)
    Data <- DGP(simu_n=simu_n,simu_p=simu_p,
                pi_alpha1=pi_alpha1,pi_beta1=pi_beta1,
                simu_pi=simu_pi,
                alpha1=alpha1,beta0=beta0,beta1=beta1,
                Sigma.MY=Sigma.MY,X_type=X_type,sigma_X=sigma_X,X_prob=X_prob)
    X <- Data$X
    M <- Data$M
    Y <- Data$Y
    H0 <- Data$H0
    if(is.oneside){
      H_sign <- Data$H_sign
      H_XM_sign <- Data$H_XM_sign
      H_MY_sign <- Data$H_MY_sign
    }else{
      H_sign <- H_XM_sign <- H_MY_sign<- NULL
    }
    
    #==== Organize Result
    Method <- c("Pmax","Pmax.Oracle","Sobel","Sobel.Oracle","JS_mix_asy","JS_mix_finite",
                "DACT","DACT.JC","DACT.Efron","MT_Comp",
                "JM.Pmax",
                "JM.NoShape",
                "JM.Product"
    )
    
    #=== Get p-value
    FitRes <- Fit_Med(X,M,Y,H_XM_sign=H_XM_sign,H_MY_sign=H_MY_sign)
    z_XM <- FitRes$z_XM
    z_MY <- FitRes$z_MY
    p_XM <- FitRes$p_XM
    p_MY <- FitRes$p_MY
    
    # To remove ties
    if((sum(diff(sort(p_XM))==0)>0)|(sum(diff(sort(p_MY))==0)>0)|
       (sum(p_XM==0)>0)|(sum(p_MY==0)>0)){
      warnings("There are ties or zeros in p-values. For running, we remove the ties to implement HDMT in Dai2020.")
      #p_XM <- Break_ties(p_XM)
      #p_MY <- Break_ties(p_MY)
      # We remove ties and zeros
      input_pvalues <- cbind(Break_ties(p_XM),Break_ties(p_MY))
    }else{
      input_pvalues <- cbind(p_XM,p_MY)
    }
    
    #=== pmax
    Pmax.Res <- qvalue(pmax(p_XM,p_MY),pi0=1)$qvalues#,pi0=1)
    Pmax.sel <- which(Pmax.Res<=trgt.fdr.level)
    
    Pmax.Oracle.Res <- Oracle_Pmax(p_XM,p_MY,H0)$fdr_est
    Pmax.Oracle.sel <- which(Pmax.Oracle.Res<=trgt.fdr.level)
    
    #=== Sobel Test
    Sobel.Res <- qvalue(P_sobel(coef_XM.est=FitRes$coef_XM.est,coef_XM.sd=FitRes$coef_XM.sd,
                                coef_MY.est=FitRes$coef_MY.est,coef_MY.sd=FitRes$coef_MY.sd,H_sign=H_sign),
                        pi0=1)$qvalues
    Sobel.sel <- which(Sobel.Res<=trgt.fdr.level)
    
    
    Sobel.Oracle.Res <- Oracle_Sobel(z_XM,z_MY,H0,H_sign)$fdr_est
    Sobel.Oracle.sel <- which(Sobel.Oracle.Res<=trgt.fdr.level)
    
    #=== MT_Comp Huang2019
    # We use the code implemented in Liu2021
    MT_Comp.Res <- qvalue(pmin(1,MT_Comp(z_XM, z_MY)),pi0=1)$qvalues
    MT_Comp.sel <- which(MT_Comp.Res<=trgt.fdr.level)
    
    #=== Dai2020
    # There is little difference in nullprop between the data with ties and removing ties
    nullprop <- null_estimation(input_pvalues)#,lambda=0.5)
    pi00 <- nullprop$alpha00;pi10 <- nullprop$alpha10;pi01 <- nullprop$alpha01;
    pi1 <- nullprop$alpha1;pi2 <- nullprop$alpha2
    
    JS_mix_asy.Res <- fdr_est(pi00,pi01,pi10,pi1,pi2,input_pvalues,
                              exact = 0)
    JS_mix_asy.sel <- which(JS_mix_asy.Res<= trgt.fdr.level)
    
    JS_mix_finite.Res <- fdr_est(alpha00=pi00,alpha01=pi01,alpha10=pi10,
                                 alpha1=pi1,alpha2=pi2,input_pvalues=input_pvalues,
                                 exact = 1)
    JS_mix_finite.sel <- which(JS_mix_finite.Res<= trgt.fdr.level)
    
    #=== DACT Liu2021
    DACT.pvalue <- DACT.z_stat(p_a=p_XM,p_b=p_MY,Z_a=z_XM,Z_b=z_MY)
    DACT.Res <- qvalue(DACT.pvalue,pi0=1)$qvalues
    DACT.sel <- which(DACT.Res<=trgt.fdr.level)
    
    DACT.JC.pvalue <- DACT.z_stat(p_a=p_XM,p_b=p_MY,Z_a=z_XM,Z_b=z_MY,correction="JC")
    DACT.JC.Res <- qvalue(DACT.JC.pvalue,pi0=1)$qvalues
    DACT.JC.sel <- which(DACT.JC.Res<=trgt.fdr.level)
    
    DACT.Efron.pvalue <- DACT.z_stat(p_a=p_XM,p_b=p_MY,Z_a=z_XM,Z_b=z_MY,correction="Efron")
    DACT.Efron.Res <- qvalue(DACT.Efron.pvalue,pi0=1)$qvalues
    DACT.Efron.sel <- which(DACT.Efron.Res<=trgt.fdr.level)
    
    
    ## JM procedures
    init.thred <- 0.3
    JM.Pmax.Res <- JointMirror.R(cbind(p_XM,p_MY),init.thred=init.thred,
                                 offset=offset,trgt.fdr.level=trgt.fdr.level,
                                 rank.Mode = "Pmax")
    JM.Pmax.sel <- JM.Pmax.Res$selected
    
    JM.Product.Res <- JointMirror.R(cbind(p_XM,p_MY),init.thred=init.thred,
                                    offset=offset,trgt.fdr.level=trgt.fdr.level,
                                    rank.Mode = "Product")
    JM.Product.sel <- JM.Product.Res$selected
    
    
    
    JM.NoShape.Res <- JointMirror.R(cbind(p_XM,p_MY),init.thred=init.thred,
                                    offset=offset,trgt.fdr.level=trgt.fdr.level,
                                    rank.Mode = "Other")
    JM.NoShape.sel <- JM.NoShape.Res$selected
    
    Sel.Res <- sapply(Method,
                      function(x){
                        Sel <- get(paste0(x,".sel"));
                        c(adjfdp(Sel,cbind(Data$alpha_vec,Data$beta_vec)),fdp(Sel,H0),Pow(Sel,H0))})
    cat("Step: ",jj,"\n")
    print(Sel.Res)
    return(list(AdjFDP=Sel.Res[1,],
                FDP=Sel.Res[2,],
                POW=Sel.Res[3,]))
  }
  
  
  registerDoMC(cores = num_cores)
  
  rr <- foreach(jj = 1:nsim) %dopar% one_step(jj)
  
  return(rr)
}



main.DGP.Rep.K.parallel <- function(seeds,nsim,
                                   M = 10000, n = 3, 
                                   r = 3, alternative.frac = 0.1, 
                                   all.zero.frac = 0.8,rho=0.5,
                                   block.size =100,
                                   mu = c(3, 4, 5),
                                   decay.last=1,
                                   trgt.fdr.level=0.2,
                                   num_cores=12){
  one_step <- function(jj=1){
    seed1 <- seeds[jj]
    set.seed(seed1)
    K <- n
    #==== Organize Result
    Method <- c("BonferroniBH",
                "FisherBH",
                "adaFilter",
                "JM.Pmax",
                "JM.NoShape",
                "repfdr",
                "JM.Product",
                "cPC","ssa"
    )
    #,"Emp","Emp.Con") #In cluster, the empirical Bayes based method doesn't work
    #library(adaFilter);library(repfdr); K=3;alternative.frac = 0.1; all.zero.frac = 0.8;rho=0.5;block.size = 10;M=5000;trgt.fdr.level=0.05
    data <- DGP_GenPMat_K(M = M, n = K, 
                    r = K, alternative.frac = alternative.frac, 
                    all.zero.frac = all.zero.frac,rho=rho,
                    mu=mu,
                    decay.last = decay.last,
                    block.size =block.size)
    Pval <- data$pvalue.mat
    H0 <- data$truth.pc==0
    mumat <- data$mean.matrix
    zmat <- data$zmat
    
    adaFilter.Res <- adaFilter(Pval, r = K, type.I.err = "FDR",alpha=trgt.fdr.level)
    adaFilter.sel <- which(adaFilter.Res$decision==1)
    
    ## standard BH adjustment on Bonferroni PC p-values
    BonferroniBH.Res <- ClassicalMTP(data$pvalue.mat, r = K, alpha = trgt.fdr.level, 
                                     method = "Bonferroni")
    BonferroniBH.sel <- which(BonferroniBH.Res$decision)
    
    ## standard BH adjustment on Fisher PC p-values 
    FisherBH.Res <- ClassicalMTP(data$pvalue.mat, r = K, alpha = trgt.fdr.level, 
                                     method = "Fisher")
    FisherBH.sel <- which(FisherBH.Res$decision)
    
    ## repfdr
    ztobins_res = ztobins(zmat,n.association.status = 2) 
    
    #run repfdr
    if(all(ztobins_res$proportions[1,]==1)){
      repfdr.sel = integer(0)
    }else{
    repfdr.Res = repfdr(ztobins_res$pdf.binned.z,
                        ztobins_res$binned.z.mat,
                        non.null = 'replication',
                        control=em.control(verbose=F))

    repfdr.sel = which(repfdr.Res$mat[,2] <= trgt.fdr.level)
    }
    # simultaneous signal detection
    # too conservative
    ssa.sel <- nfsdr(abs(zmat), trgt.fdr.level)
    
    # conditional PC
    cPC.Res <- conditionalPC(Pval,gamma=K,s=K,
                             sign.level = trgt.fdr.level,
                             ErrType="FDR",
                             partial.adj = "Fisher",tau=0.1)
    cPC.sel <- cPC.Res$final.sel
    
    # Joint Mirror
    init.thred <- 0.3
    JM.Pmax.Res <- JointMirror.R(Pval,init.thred=init.thred,
                                 offset=1,trgt.fdr.level=trgt.fdr.level,
                                 rank.Mode = "Pmax")
    JM.Pmax.sel <- JM.Pmax.Res$selected
    
    JM.Product.Res <- JointMirror.R(Pval,init.thred=init.thred,
                                    offset=1,trgt.fdr.level=trgt.fdr.level,
                                    rank.Mode = "Product")
    JM.Product.sel <- JM.Product.Res$selected
    
    JM.NoShape.Res <- JointMirror.R(Pval,init.thred=init.thred,
                                    offset=1,trgt.fdr.level=trgt.fdr.level,
                                    rank.Mode = "Other")
    JM.NoShape.sel <- JM.NoShape.Res$selected
    
    Sel.Res <- sapply(Method,
                      function(x){
                        Sel <- get(paste0(x,".sel"));
                        c(adjfdp(Sel,mumat),fdp(Sel,H0),Pow(Sel,H0))})
    cat("Step: ",jj,"\n")
    print(Sel.Res)
    return(list(AdjFDP=Sel.Res[1,],
                FDP=Sel.Res[2,],
                POW=Sel.Res[3,]))
  }
  
  
  registerDoMC(cores = num_cores)
  
  rr <- foreach(jj = 1:nsim) %dopar% one_step(jj)
  
  return(rr)
}
