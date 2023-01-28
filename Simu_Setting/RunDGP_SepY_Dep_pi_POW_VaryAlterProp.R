source("InitFiles/Init_Simu.R", echo=TRUE)

######################################################################
alpha_base <- 0.2
beta_base <- 0.3

# Varying Sparsity
pi_total_coef <- matrix(
  c( 0.88,0.12, # Sparse alternatives
    0.4,0.6 # Dense alternatives
    ), 
  byrow=T,ncol=2
)
tau <- 1.25 # Corresponding tau is 0.25 and 0.375
for(i in 1:nrow(pi_total_coef)){
  simu_total_pi <- pi_total_coef[i,1:2]
  for(alter_prop in seq(1,0,length.out=10)){
  #  for(tau in seq(1,2,length.out=3)){
    simu_sep_pi <- c((1-alter_prop)/2,(1-alter_prop)/2,alter_prop)
    simu_pi <- c(simu_total_pi[1],simu_sep_pi*simu_total_pi[2])
    for(trgt.fdr.level in c(0.05,0.2)){
      alpha1 <- alpha_base*tau
      beta1 <- beta_base*tau
      save.image("Tmp.RData")
      results <- main.DGP.SepY.parallel(seeds=seeds,nsim=nsim,simu_n=simu_n,simu_p=simu_p,
                                        #pi_alpha1,pi_beta1,
                                        simu_pi=simu_pi,
                                        alpha1=alpha1,beta0=beta0,beta1=beta1,
                                        sigma.M=sigma.M,sigma.Y=sigma.Y,Sigma.MY=Sigma.MY,
                                        X_type=X_type,sigma_X=sigma_X,X_prob=X_prob,
                                        is.oneside=is.oneside,
                                        prop_Split=prop_Split,num_Split=num_Split,
                                        trgt.fdr.level=trgt.fdr.level,
                                        DGP=DGP2,num_cores=num_cores)
      file.save <- paste0(paste0(dir.subfold),"/",
                          # DGP_name,"_",H1_type,"_n_",simu_n,"_p_",simu_p,"_a_",
                          "pi_",paste0(simu_pi,collapse = "_"),
                          "_a_",alpha1,"_b_",beta1,"_fdr_",trgt.fdr.level*100,".Rdata")
      save.image(file.save)
    }
  }
}