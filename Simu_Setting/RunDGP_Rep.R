source("InitFiles/Init_Simu_Rep.R", echo=TRUE)

######################################################################
#

rho <- 0.5
mu.sam <- c( 3, 4, 5)

simuK.seq <- c(2,4,6,8)

alternative.frac.seq <- alter_frac#0.01,0.03
all.zero.frac.seq <- c(0.8,0.95)
block.size.seq <- c(100,10)
decay.last.seq <- c(1,0.5)
trgt.fdr.level.seq <- c(0.05,0.2)

config.mat <- expand.grid(simuK=simuK.seq,
                          alternative.frac= alternative.frac.seq,
                          all.zero.frac= all.zero.frac.seq,
                          block.size=block.size.seq,
                          decay.last=decay.last.seq,
                          trgt.fdr.level=trgt.fdr.level.seq)
#config.mat <- config.mat[1:2,]
for(i in 1:nrow(config.mat)){
  config.item <- config.mat[i,]
  for( x in colnames(config.item)){assign(x,config.item[[x]])}
  #alpha1 <- alpha_base*tau
  #beta1 <- beta_base*tau
  save.image("Tmp.RData")
  results <- main.DGP.Rep.K.parallel(seeds,nsim,
                                     M = simuM, n = simuK, 
                                     r = simuK, alternative.frac = alternative.frac, 
                                     all.zero.frac = all.zero.frac,rho=rho,
                                     block.size =block.size,
                                     mu = mu.sam,decay.last=decay.last,
                                     trgt.fdr.level=trgt.fdr.level,
                                     num_cores=num_cores)
  file.save <- paste0(paste0(dir.subfold),"/",
                      # DGP_name,"_",H1_type,"_n_",simu_n,"_p_",simu_p,"_a_",
                      paste0(config.item,collapse = "_"),".Rdata")
  save.image(file.save)
  
}