#=================================================================================
# To implement Joint Mirror Procedure.
# We also provide a R package for readers.
#===============================================================================

JointMirror.R <- function(Pval,init.thred=0.2,
                          offset=1,trgt.fdr.level=0.1,
                          rank.Mode = "Product",
                          h=1,#is.adapt.kernel =F,
                          Ker.base = "InProjPval",
                          Ker.BW = "Silverman", # Rule of Thumb
                          is.gc = T){
  ProjPval <- pmin(Pval,1-Pval)
  
  if(Ker.base=="Pval"){
    BinP <- Pval
  }else if(Ker.base=="ProjPval"){
    BinP <- ProjPval
  }else if(Ker.base=="InProjPval"){
    ind.Ex <- which(rowSums(Pval>=0.5)>1)
    BinP <- ProjPval[-ind.Ex,]
  }else if(Ker.base=="InPval"){
    ind.Ex <- which(rowSums(Pval>=0.5)>1)
    BinP <- Pval[-ind.Ex,]
  }
  
  if(Ker.BW == "Hpi"){
    # Wand & Jones (1995) and Chacon, J.E. & Duong, T. (2010)
    Ker.mat <- Hpi(BinP)
    Hker.inv.chol <- t(solve(chol(Ker.mat)))
    if(is.gc){
      gc()
    }
    #print(Hker.inv.chol)
  }else{
    n <- dim(BinP)[1]
    d <- dim(BinP)[2]
    if( Ker.BW == "Silverman"){
      # 
      #SqrtKer.diag <-(4/(d+2))^(1/(d+4))*n^(-1/(d+4)) * apply(BinP,2,sd)
      # library(kernelboot)
       Ker.mat <- bw.silv(Pval, na.rm = FALSE)
    }
    if(Ker.BW =="Scott"){
      #SqrtKer.diag <-n^(-1/(d+4)) * apply(BinP,2,sd)
      Ker.mat <-  bw.scott(Pval,na.rm = FALSE)
    }
    #Hker.inv.chol <- diag(1/SqrtKer.diag)
    Hker.inv.chol <- t(solve(chol(Ker.mat)))
  }
  #DistMat <- as.matrix(dist(ProjPval %*% Hker.inv.chol))
  rank_Mode <- switch(rank.Mode,"Product"=1,
                      "Pmax"=2,"Max"=2,3)
  JMirrorSolver <- new(JointMirror,Pval,Hker.inv.chol)
  cat("Run JM procedure...... \n")
  JMirrorSolver$InitPara(offset_=offset,fdr_level_=trgt.fdr.level,
                         init_p_cut_=init.thred,rank_Mode_=rank_Mode)
  JMirrorSolver$runJM()
  
  # Add 1 to make the index consistent
  selected <- JMirrorSolver$getRejInd()+1
  return(list(selected=selected,
              Hker.inv.chol = Hker.inv.chol,
              JMirrorS=JMirrorSolver))
}

JointMirror.Qvalue <- function(Pval,JMirrorS){
  ## Provide ``Qvalue"
  m <- dim(Pval)[1]
  is.RejSide <- apply(Pval,1,function(x){min(x<0.5)})
  
  FDP.est.val <- JMirrorS$getFDPest()
  # In.Ind is the indexes of unmasked hypothesis
  # In.Ind also includes the order of unmasking
  In.Ind <- JMirrorS$getunMaskInd()+1
  In.Ind <- In.Ind[In.Ind!=m+1]
  Ex.Ind <- (1:m)[!(1:m %in% In.Ind)]
  
  # To ensure FDP decreases as rejection region shrinks
  FDP.est.In <- FDP.est.val[In.Ind]
  FDP.est.In <- cummin(FDP.est.In)
  # FDP.est.In <- sapply(1:length(In.Ind), 
  #                      function(i){
  #                        min(FDP.est.val[In.Ind[1:i]])})
  
  FDP.est.Total <- rep(1,m)
  FDP.est.Total[In.Ind] <-  FDP.est.In
  FDP.est.Total <- pmin(1,FDP.est.Total)
  return(list(Unmask.Order = rev(In.Ind),
              is.RejSide=is.RejSide,
              FDP.est.Total=FDP.est.Total))
}
