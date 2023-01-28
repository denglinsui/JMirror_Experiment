## Analysis Function for Crohn’s disease
## Output the selection of each methods
CDAnalysis <- function(Pval,init.thred=0.3,trgt.fdr.level=0.05){
  ## Conduct JM procedure
  #-- To save time, we use init.thred to avoid the initialized rejection region being too large
  #-- We use trgt.fdr.level = 0 to obtain the full order of unmasking
  simuK <- dim(Pval)[2]
  cat("Conducting Joint Mirror Procedure......\n")
  Res.Product <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                               offset=1,trgt.fdr.level=0,
                               rank.Mode = "Product",
                               is.gc = T)
  
  Res.EmptyPoset <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                               offset=1,trgt.fdr.level=0,
                               rank.Mode = "EmptyPoset",
                               is.gc = T)
  Res.Max <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                            offset=1,trgt.fdr.level=0,
                            rank.Mode = "Max",
                            is.gc = T)
  
  ## Obtain `Q-value` to save time for recalculating
  Qval.Product <- JointMirror.Qvalue(Pval,Res.Product$JMirrorS)
  Qval.EmptyPoset <- JointMirror.Qvalue(Pval,Res.EmptyPoset$JMirrorS)
  Qval.Max <- JointMirror.Qvalue(Pval,Res.Max$JMirrorS)
  
  #save.image("ResRealData/CD_Med_Res.RData")
  ## Obtain the selected indexes
  Sel.Product <- which(Qval.Product$is.RejSide==1 & Qval.Product$FDP.est.Total<=trgt.fdr.level)
  Sel.EmptyPoset <- which(Qval.EmptyPoset$is.RejSide==1 & Qval.EmptyPoset$FDP.est.Total<=trgt.fdr.level)
  Sel.Max <- which(Qval.Max$is.RejSide==1 & Qval.Max$FDP.est.Total<=trgt.fdr.level)
  
  ## Obtain the result for AdaFilter
  
  cat("Conducting adaFilter Procedure......\n")
  adaFilter.Res <- adaFilter(Pval, r = simuK,
                             type.I.err = "FDR",alpha=trgt.fdr.level)
  Sel.adaFilter <- which(adaFilter.Res$decision==1)
  
  
  # Try cPC
  cat("Conducting cPC Procedure......")
  cPC.Res <- conditionalPC(Pval,gamma=simuK,s=simuK,
                           sign.level = trgt.fdr.level,
                           ErrType="FDR",
                           partial.adj = "Fisher",tau=0.1)
  Sel.cPC <- cPC.Res$final.sel
  
  return(list(Sel.adaFilter = Sel.adaFilter,
              Sel.cPC=Sel.cPC,
              Sel.Product = Sel.Product,
              Sel.EmptyPoset = Sel.EmptyPoset,
              Sel.Max = Sel.Max
  ))
}

## For leave one out study
## Output the number of candidates and the rejection order of candidates for three JM procedures
## The order of JM.Pmax is exactly the order of maximum p-values across different experiments
## To avoid inconsistency, we order it by maximum of p-values

CDAnalysisValid <- function(Pval,
                            Val.Exp.P=NULL,
                            Val.Exp="belge",init.thred=0.3){
  # Val.Exp <- "belge" "cedar2" "adolescent" "cedar1" "niddkj" "german" "niddknj" "wtccc" 
  if(is.null(Val.Exp.P)){
    Val.Exp.P <- paste0(Val.Exp,".P")
  }
  Pval.Val.Rem <- Pval[,which(!(colnames(Pval)%in%Val.Exp.P))]
  Pval.Val <- Pval[,which((colnames(Pval)%in%Val.Exp.P)),drop=F]
  #-- We use trgt.fdr.level = 0 to obtain the full order of unmasking
  cat("Obtaining the rejection order for Joint Mirror Procedure......\n")
  Res.Product <- JointMirror.R(Pval=Pval.Val.Rem,init.thred=init.thred,
                               offset=1,trgt.fdr.level=0,
                               rank.Mode = "Product",
                               is.gc = T)
  
  Res.EmptyPoset <- JointMirror.R(Pval=Pval.Val.Rem,init.thred=init.thred,
                               offset=1,trgt.fdr.level=0,
                               rank.Mode = "EmptyPoset",
                               is.gc = T)
  Res.Max <- JointMirror.R(Pval=Pval.Val.Rem,init.thred=init.thred,
                            offset=1,trgt.fdr.level=0,
                            rank.Mode = "Max",
                            is.gc = T)
  
  ## Obtain `Q-value` to save time for recalculating
  Qval.Product <- JointMirror.Qvalue(Pval.Val.Rem,Res.Product$JMirrorS)
  Qval.EmptyPoset <- JointMirror.Qvalue(Pval.Val.Rem,Res.EmptyPoset$JMirrorS)
  Qval.Max <- JointMirror.Qvalue(Pval.Val.Rem,Res.Max$JMirrorS)
  
  ## Obtain Max Statistics: also reflect the order in AdaFilter
  p.max <- apply(Pval.Val.Rem,1,max)
  
  num.candidate <- sum(p.max<=init.thred)
  
  Max.Order <- Qval.Max$Unmask.Order[Qval.Max$Unmask.Order %in% 
                                         which(Qval.Max$is.RejSide==1)][1:num.candidate]
  
  Product.Order <- Qval.Product$Unmask.Order[Qval.Product$Unmask.Order %in% 
                                               which(Qval.Product$is.RejSide==1)][1:num.candidate]
  
  EmptyPoset.Order <- Qval.EmptyPoset$Unmask.Order[Qval.EmptyPoset$Unmask.Order %in% 
                                               which(Qval.EmptyPoset$is.RejSide==1)][1:num.candidate]
  
  # The candidates with num.candidate smallest pmax
  Max.Order <- order(p.max)[1:num.candidate]   
  
  Pval.Max.Order <- Pval.Val[Max.Order,,drop=F]
  Pval.Product.Order <-  Pval.Val[Product.Order,,drop=F]
  Pval.EmptyPoset.Order <-  Pval.Val[EmptyPoset.Order,,drop=F]
  
  return(list(Pval.Product.Order=Pval.Product.Order,
              Pval.EmptyPoset.Order=Pval.EmptyPoset.Order,
              Pval.Max.Order=Pval.Max.Order,
              num.candidate=num.candidate))
}


## For leave one out study
## Integrative Function: Leave one out for each experiment, integrate the results for plotting
CDAnalysisValidVisual <- function(Pval,init.thred=0.3){
    # We will compare the cumulative p-values of Product and Max
    Res <- NULL
    for(Val.Exp.P in colnames(Pval)){
      CDAnalysisValid.Res <- 
        CDAnalysisValid(Pval, Val.Exp.P=Val.Exp.P,init.thred=init.thred)
      num.candidate  <- CDAnalysisValid.Res$num.candidate
      Product.tmp = data.table::melt(data.table(ID=1:num.candidate,
                                                apply(CDAnalysisValid.Res$Pval.Product.Order,2,cumsum)),
                                     id.vars = "ID",
                                     value.name = "CumulatePval(Product)",variable.name="Experiment")
      Max.tmp = data.table::melt(data.table(ID=1:num.candidate,
                                             apply(CDAnalysisValid.Res$Pval.Max.Order,2,cumsum)),
                                  id.vars="ID",
                                  value.name = "CumulatePval(Max)",variable.name="Experiment")
      Res.tmp <- full_join(Product.tmp,Max.tmp) %>% 
        mutate( Val.Exp=paste0(unlist(strsplit(Val.Exp.P,".P")),collapse = "+"))
      Res <- rbind(Res,Res.tmp)
    }
    #Res$Experiment = unlist(strsplit(Res$Experiment,".P"))
    levels(Res$Experiment) = unlist(strsplit(levels(Res$Experiment),".P"))
    return(Res)
  }


## We consider applying AdaFilter to different settings of  Partial Conjunction Hypothesis
CDAnalysisAdaFilter <- function(Pval,trgt.fdr.level=0.05){
  simuK <- dim(Pval)[2]
  adaRes <- list()
  cat("Conducting adaFilter Procedure......\n")
  for(simur in 1:simuK){
    
    adaFilter.Res <- adaFilter(Pval, r = simur,
                               type.I.err = "FDR",alpha=trgt.fdr.level)
    Sel.adaFilter <- which(adaFilter.Res$decision==1)
    
    adaRes[[simur]] <- Sel.adaFilter
  }
  
  
  return(adaRes)
}
