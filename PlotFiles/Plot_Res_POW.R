library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)
library(latex2exp)
library(xtable)

source("PlotFiles/Plot_Tools.R")

#######===========================
#######========= 
#######========= Plot 
#######========= 
#######===========================

BigFold.names <- c("Res/Separate_Y_Dep_pi_POW_TwoSide")

## Extract DataSet
Res.DF <- NULL

for(BigFold.name in BigFold.names){
Fold.name <- list.files(BigFold.name,full.names = T)
for(fold.name in Fold.name){
  File.name <- list.files(fold.name,full.names = T)
  Res.df <- NULL
  for(file.name in File.name){
    load(file.name)
    AdjFDP.Mat <- NULL; FDP.Mat <- NULL;POW.Mat <- NULL
    for(i in 1:nsim){
      AdjFDP.Mat <- rbind(AdjFDP.Mat,results[[i]]$AdjFDP)
      FDP.Mat <- rbind(FDP.Mat,results[[i]]$FDP)
      POW.Mat <- rbind(POW.Mat,results[[i]]$POW)
    }
    
    final_res.mean <- rbind(AdjFDP=colMeans(AdjFDP.Mat),
                            FDP=colMeans(FDP.Mat),
                            POW=colMeans(POW.Mat))
    
    final_res.se <- rbind(AdjFDP=apply(AdjFDP.Mat,2,sd),
                          FDP=apply(FDP.Mat,2,sd),
                          POW=apply(POW.Mat,2,sd))/sqrt(simu_n)
    
    colnames(final_res.mean)[colnames(final_res.mean)=="JM.NoShape"] ="JM.EmptyPoset"
    colnames(final_res.mean)[colnames(final_res.mean)=="JM.Pmax"] ="JM.Max"
    
    colnames(final_res.mean)[colnames(final_res.mean)=="JS_mix_asy"] ="JS.Mix.Asy"
    colnames(final_res.mean)[colnames(final_res.mean)=="JS_mix_finite"] ="JS.Mix.Finite"
    colnames(final_res.mean)[colnames(final_res.mean)=="MT_Comp"] ="MT.Comp"
    
    
    colnames(final_res.se)[colnames(final_res.se)=="JM.NoShape"] ="JM.EmptyPoset"
    colnames(final_res.se)[colnames(final_res.se)=="JM.Pmax"] ="JM.Max"
    colnames(final_res.se)[colnames(final_res.se)=="JS_mix_asy"] ="JS.Mix.Asy"
    colnames(final_res.se)[colnames(final_res.se)=="JS_mix_finite"] ="JS.Mix.Finite"
    colnames(final_res.se)[colnames(final_res.se)=="MT_Comp"] ="MT.Comp"
    
    res.df <- merge(melt(final_res.mean),melt(final_res.se),by=c("Var1","Var2"))
    colnames(res.df) <- c("Criterion","Method","Mean","SE")
    res.df$tau <- tau
    res.df$trgt.fdr.level <- trgt.fdr.level
    res.df$alter_type <- Transfer_pi(simu_pi)
    res.df$alter_side <- H1_type
    Res.df <- rbind(Res.df,res.df)
  }
  Res.DF <- rbind(Res.DF,Res.df)
}
}

## Extract Table for FDP
Res.df$alter_type <- factor(Res.df$alter_type,#levels=c("Global Null","Sparse Null","Dense Null",
                            #        "Sparse Alternative","Dense Alternative"))
                            levels=c("GNull","SNull","DNull",
                                     "SAlter","DAlter"))
FDPSummary005Mean <- Res.df %>% filter(Criterion=="FDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

FDPSummary020Mean <-Res.df %>% filter(Criterion=="FDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

FDPSummary005Se <- Res.df %>% filter(Criterion=="FDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type)

FDPSummary020Se <-Res.df %>% filter(Criterion=="FDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type) 

OutTable <- SummaryTable(FDPSummary005Mean,FDPSummary020Mean,
                         FDPSummary005Se,FDPSummary020Se)
print(xtable(OutTable),include.rownames = FALSE)



## Extract Table for mFDP
Res.df$alter_type <- factor(Res.df$alter_type,#levels=c("Global Null","Sparse Null","Dense Null",
                            #        "Sparse Alternative","Dense Alternative"))
                            levels=c("GNull","SNull","DNull",
                                     "SAlter","DAlter"))
mFDPSummary005Mean <- Res.df %>% filter(Criterion=="AdjFDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

mFDPSummary020Mean <-Res.df %>% filter(Criterion=="AdjFDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

mFDPSummary005Se <- Res.df %>% filter(Criterion=="AdjFDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type)

mFDPSummary020Se <-Res.df %>% filter(Criterion=="AdjFDP") %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type) 

OutTable_mFDP <- SummaryTable(mFDPSummary005Mean,mFDPSummary020Mean,
                         mFDPSummary005Se,mFDPSummary020Se)
print(xtable(OutTable_mFDP),include.rownames = FALSE)


## Extract Table for POW
POWSummary005Mean <- Res.df %>% filter(Criterion=="POW"& alter_type %in% c("SAlter","DAlter") ) %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

POWSummary020Mean <-Res.df %>% filter(Criterion=="POW"& alter_type %in% c("SAlter","DAlter")) %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="Mean") %>% 
  arrange(alter_type)

POWSummary005Se <- Res.df %>% filter(Criterion=="POW"& alter_type %in% c("SAlter","DAlter")) %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.05) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type)

POWSummary020Se <-Res.df %>% filter(Criterion=="POW"& alter_type %in% c("SAlter","DAlter")) %>% 
  filter(Method %in% c("JM.Max","JM.EmptyPoset","JM.Product","DACT","JS.Mix.Asy","JS.Mix.Finite","MT.Comp")) %>%
  filter(trgt.fdr.level == 0.2) %>%
  dcast(alter_type~Method,value.var="SE") %>% 
  arrange(alter_type) 

OutTable_POW <- SummaryTable(POWSummary005Mean,POWSummary020Mean,
                         POWSummary005Se,POWSummary020Se)
print(xtable(OutTable_POW),include.rownames = FALSE)