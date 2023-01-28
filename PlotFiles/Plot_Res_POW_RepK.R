library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)
library(latex2exp)
library(xtable)
library(cowplot)

library(ggh4x)
source("PlotFiles/Plot_Tools.R")
#######===========================
#######========= 
#######========= Plot 
#######========= 
#######===========================

BigFold.names <- c("Res/Rep_K1_TwoSide_0.01")
BigFold.names <- c("Res/Rep_K3_TwoSide_0.03")
## Draw plot
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
                            POW=apply(POW.Mat,2,sd))/sqrt(dim(AdjFDP.Mat)[1])
      
      
      colnames(final_res.mean)[colnames(final_res.mean)=="JM.NoShape"] ="JM.EmptyPoset"
      colnames(final_res.mean)[colnames(final_res.mean)=="JM.Pmax"] ="JM.Max"
      
      colnames(final_res.se)[colnames(final_res.se)=="JM.NoShape"] ="JM.EmptyPoset"
      colnames(final_res.se)[colnames(final_res.se)=="JM.Pmax"] ="JM.Max"
      
      res.df <- merge(melt(final_res.mean),melt(final_res.se),by=c("Var1","Var2"))
      colnames(res.df) <- c("Criterion","Method","Mean","SE")
      res.df <- cbind(res.df,config.item)
      res.df$trgt.fdr.level <- paste0("fdrlevel=",trgt.fdr.level)
      # res.df$tau <- tau
      # res.df$trgt.fdr.level <- trgt.fdr.level
      # res.df$alter_type <- Transfer_pi(simu_pi)
      # res.df$alter_side <- H1_type
      Res.df <- rbind(Res.df,res.df)
    }
    Res.DF <- rbind(Res.DF,Res.df)
    if(!is.null(Res.df)){
      Res.df$trgt.fdr.level <- factor(Res.df$trgt.fdr.level)
      Res.df$block.size <- factor(Res.df$block.size)
      Res.df$simuK <- factor(Res.df$simuK)
      Res.df$all.zero.frac <- factor(Res.df$all.zero.frac)
      Res.df$decay.last <- factor(Res.df$decay.last)
      
      
     
      
      levels(Res.df$block.size) <- c(`10` = TeX("$b=10$"), 
                                     `100` = TeX("$b=100$"))
      levels(Res.df$decay.last) <- c(`0.5` = TeX("$w_0=0.5$"), 
                                     ` 1` = TeX("$w_0=1$"))
      levels(Res.df$all.zero.frac) <- c(`0.8` = TeX("$\\pi^{global}_{0}=0.8$"), 
                                     `0.95` = TeX("$\\pi^{global}_{0}=0.95$"))
      
      data.inter <- unique(Res.df[,c("block.size","trgt.fdr.level",
                                     "all.zero.frac","decay.last")])
      data.inter$trgt.fdr.level.val <- data.inter$trgt.fdr.level
      
      yinter.val <- as.numeric(unlist(strsplit(as.character(data.inter$trgt.fdr.level),"=")))
      data.inter$yinter <- yinter.val[!is.na(yinter.val)]
      
      levels(Res.df$trgt.fdr.level) <- c(`fdrlevel=0.05` = TeX("$q=0.05$"), 
                                         ` fdrlevel=0.2` = TeX("$q=0.2$"))
      levels(data.inter$trgt.fdr.level) <- c(`fdrlevel=0.05` = TeX("$q=0.05$"), 
                                         ` fdrlevel=0.2` = TeX("$q=0.2$"))
      Res.plot.df.total <- NULL
      for(criterion in c("FDP","POW","AdjFDP")){
        Res.plot.df <- Res.df %>% as.data.table() %>% 
          # filter(tau>=1.0&tau<=2.0) %>%
          filter(!(Method%in%c("repfdr","FisherBH"))) %>%
          filter(Criterion == criterion) 
        
        
        Res.plot.df$Method <- factor(Res.plot.df$Method,
                                     levels = c("BonferroniBH","cPC","ssa","adaFilter",
                                                #"JM.Pmax","JM.NoShape","JM.Product"))
                                                "JM.Max","JM.EmptyPoset","JM.Product"))
        p <- ggplot(Res.plot.df,
                    aes(x=simuK,y=Mean,
                                    fill=Method))+
          geom_bar(stat = 'identity',position="dodge")+
          scale_fill_manual(values=RepAna.color.Panel)+
          geom_errorbar(aes(ymin = Mean-SE, ymax = Mean+SE),position = position_dodge(0.9), width = 0.5)+
          ylab(ylab.change(criterion))+
          xlab(TeX("$K$"))+
          #facet_nested(trgt.fdr.level+decay.last~all.zero.frac+block.size,
          #             scales="free_y")+
          facet_nested(rows=vars(trgt.fdr.level,decay.last),
                     cols=vars(all.zero.frac,block.size),
                     labeller=label_parsed, scales="free_y")+
          theme_bw() +
          theme_Rep()
        if(criterion%in%c("AdjFDP","FDP")){
          
          p <- p+geom_hline(aes(yintercept=yinter),data.inter,linetype=2)
        }
        print(p)
        
        assign(paste0("p.",criterion),p)
        fig.eps.name <- paste0("Figure/",DGP_name,"_",H1_type,
                               "_",paste0(config.item,collapse = "_"),"_",criterion,".eps")
        fig.pdf.name <- paste0("Figure/",DGP_name,"_",H1_type,
                               "_",paste0(config.item,collapse = "_"),"_",criterion,".pdf")
        
        if(criterion%in%c("AdjFDP","FDP")){
         # ggsave(fig.pdf.name,height=5,width=10)
        #  ggsave(fig.eps.name,height=5,width=10)
        }else{
         # ggsave(fig.pdf.name,height=6,width=10)
        #  ggsave(fig.eps.name,height=6,width=10)
        }
        Res.plot.df.total <- rbind(Res.plot.df.total,
                                   Res.plot.df %>% mutate(Criterion = criterion))
      }
      
      
      
      grobs <- get_legend(p.FDP)
      
      pgrid <- plot_grid(p.FDP+theme(legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.ticks.x  = element_blank())+
                           labs(x = NULL),
                         p.AdjFDP+theme(legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.ticks.x  = element_blank(),
                                     strip.text.x = element_blank())+
                           labs(x = NULL),
                         p.POW+theme(legend.position = "none",
                                     strip.text.x = element_blank()
                                     ), labels = c('A', 'B','C'),
                         nrow=3,
                         rel_heights = c(0.9,0.7, 1),
                         label_size = 12)
      p_grid <- plot_grid(pgrid, grobs, nrow = 2, rel_heights = c(1, .05))
      
      fig.eps.name <- paste0("Figure/",DGP_name,"_",H1_type,
                             "_",paste0(config.item,collapse = "_"),".eps")
      fig.pdf.name <- paste0("Figure/",DGP_name,"_",H1_type,
                             "_",paste0(config.item,collapse = "_"),".pdf")
      
      
      save_plot(fig.eps.name,p_grid,base_height=13,base_width=10)
      save_plot(fig.pdf.name,p_grid,base_height=13,base_width=10)
      
    }
  }
}
