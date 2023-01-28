library(reshape2)
library(ggplot2)
library(dplyr)
library(data.table)
library(latex2exp)
library(xtable)
library(facetscales)
library(cowplot)

source("PlotFiles/Plot_Tools.R")
#######===========================
#######========= 
#######========= Plot 
#######========= 
#######===========================

BigFold.names <- c("Res/Separate_Y_Dep_pi_POW_VaryAlterProp_TwoSide")
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
                            POW=apply(POW.Mat,2,sd))/sqrt(simu_n)
      
      colnames(final_res.mean)[colnames(final_res.mean)=="JS_mix_asy"] ="JS.Mix.Asy"
      colnames(final_res.mean)[colnames(final_res.mean)=="JS_mix_finite"] ="JS.Mix.Finite"
      colnames(final_res.mean)[colnames(final_res.mean)=="MT_Comp"] ="MT.Comp"
      colnames(final_res.mean)[colnames(final_res.mean)=="JM.NoShape"] ="JM.EmptyPoset"
      colnames(final_res.mean)[colnames(final_res.mean)=="JM.Pmax"] ="JM.Max"
      
      colnames(final_res.se)[colnames(final_res.se)=="JS_mix_asy"] ="JS.Mix.Asy"
      colnames(final_res.se)[colnames(final_res.se)=="JS_mix_finite"] ="JS.Mix.Finite"
      colnames(final_res.se)[colnames(final_res.se)=="MT_Comp"] ="MT.Comp"
      colnames(final_res.se)[colnames(final_res.se)=="JM.NoShape"] ="JM.EmptyPoset"
      colnames(final_res.se)[colnames(final_res.se)=="JM.Pmax"] ="JM.Max"
      
      res.df <- merge(melt(final_res.mean),melt(final_res.se),by=c("Var1","Var2"))
      colnames(res.df) <- c("Criterion","Method","Mean","SE")
      res.df$tau <- tau
      res.df$trgt.fdr.level <- trgt.fdr.level
      res.df$alter_type <-simu_pi[1]
      res.df$alter_side <- H1_type
      res.df$alter_prop <- alter_prop
      Res.df <- rbind(Res.df,res.df)
    }
    Res.DF <- rbind(Res.DF,Res.df)
    if(!is.null(Res.df)){
      Res.df$alter_type <- factor(Res.df$alter_type)
      Res.df$trgt.fdr.level.val <- Res.df$trgt.fdr.level
      Res.df$trgt.fdr.level <- factor(Res.df$trgt.fdr.level)
      levels(Res.df$alter_type) <- c(`0.4` = TeX("$\\pi_{\"00\"}=0.4$"), 
                                     `0.88` = TeX("$\\pi_{\"00\"}=0.88$"))
      levels(Res.df$trgt.fdr.level) <- c(`0.05` = TeX("$q=0.05$"), 
                                     `0.2` = TeX("$q=0.2$"))
      for(criterion in c("FDP","POW","AdjFDP")){
        Res.plot.df <- Res.df %>% as.data.table() %>% 
          # filter(tau>=1.0&tau<=2.0) %>%
          
          #filter(alter_prop>0) %>%
          filter(Criterion == criterion) %>%
          filter(!(Method %in% c("DACT.JC",
                                 "DACT.Efron","Emp",#"MT.Comp",
                                 #"Sym2",
                                 #"JS_mix_asy",
                                 "Sym.pval.OneStep.0","Sym.pval.OneStep.1",
                                 "Sym.pval.0","Sym.pval.1",
                                 "Pmax","Sobel",
                                 "JM.Pmax.AdaKer",
                                 "TwoDMirror.try.0","TwoDMirror.try.1",
                                 "TwoDMirror.0","TwoDMirror.1",
                                 "TwoDMirror.finite.0","TwoDMirror.asy.0","TwoDMirror.Ker.0",
                                 "Mirror_mix_asy","Mirror_mix_finite",
                                 "Pmax.Oracle","Sobel.Oracle")))
        if(criterion%in%c("AdjFDP","FDP")){
          Res.plot.df$Mean <- pmin(Res.plot.df$Mean,Res.plot.df$trgt.fdr.level.val*2)
        }
        p <- ggplot(Res.plot.df,aes(x=alter_prop,y=Mean,
                                    group=Method,
                                    color=Method,
                                    linetype=Method,
                                    shape=Method))+
          geom_line()+
          geom_point()+
          #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.02)+
          ylab(ylab.change(criterion))+
          xlab(TeX("$\\tilde{\\pi}_{1}$"))+
          scale_color_manual(values=MedAna.color.Panel)+
          scale_shape_manual(values=MedAna.shape.Panel)+
          #scale_color_brewer(palette = "Paired")+
          #scale_fill_brewer(palette = "Paired")+
          theme_bw() +theme_Med()
        if(criterion=="POW"){
          p <- p +
            facet_grid(rows=vars(trgt.fdr.level),
                       cols=vars(alter_type),
                       labeller=label_parsed, scales="free_y")
        }
        if(criterion%in%c("AdjFDP","FDP")){
          data.inter <- expand.grid(alter_type=unique(Res.plot.df$alter_type),
                                    trgt.fdr.level=unique(Res.plot.df$trgt.fdr.level))
          map.df <- unique(Res.plot.df%>%dplyr::select(trgt.fdr.level,trgt.fdr.level.val))
          data.inter <- data.inter %>% right_join(map.df)
          #data.inter$yinter <- yinter.val[!is.na(yinter.val)]
          p <- p+
            geom_hline(aes(yintercept=trgt.fdr.level.val),data.inter,linetype=2)+
            facet_grid_sc(rows=vars(trgt.fdr.level),
                          cols=vars(alter_type),
                          scales = list(y = scales_y_MedAna),
                          labeller=label_parsed)
        }
        print(p)
        assign(paste0("p.",criterion),p)
        fig.name <- paste0("Figure/",DGP_name,"_",H1_type,
                           "_",strsplit(fold.name,"/")[[1]][3],"_",criterion)
        fig.eps.name <- paste0(fig.name,".eps")
        fig.pdf.name <- paste0(fig.name,".pdf")
          if(criterion%in%c("AdjFDP","FDP")){
      #  ggsave(fig.pdf.name,height=4,width=10)
      #  ggsave(fig.eps.name,height=4,width=10)
          }else{
      #      ggsave(fig.pdf.name,height=5,width=10)
      #      ggsave(fig.eps.name,height=5,width=10)
          }
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
                                     strip.text.x = element_blank()), labels = c('A', 'B','C'),
                         nrow=3,
                         rel_heights = c(0.9, 0.7,1),
                         label_size = 12)
      p_grid <- plot_grid(pgrid, grobs, nrow = 2, rel_heights = c(1, .07))
      
      
      fig.name <- paste0("Figure/",DGP_name,"_",H1_type,
                         "_",strsplit(fold.name,"/")[[1]][3])
      fig.eps.name <- paste0(fig.name,".eps")
      fig.pdf.name <- paste0(fig.name,".pdf")
      
      save_plot(fig.eps.name,p_grid,base_height=12,base_width=10)
      save_plot(fig.pdf.name,p_grid,base_height=12,base_width=10)
    }
  }
}
