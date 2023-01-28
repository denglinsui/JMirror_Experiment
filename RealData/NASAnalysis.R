#######################
## NAS data analysis 
#####################
library(qqman)
library(locfdr)
library(qvalue)
library(HDMT)
library(ks)
library(cowplot)
library(kernelboot)
library(ggplot2)
library(dplyr)
library(data.table)
library(latex2exp)

# Load Methods
Rcpp::sourceCpp("Methods/JointMirror.cpp")
Rcpp::sourceCpp("Methods/modules.cpp")
source("Methods/loadJointMirror.R")
source("Methods/JointMirrorR.R")
source("Methods/JointMirrorR.R")
source("RealData/RealData_Tools.R")

#########
## We use partial results of the supplementary from Liu2022 
## (Avaiable at https://www.tandfonline.com/doi//10.1080/01621459.2021.1914634)
# source("RealData/EstNull.func.R")
# source("RealData/epsest.func.R")
# source("RealData/myManhattan.R")
# source("RealData/myQQNAS.R")
load("ResRealData/DACTLiu2022.RData") ## Results after running QQforNASdata.R

##################################
## Load Summary Data
med = read.csv("RealData/mval_mediator_sorted.csv",header=T)
outcome = read.csv("RealData/mval_outcome_sorted.csv",header=T)
dat = merge(med,outcome, by="probe") ## crucial step
dat$p_DACT <- DACT_correct
n_cpg = dim(dat)[1]
Z_a = dat$z.x
Z_b = dat$z.y
p_a = dat$Pval.x
p_b = dat$Pval.y
dat$nie = dat$coef.x * dat$coef.y
dat$nie_se = sqrt(dat$coef.x^2*dat$se_coef.y^2 + dat$coef.y^2*dat$se_coef.x^2 + dat$se_coef.x^2*dat$se_coef.y^2)

dat$P_sobel = Sobel(Z_a,Z_b)

p.mat = cbind(p_a,p_b)

Pval = p.mat

## Conduct JM procedure
#-- To save time, we use init.thred = 0.1 to avoid the initialized rejection region being too large
#-- We use trgt.fdr.level = 0 to obtain the full order of unmasking
## Actually, in this study, we can just run JointMirror.R.
init.thred <- 0.1
Res.Product <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                             offset=1,trgt.fdr.level=0,
                             rank.Mode = "Product")

Res.EmptyPoset <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                             offset=1,trgt.fdr.level=0,
                             rank.Mode = "NoShape")

Res.Max <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                          offset=1,trgt.fdr.level=0,
                          rank.Mode = "Pmax")
## Run JS_mix_method
nullprop <- null_estimation(Pval)#,lambda=0.5)
pi00 <- nullprop$alpha00;pi10 <- nullprop$alpha10;pi01 <- nullprop$alpha01;
pi1 <- nullprop$alpha1;pi2 <- nullprop$alpha2
c(pi00,pi10,pi01,pi1,pi2)
# [1] 1.0000000 0.0000000 0.0000000 0.9538126 1.0000000
# pi00+pi01!=pi1 So the result is untrustable, so we do not compare with it
# pi10=pi01=0, so finite sample adjustment is not necessary.
# JS_mix_asy.Res <- fdr_est(pi00,pi01,pi10,pi1,pi2,Pval,
#                           exact = 0)
# 
# JS_mix_finite.Res <- fdr_est(alpha00=pi00,alpha01=pi01,alpha10=pi10,
#                              alpha1=pi1,alpha2=pi2,input_pvalues=Pval,
#                              exact = 1)
# MT_Comp faces severe FDR inflation, so we do not compare with it, either
# MT_Comp.Res <- qvalue(pmin(1,MT_Comp(dat$z.x, dat$z.y)),pi0=1)$qvalues

## Obtain `Q-value` to save time for recalculating
Qval.Product <- JointMirror.Qvalue(Pval,Res.Product$JMirrorS)
Qval.EmptyPoset <- JointMirror.Qvalue(Pval,Res.EmptyPoset$JMirrorS)
Qval.Max <- JointMirror.Qvalue(Pval,Res.Max$JMirrorS)

## Obtain the selected indexes
## Actually, in this study, we can just run JointMirror.R.
trgt.fdr.level <- 0.2
Sel.Product <- which(Qval.Product$is.RejSide==1 & Qval.Product$FDP.est.Total<=trgt.fdr.level)
Sel.EmptyPoset <- which(Qval.EmptyPoset$is.RejSide==1 & Qval.EmptyPoset$FDP.est.Total<=trgt.fdr.level)
Sel.Max <- which(Qval.Max$is.RejSide==1 & Qval.Max$FDP.est.Total<=trgt.fdr.level)

## Obtain the result for HDMT

# There is little difference in nullprop between the data with ties and removing ties
# Sel.JS_mix_asy <- which(JS_mix_asy.Res<= trgt.fdr.level)
# Sel.JS_mix_finite <- which(JS_mix_finite.Res<= trgt.fdr.level)
# Sel.MT_Comp <- which(MT_Comp.Res<=trgt.fdr.level)

## Obtain the result for Pmax,Sobel,DACT 
## This step is exactly the same as Liu2021 but change trgt.fdr.level to 0.2
## If one take trgt.fdr.level=0.05, they can obtain the same result in Liu2021
qobj = qvalue(DACT_correct)
Sel.DACT <- which(qobj$qvalues < trgt.fdr.level)

## Examine the results
library(gplots)
# Venn plot
venn(data = list(JM.Max=Sel.Max,
                 JM.Product=Sel.Product,
                 JM.EmptyPoset=Sel.EmptyPoset,
                 DACT = Sel.DACT))

c(length(Sel.DACT),length(Sel.Max),length(Sel.Product),length(Sel.EmptyPoset))
#  21 27 32 37


### Visualization
## To get the indexes of roots: to plot the rejection region
ProjPval <- pmin(Pval,1-Pval)
Res.Product.Ribbon <- JointMirror.R(Pval=Pval,init.thred=init.thred,
                                    offset=1,trgt.fdr.level=0.2,
                                    rank.Mode = "Product")
Ribbon.cand <- ProjPval[Res.Product.Ribbon$JMirrorS$getRootInd()+1,]

Ribbon.cand <- rbind(Ribbon.cand,c(0,max(Ribbon.cand[,2])))
Ribbon.cand <- cbind(Ribbon.cand,p_c=0)
pmax.vec <- apply(Pval,1,max)
Pmax.thred <- mean(c(max(pmax.vec[Sel.Max]),min(pmax.vec[-Sel.Max])))
Ribbon.cand2 <- cbind(p_a=c(0,Pmax.thred),
                      p_b =rep(Pmax.thred,2),p_c=rep(0,2))
Ribbon.cand.full <- data.frame(rbind(Ribbon.cand,Ribbon.cand2)) %>%
  mutate(`Rejection Region` = c(rep("JM.Product",dim(Ribbon.cand)[1]),
                                rep("JM.Max",dim(Ribbon.cand2)[1])))


## Annotate
Pval.Method.df <- data.table(
  probeID = dat$probe,
  Pval)
Pval.Method.df$Method <- "Null"
Pval.Method.df$Method[Intersect.list(list(Sel.DACT,Sel.EmptyPoset,Sel.Product,Sel.Max))] <- "DACT+JM.Max+JM.EmptyPoset+JM.Product"
Pval.Method.df$Method[setdiff(Intersect.list(list(Sel.DACT,Sel.EmptyPoset,Sel.Product)),
                              Sel.Max)] <- "DACT+JM.EmptyPoset+JM.Product"
Pval.Method.df$Method[setdiff(Intersect.list(list(Sel.DACT,Sel.Max,Sel.Product)),
                              Sel.EmptyPoset)] <- "DACT+JM.Max+JM.Product"
Pval.Method.df$Method[setdiff(intersect(Sel.EmptyPoset,Sel.Product),
                              c(Sel.DACT,Sel.Max))] <- "JM.EmptyPoset+JM.Product"
Pval.Method.df$Method[setdiff(c(Sel.EmptyPoset),
                              c(Sel.DACT,Sel.Max,Sel.Product))] <- "JM.EmptyPoset"
Pval.Method.df$Method[setdiff(intersect(Sel.EmptyPoset,Sel.Max),
                              c(Sel.DACT,Sel.Product))] <- "JM.EmptyPoset+JM.Max"
Pval.Method.df$Method <- factor(Pval.Method.df$Method,levels=(c( 
  "Null","DACT+JM.Max+JM.EmptyPoset+JM.Product",
  "DACT+JM.EmptyPoset+JM.Product",
  "DACT+JM.Max+JM.Product",
  "JM.EmptyPoset+JM.Product",
  "JM.EmptyPoset+JM.Max",
  "JM.EmptyPoset")))
Pval.Method.df <- Pval.Method.df[order(Pval.Method.df$Method),]
Pval.Method.df.pval2dplot <- Pval.Method.df %>% filter(p_a<= 0.002& p_b<=0.1)

## To plot two dimensional p-value
## present the p-values of selected SNPs
p_sel <- ggplot(Pval.Method.df.pval2dplot,aes(x=-log10(p_a),y=-log10(p_b),
                                              color=Method,
                                              shape=Method,
                                              size=Method))+
  geom_point()+
  geom_ribbon(aes(x=-log10(p_a),
                  ymin = -log10(p_b), ymax = -log10(p_c),
                  fill = `Rejection Region`),
              inherit.aes = F,alpha  = 0.1,data=Ribbon.cand.full)+
  geom_line(aes(x=-log10(p_a),y=-log10(p_a)),inherit.aes = F)+
  annotate("text", x = 10, y = 10.5, 
           label = TeX("$p_{\\gamma}=p_{\\beta}$"), parse = TRUE,angle=35,size=5)+
  xlab(TeX("$-\\log_{10}(p_{\\gamma})$"))+
  ylab(TeX("$-\\log_{10}(p_{\\beta})$"))+
  scale_fill_manual(values=NAS.fill.Panel)+
  scale_color_manual(values=NAS.color.Panel)+
  scale_shape_manual(values=NAS.shape.Panel)+
  scale_size_manual(values=NAS.size.Panel)+
  theme_bw()+
  themeRealData()+
  theme(legend.direction = "horizontal", legend.box = "vertical",
        legend.title = element_text(size=13), #change legend title font size
        legend.text = element_text(size=13))+
  guides(color = guide_legend(nrow = 3))


## To plot volcano figure
Pval.Method.df <- Pval.Method.df %>% left_join(dat, by=c("probeID"="probe") ,keep = TRUE)

p_vol <- ggplot(Pval.Method.df,aes(x=nie,y=-log10(P_sobel),
                                   color=Method,
                                   shape=Method,
                                   size=Method))+
  geom_point()+
  ylab(TeX("$-\\log_{10} (p)$"))+
  xlab("Mediation Effect Size")+
  scale_color_manual(values=NAS.color.Panel)+
  scale_shape_manual(values=NAS.shape.Panel)+
  scale_size_manual(values=NAS.size.Panel)+
  theme_bw()+
  themeRealData()

# zoomin figure
p_vol_zoomin <-  p_vol+xlim(-0.2,-0.09)+ylim(1,2.25)+  
  theme(legend.position="none",
        axis.text = element_text(),panel.border = element_rect(colour = "black", 
                                                               fill = NA),
        axis.title.x=element_blank(),axis.title.y=element_blank())

p_vol <- p_vol + 
  annotate("rect", xmin = -0.2, xmax = -0.09, ymin = 1, ymax = 2.25,
           color="black",fill=NA,linetype=2)+
  annotation_custom(
    grob = ggplotGrob(p_vol_zoomin),
    xmin = -0.35,
    xmax = 0.15,
    ymin = 4.4,
    ymax = 8.5
  )             


## merge the figures
grobs <- get_legend(p_sel)

pgrid <- plot_grid(p_sel+theme(legend.position = "none"),
                   p_vol+theme(legend.position = "none"), labels = c('A', 'B'), 
                   label_size = 12)
p_grid <- plot_grid(pgrid, grobs, nrow = 2, rel_heights = c(1, .25))
save_plot("ResRealData/NAS_p_vol.pdf",p_grid,base_height=7,base_width=11)

save_plot("ResRealData/NAS_p_vol.eps",p_grid,base_height=7,base_width=11)

save_plot("ResRealData/NAS_p_vol.jpg",p_grid,base_height=7,base_width=11,bg = "white")

## Extract the information of the selected genes
annot_file = read.table("RealData/cpg_chr_pos.txt",header=T)
Seleted.df <- Pval.Method.df.pval2dplot %>% 
  filter(Method!="Null") %>%
  left_join(dat, by=c("probeID"="probe","p_a"="Pval.x","p_b"="Pval.y")) %>%
  left_join(annot_file, by=c("probeID"="SiteID")) %>%
  select(-c("p_a","p_b"))

write.table(Seleted.df,file="ResRealData/SelectedProbes.txt")
#data.table(dat[Qval.Product$Unmask.Order[Qval.Product$Unmask.Order %in% Sel.Product],c(1,2,5,6,9,11)],dat[order(DACT_correct)[1:32],c(1,2,5,6,9,11)])