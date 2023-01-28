library(ks)
library(adaFilter)
library(kernelboot)
library(qvalue)
library(data.table)
library(ggh4x)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)

# Load Methods
Rcpp::sourceCpp("Methods/JointMirror.cpp")
Rcpp::sourceCpp("Methods/modules.cpp")
source("Methods/loadJointMirror.R")
Rcpp::loadModule("JOINTMIRROR", TRUE)
source("Methods/JointMirrorR.R")
source("Methods/conditionPC.R")
source("RealData/RealData_Tools.R")
source("RealData/CDAnalysisFun.R")


## Organize Data
## Download from https://www.ibdgenetics.org
CD.tab <- read.table("RealData/cd-meta.txt",header=T,sep="\t")

# remove the SNP is NA
exp8 <- c( "belge","cedar2","adolescent","cedar1",
           "niddkj","german","niddknj","wtccc" )
exp8.P <- paste(exp8,".P",sep="")
ind.withNA <- which(apply(CD.tab[,exp8],1,function(x){any(is.na(x))}))

# Orgnize Data
simuK <- 8

CD.tab <- CD.tab[-ind.withNA,];ind.withNA <- integer(0)
CD.pval <- matrix(NA,dim(CD.tab)[1],simuK)
for(expri in 1:simuK){
  CD.pval[,expri] <- 2*(1-pnorm(abs(CD.tab[[exp8[expri]]])))
}
colnames(CD.pval) <- exp8.P
CD.tab <- cbind(CD.tab,CD.pval)
Pval <- CD.pval

## Run methods and get the selections
CD.Ana.Res <- CDAnalysis(Pval,init.thred=0.3,trgt.fdr.level=0.05)
Sel.adaFilter <- CD.Ana.Res$Sel.adaFilter
Sel.cPC <- CD.Ana.Res$Sel.cPC
Sel.EmptyPoset <- CD.Ana.Res$Sel.EmptyPoset
Sel.Max <- CD.Ana.Res$Sel.Max
Sel.Product <- CD.Ana.Res$Sel.Product
## Number of Selection
c(length(Sel.adaFilter),length(Sel.cPC),length(Sel.EmptyPoset),
  length(Sel.Max),length(Sel.Product))

## Visualize the discoveries
# Venn plot
library(gplots)

## Venn plot
venn(data = list(JM.Max=Sel.Max,
                 JM.Product=Sel.Product,
                 JM.EmptyPoset=Sel.EmptyPoset,
                 adaFilter = Sel.adaFilter,
                 cPC = Sel.cPC))

## Upset plot
CD.upset <- fromList(list(JM.Max=Sel.Max,
              JM.Product=Sel.Product,
              JM.EmptyPoset=Sel.EmptyPoset,
              adaFilter = Sel.adaFilter,
              cPC = Sel.cPC))

pdf(file="ResRealData/CDAnalysisUpSet.pdf", 
    width=8,height=5, onefile=FALSE) # or other device
upset(CD.upset,intersect=colnames(CD.upset),
      sort_intersections_by='degree',
      set_sizes=(
  upset_set_size()
  + geom_text(aes(label=..count..), hjust=1.1, stat='count')
  + expand_limits(y=140)
  # you can also add annotations on top of bars:
))
dev.off()

setEPS()
postscript("ResRealData/CDAnalysisUpSet.eps", width=8,height=5, onefile=FALSE)

upset(CD.upset,intersect=colnames(CD.upset),
      sort_intersections_by='degree',
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')
        + expand_limits(y=140)
        # you can also add annotations on top of bars:
      ))
dev.off()


### Cumulative p-value study: Leave one study out as validation
Res.CDValidate <- CDAnalysisValidVisual(Pval,init.thred=0.3)
ggplot(Res.CDValidate,aes(x=`CumulatePval(Max)`,y=`CumulatePval(Product)`))+
  # If the validation set has more than one experiment, we can add color to distinguish them
  #geom_point(aes(color=Experiment),size=0.1)+ 
  geom_point(size=0.1)+
  geom_line(aes(x=`CumulatePval(Product)`,y=`CumulatePval(Product)`),color="#FB9A99")+
  facet_wrap(vars(Val.Exp),ncol=4,scale="free")+
  theme_bw()+
  xlab(TeX("Cumulative Sum of P-values Ordered by $p_{\\max}$"))+
  ylab("Cumulative Sum of P-values Ordered by JM.Product")+
  themeRealData()

ggsave("ResRealData/CDAnalysisValidation.pdf",width=15,height=8)

ggsave("ResRealData/CDAnalysisValidation.eps",width=15,height=8)


## Use adaFilter to do partial conjuction hypothesis testing.
## AdaFilterPCH <- CDAnalysisAdaFilter(Pval,trgt.fdr.level=0.05)
## JM.Product and JM.EmptyPoset has 124 discoveries, between 127 and 50 
lapply(AdaFilterPCH,length) %>% unlist()
# [1] 953154    578    366    197    137     50     37     15


### Some Additional Vsialization
# Pval.Sel.df <- data.table(cbind(#FDP.est = Qval.Res$FDP.est.Total,
#   SNP = CD.tab$SNP,
#   CHR = CD.tab$CHR,
#   POS = CD.tab$POS,
#   Pval))
# Pval.Sel.df$Sel <- "Null"
# Pval.Sel.df$Sel[Sel.cPC] <- "cPC+adaFilter+JM.Max+JM.EmptyPoset+JM.Product"
# Pval.Sel.df$Sel[setdiff(Sel.adaFilter,Sel.cPC)] <- "adaFilter+JM.Max+JM.EmptyPoset+JM.Product"
# Pval.Sel.df$Sel[setdiff(Intersect.list(list(Sel.Max,Sel.Product,Sel.EmptyPoset)),
#                         Sel.adaFilter)] <- "JM.Max+JM.EmptyPoset+JM.Product"
# Pval.Sel.df$Sel[setdiff(intersect(Sel.Max,Sel.Product),
#                         Sel.EmptyPoset)] <- "JM.Max+JM.Product"
# Pval.Sel.df$Sel[setdiff(intersect(Sel.Product,Sel.EmptyPoset),
#                         Sel.Max)] <- "JM.EmptyPoset+JM.Product"
# Pval.Sel.df$Sel[setdiff(c(Sel.Max),
#                         c(Sel.EmptyPoset,Sel.Product))] <- "JM.Max"
# Pval.Sel.df$Sel[setdiff(c(Sel.Product),
#                         c(Sel.Max,Sel.EmptyPoset))] <- "JM.Product"
# Pval.Sel.df$Sel[setdiff(c(Sel.EmptyPoset),
#                         c(Sel.Max,Sel.Product))] <- "JM.EmptyPoset"
# Pval.Sel.df$Sel <- factor(Pval.Sel.df$Sel,levels=c("cPC+adaFilter+JM.Max+JM.EmptyPoset+JM.Product",
#                                                    "adaFilter+JM.Max+JM.EmptyPoset+JM.Product",
#                                                    "JM.Max+JM.EmptyPoset+JM.Product",
#                                                    "JM.Max+JM.EmptyPoset",
#                                                    "JM.EmptyPoset+JM.Product",
#                                                    "JM.Max+JM.Product",
#                                                    "JM.Max",
#                                                    "JM.EmptyPoset",
#                                                    "JM.Product",
#                                                    "Null"))
# 
# selected <- which(apply(Pval,1,max)<init.thred)
# selected <- c(Sel.Max,Sel.Product,Sel.EmptyPoset)
# Pval.Sel.df <- Pval.Sel.df[selected,]
# Pval.Sel.df <- melt(Pval.Sel.df,id=c("SNP","CHR","POS","Sel"),
#                     variable.name="Experiment",value.name="P-value")
# Pval.Sel.df <- Pval.Sel.df[order(-Pval.Sel.df$Sel),]
# ggplot(Pval.Sel.df,aes(x=Experiment,y=-log(`P-value`,10),group=SNP,
#                        color=Sel))+
#   scale_colour_brewer(palette = "Set1")+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#   facet_grid(Sel~.)+
#   theme(axis.text.x = element_blank(),
#         legend.position="bottom")
# #Pval.Sel.df$position = paste(letters[Pval.Sel.df$CHR],Pval.Sel.df$POS)
# p <- ggplot(Pval.Sel.df,aes(x=POS,y=-log(`P-value`,10),
#                             color=Sel))+
#   geom_point()+
#   facet_grid(Experiment~CHR,scales="free")+
#   scale_colour_brewer(palette = "Set1")+
#   force_panelsizes(cols = c(1, 1,1,1,0.3,0.3,1,1,0.3,0.3,1)) +
#   theme_light()+
#   theme_bw()+
#   theme_minimal()+
#   theme(axis.text.x = element_blank(),
#         legend.position="bottom")
# 
# ggsave(p,"ResRealData/MaPlot.jpg",height=10,width=18)
# #boxplot(Pval[apply(Pval,1,max)<=0.3,])
# 
# Pval.df <- melt(data.table(cbind(SNP = CD.tab$SNP,
#                                  Pval)),id=c("SNP"),
#                 variable.name="Experiment",value.name="P-value")
# Pval.df <- Pval.df[Pval.df$`P-value`>0.5,]
# ggplot(Pval.df,aes(x=`P-value`,color=Experiment))+
#   geom_density(bounds = c(0.5, 1))