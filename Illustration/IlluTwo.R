set.seed(100)

## Provide an example illustrating the rejection regions.

Rcpp::sourceCpp("Methods/JointMirror.cpp")
Rcpp::sourceCpp("Methods/modules.cpp")
source("Methods/loadJointMirror.R", echo=TRUE)
source("Methods/Visualize.R", echo=TRUE)
source("Methods/JointMirrorR.R", echo=TRUE)

library(kernelboot)
library(ggplot2)
library(patchwork)
library(cowplot)
m <- 2000
mu1H10 <- 1.5
mu1H11 <- 2
mu2H01 <- 2.5
mu2H11 <- 3
pi.seq <- c(0.4,0.2,0.2,0.2)
#pi.seq <- c(0.8,0,0,0.2)
H <- rmultinom(m,1,pi.seq)
mu1 <- mu2 <- rep(0,m)

mu1[H[2,]==1] <- mu1H10
mu1[H[4,]==1] <- mu1H11
mu2[H[3,]==1] <- mu2H01
mu2[H[4,]==1] <- mu2H11

X1 <- rnorm(m,mu1)
X2 <- rnorm(m,mu2)

p1 <- 2*(pnorm(-abs(X1)))
p2 <- 2*(pnorm(-abs(X2)))

H_type <- sapply(1:m,function(i){
  if(which(H[,i]==1)==1){return("H00")}
  if(which(H[,i]==1)==2){return("H10")}
  if(which(H[,i]==1)==3){return("H01")}
  if(which(H[,i]==1)==4){return("H11")}
}
)

Res.Product <- Visual_JMirror(p1,p2,fdr.level.plot=c(0.05,0.1,0.2,0.3,0.4),
                              init.thred = 0.3,
                              rank.Mode = "Product",
                              H_type=H_type)
p.Product <- Res.Product$p+
  scale_x_continuous(expand=c(1e-3,1e-3))+
  scale_y_continuous(expand=c(1e-3,1e-3))+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p.Product 
ggsave("Figure/Example_AdaRegion.pdf", width = 7, height = 7)

Res.Pmax <- Visual_JMirror(p1,p2,fdr.level.plot=c(0.05,0.1,0.2,0.3,0.4),
                              init.thred = 0.3,
                              rank.Mode = "Pmax",
                              H_type=H_type)

#### Plot Region
region.thred = 0.5
p.Pmax <- Res.Pmax$p+
  geom_segment(aes(x=region.thred,xend=region.thred,y=region.thred,yend=1),color="#FF7F00")+
  geom_segment(aes(y=region.thred,yend=region.thred,x=region.thred,xend=1),color="#FF7F00")+
  geom_segment(aes(x=0,xend=0,y=0,yend=1),color="#FF7F00")+
  geom_segment(aes(y=0,yend=0,x=0,xend=1),color="#FF7F00")+
  geom_segment(aes(x=1,xend=1,y=0,yend=region.thred),color="#FF7F00")+
  geom_segment(aes(y=1,yend=1,x=0,xend=region.thred),color="#FF7F00")+
  scale_x_continuous(expand=c(1e-3,1e-3))+
  scale_y_continuous(expand=c(1e-3,1e-3))
p.Pmax 
ggsave("Figure/Example_PmaxRegion.pdf", width = 7, height = 7)
#p.Pmax <- Res.Pmax$p
p.Pmax 
p.Pmax+p.Product+ plot_layout(guides = 'collect')+
  plot_annotation(tag_levels="A")
ggsave("Figure/Example_TwoRegion.pdf", width = 14, height = 7)
  #plot_annotation(title = "combine", tag_levels = "A")

qplot(p1,p2,color=Res.Product$JM.Solver$getProbInRej())

## Additional Result
## Obtain Rejection and FDP
Qval.Product <- JointMirror.Qvalue(cbind(p1,p2),Res.Product$JM.Solver) 
Qval.Pmax <- JointMirror.Qvalue(cbind(p1,p2),Res.Pmax$JM.Solver) 

fdp.H0.eval <- function(selected,H_type){
  (sum(H_type[selected]=="H00")+sum(H_type[selected]=="H01")+
     sum(H_type[selected]=="H10"))/max(length(selected),1)
}

fdp.adj.eval <- function(selected,H_type){
  (2*sum(H_type[selected]=="H00")+sum(H_type[selected]=="H01")+
     sum(H_type[selected]=="H10"))/max(length(selected),1)
}

Pow.eval<- function(selected,H_type){
  sum(H_type[selected]=="H11")
}
eval.mat <- NULL  
for(fdr.t.val in c(0.05,0.1,0.2,0.3)){
  Rej.Pmax <- which(Qval.Pmax$FDP.est.Total<fdr.t.val&Qval.Pmax$is.RejSide==1)
  Rej.Product <- which(Qval.Product$FDP.est.Total<fdr.t.val&Qval.Product$is.RejSide==1)
  eval.mat <- rbind(eval.mat,c(fdr.t.val,
                               fdp.H0.eval(Rej.Pmax,H_type),
                               fdp.adj.eval(Rej.Pmax,H_type),
                               Pow.eval(Rej.Pmax,H_type),
                               fdp.H0.eval(Rej.Product,H_type),
                               fdp.adj.eval(Rej.Product,H_type),
                               Pow.eval(Rej.Product,H_type)))
}

colnames(eval.mat) <- c("fdr.level","fdp","adj.fdp","true discovery","fdp","adj.fdp","true discovery")
xtable(eval.mat)

# is.RejSide <- p1<0.5&p2<0.5
# FDP.est.val <- Res.Product$JM.Solver$getFDPest()
# In.Ind <- Res.Product$JM.Solver$getunMaskInd()+1
# In.Ind <- In.Ind[In.Ind!=m+1]
# Ex.Ind <- (1:m)[!(1:m %in% In.Ind)]
# FDP.est.In <- sapply(1:length(In.Ind), function(i){min(FDP.est.val[In.Ind[1:i]])})
# FDP.true.In <- rep(1,m)
# FDP.true.In[In.Ind] <- sapply(1:length(In.Ind), 
#                               function(i){sum(H_type[In.Ind[i:length(In.Ind)]]!="H11"&
#                                                 is.RejSide[In.Ind[i:length(In.Ind)]] )/max(sum(
#                                                                                                  is.RejSide[In.Ind[i:length(In.Ind)]] ),1)}) 
# 
# FDP.true.In <- sapply(1:length(In.Ind), function(i){min(FDP.true.In[In.Ind[1:i]])})
# 
# qplot(FDP.est.In,Res.Product$JM.Solver$getProbInRej()[In.Ind])
# 
# plot(FDP.est.In,FDP.true.In)
# lines(FDP.est.In,FDP.est.In)
# 
# 
# FDP.est.val.Pmax <- Res.Pmax$JM.Solver$getFDPest()
# In.Ind.Pmax <- Res.Pmax$JM.Solver$getunMaskInd()+1
# In.Ind.Pmax <- In.Ind.Pmax[In.Ind.Pmax!=m+1]
# Ex.Ind.Pmax <- (1:m)[!(1:m %in% In.Ind.Pmax)]
# FDP.est.In.Pmax <- sapply(1:length(In.Ind.Pmax), function(i){min(FDP.est.val.Pmax[In.Ind.Pmax[1:i]])})
# FDP.true.In.Pmax <- rep(1,m)
# FDP.true.In.Pmax[In.Ind.Pmax] <- sapply(1:length(In.Ind.Pmax), 
#                               function(i){
#                                 sum(H_type[In.Ind.Pmax[i:length(In.Ind.Pmax)]]!="H11"&
#                                       is.RejSide[In.Ind.Pmax[i:length(In.Ind.Pmax)]])/max(sum(is.RejSide[In.Ind.Pmax[i:length(In.Ind.Pmax)]]),1)}) 
# 
# FDP.true.In.Pmax <- sapply(1:length(In.Ind.Pmax), function(i){min(FDP.true.In.Pmax[In.Ind.Pmax[1:i]])})
# 
# qplot(FDP.est.In.Pmax,Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax])
# 
# plot(FDP.est.In.Pmax,FDP.true.In.Pmax)
# points(FDP.est.In,FDP.true.In,col="red")
# lines(FDP.est.In.Pmax,FDP.est.In.Pmax)
# 
# plot(Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax],
#      Res.Product$JM.Solver$getProbInRej()[In.Ind])
# lines(Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax],
#       Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax])
# 
# plot(sort(Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax])[-c(1:2000)],
#           sort(Res.Product$JM.Solver$getProbInRej()[In.Ind])[-c(1:2000)])
# lines(sort(Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax])[-c(1:2000)],
#            sort(Res.Pmax$JM.Solver$getProbInRej()[In.Ind.Pmax])[-c(1:2000)])
