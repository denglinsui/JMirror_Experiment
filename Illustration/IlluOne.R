set.seed(100)

## Provide an example for illustration
## We show what happens when rejection region is [0,0.2]\times[0,0.2]
## Not presented in paper

Rcpp::sourceCpp("~/project/HDMediate/Methods/JointMirror.cpp")
Rcpp::sourceCpp("~/project/HDMediate/Methods/modules.cpp")
source("~/project/HDMediate/Methods/loadJointMirror.R", echo=TRUE)
source("~/project/HDMediate/Methods/Visualize.R", echo=TRUE)
source("~/project/HDMediate/Methods/JointMirrorR.R", echo=TRUE)
library(ggplot2)
library(ggpattern)
m <- 1000
mu1H10 <- 1.5
mu1H11 <- 2
mu2H01 <- 2.5
mu2H11 <- 3
pi.seq <- c(0.4,0.2,0.2,0.2)
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

data.plot <- data.frame(p1=p1,
                        p2=p2,
                        X1=X1,
                        X2=X2,
                        H_type=H_type)
data.rect <- data.frame(
  xmin=c(0,0,0.8),#,0.8),
  xmax=c(0.2,0.2,1),#,1),
  ymin=c(0,0.8,0),#0.8),
  ymax=c(0.2,1,0.2),#,1),
  Region=c("R","A1","A2")
)

ggplot(data=data.plot)+
  geom_point(aes(x=p1,y=p2,
                 color=H_type,
                 shape=H_type),size=2)+
  scale_shape_manual(values=c(3,4,20,8))+
  scale_color_manual(values=c("#CAB2D6","#1F78B4","#33A02C","#E31A1C"))+
  #geom_vline(aes(xintercept=0.2))+
  #geom_hline(aes(yintercept=0.2))+
  geom_segment(aes(x=0.2,xend=0.2,y=0.2,yend=1))+
  geom_segment(aes(y=0.2,yend=0.2,x=0.2,xend=1))+
  geom_segment(aes(x=0,xend=0,y=0,yend=1))+
  geom_segment(aes(y=0,yend=0,x=0,xend=1))+
  geom_segment(aes(x=1,xend=1,y=0,yend=0.2))+
  geom_segment(aes(y=1,yend=1,x=0,xend=0.2))+
  geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 0.2,
           alpha = 0,linejoin="round")+
  scale_x_continuous(expand=c(1e-3,1e-3))+
  scale_y_continuous(expand=c(1e-3,1e-3))+
  geom_rect(data=data.rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill=Region),alpha=0.2)+
  scale_fill_manual(values=c("#A6CEE3","#B2DF8A","#FB9A99"))+
  theme_bw()+
  theme(#legend.position = c(.99, .8),
        #legend.justification = c("right", "top"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        axis.title=element_text(size=15),
        axis.text=element_text(size=12))#, legend.box="vertical")

ggsave("Figure/Example_Region.pdf", width = 7, height = 7)

p.type <- sapply(1:m,function(i){
  p1.sub <- p1[i]
  p2.sub <- p2[i]
  if(p1.sub<0.2&p2.sub<0.2){return("LeftDown")}
  if(p1.sub<0.2&p2.sub>0.8){return("LeftUp")}
  if(p1.sub>0.8&p2.sub<0.2){return("RightDown")}
  if(p1.sub>0.8&p2.sub>0.8){return("RightUp")}
  return("unmask")
}
)

## We can see how many hypothesis locating into the rejection and control regions

Tab.Htype <- table(H_type,p.type)
Tab.Htype

## We can compare the FDPs is bounded by our FDP estimator
## The bound is not tight as we just want to show the relationship between the control and rejection regions
## In practice, the rejection region would be much smaller.
FDP.H0 <- sum(Tab.Htype[1:3,"LeftDown"])/sum(Tab.Htype[1:4,"LeftDown"])
FDP.H00 <- sum(Tab.Htype[1,"LeftDown"])/sum(Tab.Htype[1:4,"LeftDown"])
FDP.H0d <- sum(Tab.Htype[1:2,"LeftDown"])/sum(Tab.Htype[1:4,"LeftDown"])
FDP.Hd0 <- sum(Tab.Htype[c(1,3),"LeftDown"])/sum(Tab.Htype[1:4,"LeftDown"])

FDP.adj.H0 <- (2*Tab.Htype[1,"LeftDown"]+sum(Tab.Htype[2:3,"LeftDown"]))/sum(Tab.Htype[1:4,"LeftDown"])
c(FDP.H0,FDP.H00,FDP.H0d,FDP.Hd0,FDP.adj.H0)


FDP.est.H0 <- sum(Tab.Htype[1:4,c("LeftUp","RightDown")])/sum(Tab.Htype[1:4,"LeftDown"])
c(FDP.est.H0)

# Res.JM <- JointMirror.R(cbind(p1,p2),init.thred=0.2,
#                           offset=1,trgt.fdr.level=0.1,
#                           rank.Mode = "Product")
