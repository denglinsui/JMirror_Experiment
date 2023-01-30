set.seed(100)

## Provide an example showing the performance of our FDP estimator
## Including the true false discoviers(FD), FD we estimate, FD conventional joint significant test estimate
## As well as the modified FD.

library(kernelboot)
library(ggplot2)
library(patchwork)
library(cowplot)
m <- 10000
mu1H10 <- 1.5
mu1H11 <- 2
mu2H01 <- 2.5
mu2H11 <- 3
pi.seq <- c(0.4,0.2,0.2,0.2)
#pi.seq <- c(0,0.4,0.4,0.2)
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

## Quantile Quantile Plot
thred.seq <- seq(0,0.5,length.out=1000)
true.emp.fd1 <- function(tt){
  ind.exam <- which(p1<=tt & p2<=tt)
  sum(H_type[ind.exam] %in% c("H00","H01","H10"))
}

true.theo.fd1 <- function(tt){
  ind.exam <- 1:m
  sum(H_type[ind.exam] == "H00")*tt^2 +
    sum(H_type[ind.exam] == "H01")*tt*(pnorm(qnorm(tt/2),mu2H01)+
                                         1-pnorm(-qnorm(tt/2),mu2H01))+
    sum(H_type[ind.exam] == "H10")*tt*(pnorm(qnorm(tt/2),mu1H10)+
                                         1-pnorm(-qnorm(tt/2),mu1H10))
  
}

est.emp.fd1 <- function(tt){
  ind.exam <- which((p1>=1-tt & p2<=tt)|(p1<=tt & p2>=1-tt))
  sum(H_type[ind.exam] %in% c("H00","H01","H10"))
}

est.theo.fd1 <- function(tt){
  ind.exam <- 1:m
  2*sum(H_type[ind.exam] == "H00")*tt^2 +
    sum(H_type[ind.exam] == "H01")*tt*(pnorm(qnorm(tt/2),mu2H01)+
                                         1-pnorm(-qnorm(tt/2),mu2H01))+
    sum(H_type[ind.exam] == "H10")*tt*(pnorm(qnorm(tt/2),mu1H10)+
                                         1-pnorm(-qnorm(tt/2),mu1H10))+
    sum(H_type[ind.exam] == "H01")*tt*(-pnorm(qnorm((1-tt)/2),mu2H01)+
                                         pnorm(-qnorm((1-tt)/2),mu2H01))+
    sum(H_type[ind.exam] == "H10")*tt*(-pnorm(qnorm((1-tt)/2),mu1H10)+
                                         pnorm(-qnorm((1-tt)/2),mu1H10))
  #sum(H_type[ind.exam] %in% c("H00","H01","H10"))
}


true.emp.cfd1 <- function(tt){
  ind.exam <- which(p1<=tt & p2<=tt)
  2*sum(H_type[ind.exam] =="H00")+
    sum(H_type[ind.exam] =="H10")+
    sum(H_type[ind.exam] =="H01")
}

true.theo.cfd1 <- function(tt){
  ind.exam <- 1:m
  2*sum(H_type[ind.exam] == "H00")*tt^2 +
    sum(H_type[ind.exam] == "H01")*tt*(pnorm(qnorm(tt/2),mu2H01)+
                                         1-pnorm(-qnorm(tt/2),mu2H01))+
    sum(H_type[ind.exam] == "H10")*tt*(pnorm(qnorm(tt/2),mu1H10)+
                                         1-pnorm(-qnorm(tt/2),mu1H10))
  #sum(H_type[ind.exam] %in% c("H00","H01","H10"))
}

true.emp.fd1.val <- sapply(thred.seq,true.emp.fd1)
est.emp.fd1.val <- sapply(thred.seq,est.emp.fd1)
true.emp.cfd1.val <- sapply(thred.seq,true.emp.cfd1)
true.theo.cfd1.val <- sapply(thred.seq,true.theo.cfd1)
true.theo.fd1.val <- sapply(thred.seq,true.theo.fd1)
est.theo.fd1.val <- sapply(thred.seq,est.theo.fd1)

data.plot <- data.frame(thred.seq = thred.seq,
                        true.emp.fd1.val=true.emp.fd1.val,
                        est.emp.fd1.val=est.emp.fd1.val,
                        true.theo.fd1.val=true.theo.fd1.val,
                        nocomposit2=m*thred.seq*(1-pi.seq[4]),
                        est.theo.fd1.val=est.theo.fd1.val)
data.plot.theo.df <- data.table(t = thred.seq,
                                value = c(true.theo.fd1.val,
                                          m*thred.seq*(1-pi.seq[4]),
                                          est.theo.fd1.val),
                                `Theoretical False Discovery Number`= rep(c("False Discovery Number",
                                                                            "False Discovery Number Estimated by Joint Significant Test",
                                                                            "Control Number"),
                                                                          each = length(thred.seq)))
data.plot.emp.df <- data.frame(t = thred.seq,
                               value = c(true.emp.fd1.val,
                                         est.emp.fd1.val),
                               Empirical = rep(c("False Discovery Number",
                                                 "Control Number"),
                                               each = length(thred.seq)))

ggplot()+
  geom_point(data=data.plot.emp.df,aes(x=t,y=value,color=Empirical,shape=Empirical))+
  geom_line(data=data.plot.theo.df,aes(x=t,y=value,linetype=`Theoretical False Discovery Number`))+
  theme_bw()+
  #scale_shape_manual(values=)+
  theme(#panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal", 
    legend.box = "vertical",
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))



data.plot.theo.df <- data.table(t = thred.seq,
                                value = c(true.theo.fd1.val,
                                          m*thred.seq*(1-pi.seq[4]),
                                          est.theo.fd1.val),
                                `Theoretical`= rep(c("Number of False Discoveries",
                                                     "Number Estimated by Conventional Joint Significant Test",
                                                     "Number of Controls"),
                                                   each = length(thred.seq)))
data.plot.emp.df <- data.table(t = thred.seq,
                               value = c(true.emp.fd1.val,
                                         est.emp.fd1.val),
                               Empirical = rep(c("Number of False Discoveries",
                                                 "Number of Controls"),
                                               each = length(thred.seq)))

p.emp <- ggplot()+
  geom_point(data=data.plot.emp.df,aes(x=t,
                                       y=value,color=Empirical))+
  geom_line(data=data.plot.theo.df,aes(x=t,
                                       y=value,linetype=Theoretical))+
  theme_bw()+
  scale_color_manual(values=c("Number of False Discoveries"="#A6CEE3",
                              "Number of Controls"="#FB9A99"))+
  theme(#panel.grid.major = element_blank(),
    legend.position = "bottom",
    # legend.direction = "vertical", 
    # legend.box = "horizontal",
    legend.direction =  "horizontal",
    legend.box ="vertical", 
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.title = element_text(size=15), #change legend title font size
    legend.text = element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))
p.emp
ggsave("Figure/TheoreticalValueComparison.eps",height=8,width=7)
ggsave("Figure/TheoreticalValueComparison.pdf",height=8,width=7)





data.plot.cfd <- data.table(`Empirical Number of Modified False Discoveries` = true.emp.cfd1.val,
                            `Empirical Number of Controls` = est.emp.fd1.val)

p.plot.cfd <- ggplot()+
  geom_point(data=data.plot.cfd,aes(x=`Empirical Number of Modified False Discoveries`,
                                    y=`Empirical Number of Controls`))+
  geom_line(data=data.plot.cfd,aes(x=`Empirical Number of Modified False Discoveries`,
                                   y=`Empirical Number of Modified False Discoveries`),color="red")+
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))
p.plot.cfd
ggsave("Figure/TheoreticalValueComparisonCFD.eps",height=8,width=7)
ggsave("Figure/TheoreticalValueComparisonCFD.pdf",height=8,width=7)



grobs <- get_legend(p.emp)

pgrid <- plot_grid(p.emp+theme(legend.position = "none"),
                   p.plot.cfd,
                   ncol = 2,
                   label_size = 12,labels=c("A","B"))
p_grid <- plot_grid(pgrid, grobs, nrow = 2, rel_heights = c(1, .2))

save_plot("Figure/IllustrationApp.eps",p_grid,base_height=7,base_width=14)
save_plot("Figure/IllustrationApp.pdf",p_grid,base_height=7,base_width=14)
