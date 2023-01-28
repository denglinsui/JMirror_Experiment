#=================================================================================
# To implement MT.Comp in Huang, 2019
# We copy it from Liu et al., 2021 
#===============================================================================
int<-10
B<-10000; pdf<-NULL
for (i in 1:B){
  pdfi<-besselK(x=int*i/B, nu=0)
  pdf<-c(pdf, pdfi)
  #print(i); 
  flush.console()
}
#print(int)

myp<-function(cut){
  select<-(int*1:B/B)>cut
  pdf.sub<-pdf[select]
  pval<-sum(pdf.sub)/sum(pdf)
  return(pval)
}

MT_Comp<-function(a, b){
  ab<-a*b
  pp0<-sapply(abs(ab)/sqrt(1), myp)
  pp1<-sapply(abs(ab)/sd(a), myp)
  pp2<-sapply(abs(ab)/sd(b), myp)
  pp.comp<-pp1+pp2-pp0
  return(pp.comp)
}
