#===================================================
# Conduct Z test statistic DACT in Liu et al., 2021 
#===================================================
DACT.z_stat <- function (p_a, p_b,Z_a=NULL,Z_b=NULL, correction = NULL) {
  if(is.null(Z_a)|is.null(Z_b)){
    Z_a = stats::qnorm(p_a, lower.tail = F)
    Z_b = stats::qnorm(p_b, lower.tail = F)
  }
  #pi0a = nonnullPropEst(Z_a, 0, 1)
  #pi0b = nonnullPropEst(Z_b, 0, 1)
  # What DACT estimate is non null proportion
  pia = nonnullPropEst(Z_a, 0, 1)
  pib = nonnullPropEst(Z_b, 0, 1)
  
  if (pia > 1) {
    pia = 1
  }
  if (pib > 1) {
    pib = 1
  }
  pi0a = 1-pia
  pi0b = 1-pib
  
  if (pi0a > 1) {
    pi0a = 1
  }
  if (pi0b > 1) {
    pi0b = 1
  }
  p.mat = cbind(p_a, p_b)
  p3 = (apply(p.mat, 1, max))^2
  wg1 = pi0a * (1 - pi0b)
  wg2 = (1 - pi0a) * pi0b
  wg3 = pi0a * pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1, wg2, wg3)/wg.sum
  p_dact = wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * 
    p3
  if(!is.null(correction)){
    
    if (correction == "Efron") {
      p_dact = EfronCorrect(p_dact)
    }
    if (correction == "JC") {
      p_dact = JCCorrect(p_dact)
    }
  }
  return(p_dact)
}

nullParaEst<-function (x,gamma=0.1)
{
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation
  
  n = length(x)
  t = c(1:1000)/200
  
  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)
  
  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]
  
  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)
  
  return(musigma=list(mu=uhat,s=shat))
}

nonnullPropEst <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) {
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}


EfronCorrect = function(pval){
  pval[pval==0] <- min(pval[pval!=0])
  z = stats::qnorm(1-pval)
  z[is.infinite(z)] <- ifelse(z[is.infinite(z)]>0,max(z[!is.infinite(z)]),min(z[!is.infinite(z)]))
  res <- locfdr(z,nulltype = 1)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}

JCCorrect = function(pval){
  pval[pval==0] <- min(pval[pval!=0])
  z = stats::qnorm(1-pval)
  z[is.infinite(z)] <- ifelse(z[is.infinite(z)]>0,max(z[!is.infinite(z)]),min(z[!is.infinite(z)]))
  #z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}


