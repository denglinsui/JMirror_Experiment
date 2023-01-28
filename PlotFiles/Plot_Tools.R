## Additional tools for plot

ylab.change <- function(x){
  if(x=="AdjFDP"){return("cFDP")}
  if(x=="POW"){return("Percentage of Discoveries")}
  if(x=="FDP"){return("FDP")}
  return(x)
}

break.num <- 5
scales_y_MedAna <- list(
  "q * {\n    phantom() == phantom()\n} * 0.05" = 
    scale_y_continuous(breaks = round(seq(0,0.1,length.out=break.num),2),
                       labels = c(round(seq(0,0.1,length.out=break.num),2)[1:(break.num-1)],"0.1~1")),
  "q * {\n    phantom() == phantom()\n} * 0.2" = 
    scale_y_continuous(breaks = round(seq(0,0.4,length.out=break.num),2),
                       labels = c(round(seq(0,0.4,length.out=break.num),2)[1:(break.num-1)],"0.4~1"))
)
MedAna.color.Panel <- c("JM.Product.AdaKer"="#FB9A99",
                        "JM.Product"="#E31A1C",
                        "JM.NoShape.AdaKer"="#A6CEE3",
                        #"JM.NoShape"= "#1F78B4",
                        "JM.EmptyPoset"= "#1F78B4",
                        #"JM.Pmax"="#B15928",
                        "JM.Max"="#B15928",
                        "JS.Mix.Asy"="#B2DF8A",
                        "DACT"=  "#6A3D9A" ,
                        "MT.Comp"="#FDBF6F",
                        "JS.Mix.Finite"="#33A02C" )


MedAna.shape.Panel <- c("JM.Product"=2,
                        "JM.Product.AdaKer"=17,
                       # "JM.NoShape"=5,
                        "JM.EmptyPoset"=5,
                        "JM.NoShape.AdaKer"= 18,
                        #"JM.Pmax"=8,
                        "JM.Max"=8,
                        "JS.Mix.Finite"=20,
                        "JS.Mix.Asy"=1,
                        "DACT"= 3 ,
                        "MT.Comp"=4) 


theme_Med <- function(){
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="bottom",
    legend.title = element_text(size=15), #change legend title font size
    legend.text = element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))
}


RepAna.color.Panel <- c("JM.Product.AdaKer"="#FB9A99",
                        "JM.Product"="#E31A1C",
                        "JM.NoShape.AdaKer"="#A6CEE3",
                        #"JM.NoShape"= "#1F78B4",
                        "JM.EmptyPoset"= "#1F78B4",
                        #"JM.Pmax"="#B15928",
                        "JM.Max"="#B15928",
                        "JS_mix_asy"="#B2DF8A",
                        "DACT"=  "#6A3D9A" ,
                        "MT_Comp"="#FDBF6F",
                        "JS_mix_finite"="#33A02C" )


RepAna.shape.Panel <- c("JM.Product"=2,
                        "JM.Product.AdaKer"=17,
                        #"JM.NoShape"=5,
                        "JM.EmptyPoset"=5,
                        "JM.NoShape.AdaKer"= 18,
                        #"JM.Pmax"=8,
                        "JM.Max"=8,
                        "JS_mix_finite"=20,
                        "JS_mix_asy"=1,
                        "DACT"= 3 ,
                        "MT_Comp"=4 
)

theme_Rep <- function(){
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="bottom",
    legend.title = element_text(size=15), #change legend title font size
    legend.text = element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 12),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))
}


RepAna.color.Panel <- c("JM.Product"="#E31A1C",
                        #"JM.NoShape"= "#1F78B4",
                        "JM.EmptyPoset"= "#1F78B4",
                        #"JM.Pmax"="#B15928",
                        "JM.Max"="#B15928",
                        "BonferroniBH"="#B2DF8A",
                        "adaFilter"=  "#CAB2D6" ,
                        "ssa"="#FDBF6F",
                        "cPC"="#F781BF" )


Transfer_pi <- function(s_pi){
  if(all(s_pi==c(1,0,0,0))){
    return("GNull")
  }
  if(all(s_pi==c(0.6,0.2,0.2,0))){
    return("DNull")
  }
  if(all(s_pi==c(0.9,0.05,0.05,0))){
    return("SNull")
  }
  if(all(s_pi==c(0.4,0.2,0.2,0.2))){
    return("DAlter")
  }
  if(all(s_pi==c(0.88,0.05,0.05,0.02))){
    return("SAlter")
  }
  if(all(s_pi==c(0.9,0,0,0.1))){
    return("H00Alt")
  }
  if(all(s_pi==c(0.4,0,0,0.6))){
    return("H00DAlt")
  }
  if(all(s_pi==c(0.88,0,0,0.12))){
    return("H00SAlt")
  }
}


SummaryTable <- function(M1,M2,SE1,SE2){
  col.order <- c("alter_type","JS.Mix.Asy","JS.Mix.Finite","DACT",
                 "MT.Comp","JM.Max","JM.Product","JM.EmptyPoset")
  M1 <- M1[,col.order]
  M2 <- M2[,col.order]
  SE1 <- SE1[,col.order]
  SE2 <- SE2[,col.order]
  MFull <- t(cbind(M1,M2[,-1]))
  SEFull <- t(cbind(SE1,SE2[,-1]))
  
  OutFull <- MFull
  
  for(i in 2:dim(OutFull)[1]){
    for(j in 1:dim(OutFull)[2]){
      OutFull[i,j] <- sprintf("%.3f(%.3f)",as.numeric(MFull[i,j]),
                              as.numeric(SEFull[i,j])*10)
    }
  }
  OutFull <- cbind("",rownames(OutFull),OutFull)
  
  colnames(OutFull) <- OutFull[1,]
  OutFull <- OutFull[-1,]
  OutFull
}
