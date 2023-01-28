### Additional tools for real data analysis

Intersect.list <- function(Y){
  y <- intersect(Y[[1]],Y[[2]])
  for(ind in 3:(length(Y))){
    y <- intersect(y,Y[[ind]])
  }
  return(y)
}

Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
}

## color

# RColorBrewer::display.brewer.pal(9,"Set1")
# RColorBrewer::brewer.pal(9,"Set1")
# RColorBrewer::display.brewer.pal(8,"Pastel2")
# RColorBrewer::brewer.pal(8,"Pastel2")

NAS.color.Panel <- c("DACT+JM.Max+JM.EmptyPoset+JM.Product"="#377EB8",
                     "DACT+JM.EmptyPoset+JM.Product"="#A65628",
                     "DACT+JM.Max+JM.Product"="#F781BF",
                     "JM.EmptyPoset+JM.Product"="#E41A1C",
                     "JM.EmptyPoset+JM.Max"= "#FF7F00",
                     "JM.EmptyPoset"="#4DAF4A",
                     "Null"="#CCCCCC")

NAS.fill.Panel <- c( "JM.Product"="#E31A1C",
                        "JM.Max"="#1F78B4")


NAS.shape.Panel <- c("DACT+JM.Max+JM.EmptyPoset+JM.Product"=18,
                     "DACT+JM.EmptyPoset+JM.Product"=15,
                     "DACT+JM.Max+JM.Product"=17,
                     "JM.EmptyPoset+JM.Product"=4,
                     "JM.EmptyPoset+JM.Max"= 3,
                     "JM.EmptyPoset"=2,
                     "Null"=20)

NAS.size.Panel <- c("DACT+JM.Max+JM.EmptyPoset+JM.Product"=2,
                     "DACT+JM.EmptyPoset+JM.Product"=2,
                     "DACT+JM.Max+JM.Product"=2,
                     "JM.EmptyPoset+JM.Product"=2,
                     "JM.EmptyPoset+JM.Max"= 2,
                     "JM.EmptyPoset"=2,
                     "Null"=1)

themeRealData <- function(){
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
  