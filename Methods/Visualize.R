## Functions to visualze in IllTwo.R
Visual_JMirror <- function(p1,p2,
                           fdr.level.plot=c(0.05,0.1,0.2,0.3,0.4),
                           JM.Solver = NULL,
                           init.thred = 0.3,
                           rank.Mode = "Pmax",
                           H_type=NULL){
  m <- length(p1)
  cat("Conducting Joint Mirror Procedure...\n")
  if(is.null(JM.Solver)){
  JM.Res <- JointMirror.R(cbind(p1,p2),init.thred=init.thred,
                          offset=1,trgt.fdr.level=0,
                          rank.Mode = rank.Mode)
  JM.Solver <- JM.Res$JMirrorS
  }
  
  if(is.null(H_type)){
    H_type <- rep("Unknown",length(p1))
  }
  ## We can adjust FDP estimate and identify the fdr level
  
  FDP.est.val <- JM.Solver$getFDPest()
  In.Ind <- JM.Solver$getunMaskInd()+1
  In.Ind <- In.Ind[In.Ind!=m+1]
  Ex.Ind <- (1:m)[!(1:m %in% In.Ind)]
  FDP.est.In <- cummin(FDP.est.val[In.Ind])
  
  #FDP.est.In <- sapply(1:length(In.Ind), function(i){min(FDP.est.val[In.Ind[1:i]])})
  
  data.plot <- data.frame(p1 = p1[c(In.Ind,Ex.Ind)],
                          p2 = p2[c(In.Ind,Ex.Ind)],
                          H_type = H_type[c(In.Ind,Ex.Ind)],
                          FDP.est.val = c(FDP.est.In,rep(1,length(Ex.Ind))))
  data.plot$mask.p1 <- pmin(data.plot$p1,1-data.plot$p1)
  data.plot$mask.p2 <- pmin(data.plot$p2,1-data.plot$p2)
  
  ## Get information about contour
  data.plot.sub <- data.plot[1:length(In.Ind),]
  data.contour <- NULL
  
  cat("Pareparing dataset for ploting...\n")
  for(fdr.level.p in fdr.level.plot){
    if(rank.Mode=="Product"){
    check.mask.ind <- which(data.plot.sub$FDP.est.val<=fdr.level.p)
    check.p1.mask <- data.plot.sub$mask.p1[check.mask.ind]
    check.p2.mask <- data.plot.sub$mask.p2[check.mask.ind]
    reorder.ind <- order(check.p1.mask)
    check.p1.mask <- check.p1.mask[reorder.ind]
    check.p2.mask <- check.p2.mask[reorder.ind]
    
    check.p1.unmask <- data.plot.sub$mask.p1[-check.mask.ind]
    check.p2.unmask <- data.plot.sub$mask.p2[-check.mask.ind]
    check.p2.unmask <- check.p2.unmask[check.p1.unmask<=max(check.p1.mask)]
    check.p1.unmask <- check.p1.unmask[check.p1.unmask<=max(check.p1.mask)]
    un.reorder.ind <- order(check.p1.unmask)
    check.p1.unmask <- check.p1.unmask[un.reorder.ind]
    check.p2.unmask <- check.p2.unmask[un.reorder.ind]
    
    # Find the rejection region
    Support.up.p2 <- sapply(1:length(un.reorder.ind),function(i){
      if(min(check.p2.unmask[1:i])==check.p2.unmask[i]){
        return(c(check.p1.unmask[i],check.p2.unmask[i]))
      }
    })
    Support.up.p2 <- do.call("cbind",Support.up.p2)
    
    check.p1.mask.aug <- c(check.p1.mask,Support.up.p2[1,])
    check.p2.mask.aug <- c(check.p2.mask,Support.up.p2[2,])
    reorder.ind.aug <- order(check.p1.mask.aug)
    check.p1.mask.aug <- check.p1.mask.aug[reorder.ind.aug]
    check.p2.mask.aug <- check.p2.mask.aug[reorder.ind.aug]
    
    Support.up.p1.aug <- sapply(1:length(reorder.ind.aug),function(i){
      if(max(check.p2.mask.aug[i:length(reorder.ind.aug)])==check.p2.mask.aug[i]){
        return(c(check.p1.mask.aug[i],check.p2.mask.aug[i]))
      }
    })
    Support.up.p1.aug <- do.call("cbind",Support.up.p1.aug)
    Support.up.p <- Support.up.p1.aug
    #Support.up.p <- cbind(Support.up.p1,Support.up.p2)
    #Support.up.p <- Support.up.p[,order(Support.up.p[1,])]
    y.inter <- min(check.p2.unmask[check.p1.unmask<=min(check.p1.mask)])
    Support.up.p <- cbind(c(0,ifelse(is.infinite(y.inter),Support.up.p[2,1]*1.1,y.inter)),
                          Support.up.p,
                          c(min(
                            data.plot.sub$mask.p1[-check.mask.ind][data.plot.sub$mask.p1[-check.mask.ind]>max(check.p1.mask)]),
                            0))
    f <- approxfun(Support.up.p[1,],Support.up.p[2,])
    xx <- seq(0,max(Support.up.p[1,]),length.out=1000)
    data.contour <- rbind(data.contour,
                          data.frame(
                          p1=xx,
                          p2.up =f(xx),
                          p2.down=0,
                          fdr.level = fdr.level.p))
    
    }
    if(rank.Mode=="Pmax"){
      check.mask.ind <- which(data.plot.sub$FDP.est.val<fdr.level.p)
      check.p1.mask <- data.plot.sub$mask.p1[check.mask.ind]
      check.p2.mask <- data.plot.sub$mask.p2[check.mask.ind]
      xx <- seq(0,max(c(check.p1.mask,check.p2.mask)),length.out=1000)
      data.contour <- rbind(data.contour,
                            data.frame(
                              p1=c(xx,max(c(check.p1.mask,check.p2.mask))),
                              p2.up =c(rep(max(c(check.p1.mask,check.p2.mask)),1000),0),
                              p2.down=0,
                              fdr.level = fdr.level.p))
      
    }
    }
  
  data.contour$fdr.level <- sapply(data.contour$fdr.level,
                                   function(x){
                                     if(x==0.05){return("0-0.05")}
                                     if(x==0.1){return("0.05-0.1")}
                                     if(x==0.2){return("0.1-0.2")}
                                     if(x>0.3){return("0.2-0.3")}
                                     return(">0.3")})
  data.contour$fdr.level <- factor(data.contour$fdr.level,
                                   levels = c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3",">0.3"))
  data.contour$Region <- "R"
  data.contour.control1 <- data.frame(p1=1-data.contour$p1,
                                      p2.down=0,
                                      p2.up=data.contour$p2.up,
                                      fdr.level =  data.contour$fdr.level,
                                      Region = "A1")
  data.contour.control2 <- data.frame(p1=data.contour$p1,
                                      p2.down=1-data.contour$p2.up,
                                      p2.up=1,
                                      fdr.level =  data.contour$fdr.level,
                                      Region = "A2")
  data.contour <- rbind(data.contour,data.contour.control1,data.contour.control2)
  
  cat("Generating plots...\n")
  p <- ggplot(data=data.plot)+
    geom_ribbon(data=data.contour,aes(x=p1,
                                      ymin = p2.down, ymax = p2.up, 
                                      fill = Region,
                                      alpha = fdr.level)) +
    scale_alpha_manual(values=c(0.8, 0.4, 0.2,0.1,0))+
    scale_shape_manual(values=c(3,4,20,8))+
    scale_color_manual(values=c("#CAB2D6","#1F78B4","#33A02C","#E31A1C"))+
    scale_fill_manual(values=c("#A6CEE3","#B2DF8A","#FB9A99"))+
    geom_point(aes(x=p1,y=p2,
                   color=H_type,
                   shape=H_type),size=1)+
    geom_line(data=data.contour,aes(x=p1,
                                    y = p2.up,
                                    fill = Region,
                                    linetype = fdr.level),size=1) +
    geom_line(data=data.contour,aes(x=p1,
                                    y = p2.down,
                                    fill = Region,
                                    linetype = fdr.level),size=1)+
    #scale_shape(solid=F)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw()+
    theme(#legend.position = c("bottom"),
      legend.title = element_text(size=15), #change legend title font size
      legend.text = element_text(size=15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text=element_text(size=12))
  return(list(JM.Solver=JM.Solver,data.plot=data.plot,
              data.contour=data.contour,
              p=p))
}
