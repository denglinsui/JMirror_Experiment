#===== Plot example
library(ggplot2)
library(latex2exp)
library(data.table)
library(dplyr)
library(patchwork)
#library('igraph')
set.seed(100000)

num_p <- 7

## Generate \tilde{p}
p_x <- runif(num_p,0,0.5)
p_y <- runif(num_p,0,0.5)
rej_prob <- sample(seq(0.1,0.9,length.out=num_p))

data.DAG <- data.frame(p1=p_x,p2=p_y,q=rej_prob)

data.DAG <- data.DAG[order(-data.DAG$p1,data.DAG$p2),]

adj_matrix <- matrix(0,num_p,num_p)

# Generate adjacent matrix
for(i in 1:num_p){
  for(j in 1:num_p){
    if(data.DAG$p1[i]>data.DAG$p1[j] & data.DAG$p2[i]>data.DAG$p2[j]){
      adj_matrix[i,j] <- 1
    }
  }
}

data.DAG$indegree <- colSums(adj_matrix)
data.DAG$outdegree <- rowSums(adj_matrix)

adj_matrix <- adj_matrix[order(-data.DAG$outdegree),order(-data.DAG$outdegree)]
data.DAG <- data.DAG[order(-data.DAG$outdegree),]
data.DAG$label <- letters[1:num_p]
p_base <- ggplot(data.DAG,aes(x=p1,y=p2))+
  geom_point(color="lightblue",size=10)+
  geom_point(color="black",size=10,shape=1)+
  geom_text(aes(label=label),size=6)

# Transform adjacent matrix to dataframe
# Preparation for plotting
adj2df <- function(adj_matrix,data.DAG){
  datadf <- NULL
  num_p <-  dim(adj_matrix)[1]
  for(i in 1:num_p){
    for(j in 1:num_p){
      if(adj_matrix[i,j] ==1){
        datadf <- rbind(datadf,c(data.DAG$p1[i],data.DAG$p2[i],
                                 data.DAG$p1[j],data.DAG$p2[j]))
      }
    }
  }
  colnames(datadf) <- c("xstart","ystart","xend","yend")
  data.frame(datadf)
}

df.edge.full <- adj2df(adj_matrix,data.DAG)

#p_base

adj.dist <- 0.014
p_full <- ggplot(data.DAG,aes(x=p1,y=p2))+
  geom_segment(data=df.edge.full,aes(x=xstart,y=ystart,
                                     xend=xend+adj.dist,yend=yend+adj.dist),
               arrow = arrow(length=unit(0.30,"cm"),  type = "closed"))+
  geom_point(aes(color=q),size=10)+
  geom_point(color="black",size=10,shape=1)+
  geom_text(aes(label=label),size=6)+
  theme_bw()+
  xlab(TeX("$\\tilde{p}_1$"))+ylab(TeX("$\\tilde{p}_2$"))+
  guides(colour = guide_colorbar(barwidth = unit(10, "cm"))) +
  theme(#panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size=15), #change legend title font size
    legend.text = element_text(size=15),
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'))

# transitive reduction
df.edge.red <- adj2df(adj_matrix-(adj_matrix%*%adj_matrix>0),data.DAG)              
data.DAG$unmaskorder <- c(1,2,3,4,6,7,5)
p_red <- ggplot(data.DAG,aes(x=p1,y=p2))+
  geom_segment(data=df.edge.red,aes(x=xstart,y=ystart,
                                    xend=xend+adj.dist,yend=yend+adj.dist),
               arrow = arrow(length=unit(0.30,"cm"),  type = "closed"))+
  geom_point(aes(color=q),size=10)+
  geom_point(color="black",size=10,shape=1)+
  geom_text(aes(label=label),size=6)+
  geom_label(aes(x=p1+0.025,y=p2-0.005,label=unmaskorder),size=5)+
  xlab(TeX("$\\tilde{p}_1$"))+ylab(TeX("$\\tilde{p}_2$"))+
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'))

# transitive reduction at step four
p_red_step4 <- ggplot(data.DAG,aes(x=p1,y=p2))+
  geom_segment(data=data.table(df.edge.red)%>% filter(ystart<=0.18),aes(x=xstart,y=ystart,
                                                                        xend=xend+adj.dist,yend=yend+adj.dist),
               arrow = arrow(length=unit(0.30,"cm"),  type = "closed"))+
  geom_point(aes(color=q),size=10)+
  geom_point(color="gray",data = data.DAG[data.DAG$unmaskorder<4,],size=10)+
  geom_point(color="black",size=10,shape=1)+
  xlab(TeX("\\tilde{p}_1"))+ylab(TeX("\\tilde{p}_2"))+
  geom_text(aes(label=label),size=6)+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'))+
  theme_bw()+
  xlab(TeX("$\\tilde{p}_1$"))+ylab(TeX("$\\tilde{p}_2$"))+
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12))+
  geom_label(aes(x=p1+0.025,y=p2-0.005,label=unmaskorder),size=5)

## Arrange three plots
grobs <- get_legend(p_full)

pgrid <- plot_grid(p_full+theme(legend.position = "none"),
                   p_red+theme(legend.position = "none",
                               axis.ticks.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.title.y=element_blank()),
                   p_red_step4+theme(legend.position = "none",
                                     axis.ticks.y=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.title.y=element_blank()),
                   ncol = 3,
                   label_size = 12,labels=c("A","B","C"))
p_grid <- plot_grid(pgrid, grobs, nrow = 2, rel_heights = c(1, .12))


save_plot("Figure/DAG.eps",p_grid,base_height=5,base_width=13)
save_plot("Figure/DAG.pdf",p_grid,base_height=5,base_width=13)
