## Basic

#=== fdp & pow

adjfdp <- function(selected,Beta){
  if(length(selected)==0){
    adjfdp <- 0
  }else{
    NumZero <- apply(Beta,1,function(x){sum(x==0)})
    adjfdp <- sum(NumZero[selected])/length(selected)}
  return(adjfdp)
}

fdp <- function(selected,H0){
  if(length(selected)==0){
    fdp <- 0
  }else{
    fdp <- sum(H0[selected]==1)/length(selected)}
  return(fdp)
}
Pow <- function(selected,H0){
  if(sum(1-H0)==0){
    Pow <- 0
  }else{
    Pow <- sum(H0[selected]==0)/sum(1-H0)
  }
  return(Pow)
}

# break ties for p-values: to make sure the method in dai2020 can work
Break_ties <- function(x){
  x.order <- order(x)
  x.Ind <- which(diff(sort(x))==0)
  
  if(length(x.Ind)==0){
    x[x==0] <- runif(sum(x==0),0,min(x[x!=0]))
    return(x)
  }
  Sep <- list()
  X.Diff <- diff(x.Ind);X.Diff <- c(1,X.Diff,Inf)
  i <- 1
  Wait.seq <- c()
  while(i<=(length(X.Diff))){
    if(X.Diff[i]>1){
        Wait.seq <- c(Wait.seq,max(Wait.seq)+1)
        x[x.order[Wait.seq]] <- runif(length(Wait.seq),
                                     x[x.order[Wait.seq[1]]],x[x.order[max(Wait.seq)+1]])
      # Start the next sequence
      Wait.seq <- c(x.Ind[i])
    }else{
      Wait.seq <- c(Wait.seq,x.Ind[i])
    }
    i <- i+1
  }
  
  # remove zero
  x[x==0] <- runif(sum(x==0),0,min(x[x!=0]))
  x
}
