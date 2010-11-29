source("kidney_egg.R")
source("nominate_vertex.R")

validate.graph <- function(n,m,l,s1,num.iter){

    kidney.red <- c(0.4,0.2,0)
    kidney.green <- c(0.4,0.2,0)
    egg.green <- c(0.4,0.2,0)
    egg.red <- c(0.4,0.2,sqrt(s1 - 0.2))

    stux <- matrix(0,nrow = num.iter, ncol = n)
    stun <- stux

    for(i in 1:num.iter){
        g.i <- kidney.egg.rdpg(n,m,l,kidney.red,kidney.green,egg.red,egg.green)

        gi.stat <- g.statistics(g.i)

        stux[i,] <- gi.stat$t1
        stun[i,] <- gi.stat$t2

    }
    stux.mean <- apply(stux,2,sum)/num.iter
    stux.sd <- apply(stux,2,sd)/num.iter
    stun.mean <- apply(stun,2,sum)/num.iter
    stun.sd <- apply(stun,2,sd)/num.iter
    
    return(list(t1.mean = stux.mean, t1.sd = stux.sd,
                t2.mean = stun.mean, t2.sd = stun.sd))
}
        
        


glen.driver1 <- function(n,m,l,s1,num.monte.carlo.iter, method = "inverse.rdpg"){

  n <- 184 
  kidney.red <- c(0.4,0.2,0)
  kidney.green <- c(0.4,0.2,0)
  egg.green <- c(0.4,0.2,0)
  egg.red <- c(0.4,0.2,sqrt(s1 - 0.2))

  minR.noties <- seq(n-m,0,length.out = num.monte.carlo.iter)
  minR.ties <- minR.noties
  
  for(i in 1:num.monte.carlo.iter){
    g.i <- kidney.egg.rdpg(n,m,l,kidney.red,kidney.green,egg.red,egg.green)
    xyz.i <- nominate.vertex(g.i, method)

    minR.noties[i] <- which(xyz.i$order <=m)[1]
    
    minR.noties.val <- xyz.i$value[minR.noties[i]]
    minR.noties.val.idx <- xyz.i$order[which(xyz.i$value == minR.noties.val)]
    k1 <- length(which(minR.noties.val.idx > m))
    k2 <- length(minR.noties.val.idx) - k1

    minR.ties[i] <- minR.noties[i] + k1/(1 + k2)
  }
  
  return(list(vector.noties=minR.noties,
              value.noties=sum(minR.noties)/num.monte.carlo.iter,
              vector.ties=minR.ties,
              values.ties = sum(minR.ties)/num.monte.carlo.iter))
}
