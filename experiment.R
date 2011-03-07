source("kidney_egg.R")
source("nominate_vertex.R")

validate.graph <- function(n,m,l,p.vec,s.vec,num.iter){

    ## kidney.red <- c(0.4,0.2,0)
    ## kidney.green <- c(0.4,0.2,0)
    ## egg.green <- c(0.4,0.2,0)
    ## egg.red <- c(0.4,0.2,sqrt(s1 - 0.2))

    stux <- matrix(0,nrow = num.iter, ncol = n)
    stun <- stux

    for(i in 1:num.iter){
        g.i <- kidney.egg(n,m,l,p.vec,s.vec) 
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
        
glen.driver1 <- function(n,m,l,p.vec,s.vec,
                         num.monte.carlo.iter,
                         method = "inverse.rdpg"){

  n <- 184 
  
  minR.noties <- seq(n-m,0,length.out = num.monte.carlo.iter)
  minR.ties <- minR.noties
  
  for(i in 1:num.monte.carlo.iter){
    g.i <- kidney.egg(n,m,l,p.vec,s.vec)
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

## The code below run the kidney.egg.rdpg to generate random instances
## of kidney and egg graph. However, it is not necessary for our purpose.
glen.driver2 <- function(n,m,l,s1,num.monte.carlo.iter, method = "inverse.rdpg"){

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

## tn.stats <- NULL
## tx.stats <- NULL
## tf.stats <- NULL
## rdpg.stats <- NULL
## diff.stats <- NULL

## tx.stats1 <- glen.driver1(184,24,12,0.25,1000,method="tx.statistics")
## tn.stats1 <- glen.driver1(184,24,12,0.25,1000,method="tn.statistics")
## tf.stats1 <- glen.driver1(184,24,12,0.25,1000,method="tf.statistics")
## rdpg.stats1 <- glen.driver1(184,24,12,0.25,1000,method="inverse.rdpg")
## diff.stats1 <- glen.driver1(184,24,12,0.25,1000,method="diffusion.distance")

## tx.stats2 <- glen.driver1(184,24,12,0.40,1000,method="tx.statistics")
## tn.stats2 <- glen.driver1(184,24,12,0.40,1000,method="tn.statistics")
## tf.stats2 <- glen.driver1(184,24,12,0.40,1000,method="tf.statistics")
## rdpg.stats2 <- glen.driver1(184,24,12,0.40,1000,method="inverse.rdpg")
## diff.stats2 <- glen.driver1(184,24,12,0.40,1000,method="diffusion.distance")

t1.test <- function(n,m,l,p.vec,s.vec, num.mc = 1000){

    p1 <- p.vec[1]; p2 <- p.vec[2]
    s1 <- s.vec[1]; s2 <- s.vec[2]

    rk1 <- seq(0,0,length.out = num.mc)
    rk2 <- seq(0,0,length.out = num.mc)

    for( i in 1:num.mc){
        a1 <- rbinom(m-l,l,(s1 + s2))
        a2 <- rbinom(n-m,l,(p1 + p2))
        a3 <- c(a1,a2)
        
        a3.sorted <- sort(a3, decreasing = TRUE, index.return = TRUE)
        rk1[i] <- which(a3.sorted$ix <= (m-l))[1]

        val.rk1 <- a3.sorted$x[rk1[i]]
        val.idx <- which(a3 == val.rk1)
        red.length <- length(which(val.idx <= (m-l)))
        green.length <- length(val.idx) - red.length
        rk2[i] <- rk1[i] + green.length/(1 + red.length)
    }
    return(list(nobreak.ties = mean(rk1),break.ties=mean(rk2)))
}

t2.test <- function(n,m,l,p.vec,s.vec,num.mc = 1000){

    p1 <- p.vec[1]; p2 <- p.vec[2]
    s1 <- s.vec[1]; s2 <- s.vec[2]

    rk1 <- seq(0,0,length.out = num.mc)
    rk2 <- seq(0,0,length.out = num.mc)

    for( i in 1:num.mc){
        a1 <- rbinom(m-l,m-1,s1) + rbinom(m-l,n-m,p1)
        a2 <- rbinom(n-m,n-1,p1)
        a3 <- c(a1,a2)

        a3.sorted <- sort(a3, decreasing = TRUE, index.return = TRUE)
        rk1[i] <- which(a3.sorted$ix <= (m-l))[1]
        
        val.rk1 <- a3.sorted$x[rk1[i]]
        val.idx <- which(a3 == val.rk1)
        red.length <- length(which(val.idx <= (m-l)))
        green.length <- length(val.idx) - red.length
        rk2[i] <- rk1[i] + green.length/(1 + red.length)
    }
    return(list(nobreak.ties = mean(rk1),break.ties = mean(rk2)))
}

nhlee.driver1 <- function(directory, method = "inverse.rdpg", n = 12, m = 5, l = 2){

  filelist <- list.files(".", pattern = "myTVsK*")
  num.mc <- length(filelist)
  ##num.mc <- 200
  
  minR.noties <- seq(0,0,length.out = num.mc)
  minR.ties <- minR.noties
  prob.vec <- minR.noties


  for(i in 1:num.mc){
    g.i <- parse.tvks(filelist[i])
    xyz.i <- nominate.vertex(g.i, method)
  
    minR.noties[i] <- which(xyz.i$order <=m)[1]
    
    minR.noties.val <- xyz.i$value[minR.noties[i]]
    minR.noties.val.idx <- xyz.i$order[which(xyz.i$value == minR.noties.val)]
    k1 <- length(which(minR.noties.val.idx > m))
    k2 <- length(minR.noties.val.idx) - k1

    minR.ties[i] <- minR.noties[i] + k1/(1 + k2)
    prob.vec[i]<- k2/(k1 + k2)
  }

  return(list(prob.vec=prob.vec,
              prob=sum(prob.vec)/num.mc,
              vector.ties=minR.ties,
              values.ties = sum(minR.ties)/num.mc))
        
}
tn1 <- nhlee.driver1(".", method = "tn.statistics")
tx1 <- nhlee.driver1(".", method = "tx.statistics")
tf1 <- nhlee.driver1(".", method = "tf.statistics")
inverse.rdpg <- nhlee.driver1(".", method = "inverse.rdpg")


  
glen.driver1 <- function(n,m,l,p.vec,s.vec,
                         num.monte.carlo.iter,
                         method = "inverse.rdpg"){

  n <- 184 
  
  minR.noties <- seq(n-m,0,length.out = num.monte.carlo.iter)
  minR.ties <- minR.noties
  
  for(i in 1:num.monte.carlo.iter){
    g.i <- kidney.egg(n,m,l,p.vec,s.vec)
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

