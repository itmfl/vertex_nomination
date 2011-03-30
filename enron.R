enron.inverse.rdpg <- function(ggg){
  tmp1 <- length(ggg)

  non.zero <- seq(0,0,length.out = tmp1)

  for(i in 1:tmp1){
    gggi <- ggg[[i]]
    if(sum(gggi) > 0)
      non.zero[i] = 1
  }

  non.zero.idx <- which(non.zero > 0)
  
  xconfig.list <- list()

  for(i in 1:length(non.zero.idx)){
    gggi <- ggg[[non.zero.idx[i]]]
    A1 <- (gggi == 1) + 0
    A2 <- (gggi == 2) + 0
    
    x1 <- svd.extract(A1)
    x2 <- svd.extract(A2)
    xconfig <- cbind(x1, x2)
    xconfig.list[[i]] <- xconfig
  }
  return(xconfig.list)
}

ggg2 <- list()

for(i in 2:13){
  ggg.tmp <- ggg[[4*i -3]]
  tmp1 <- ggg[[4*i - 2]]
  tmp1.idx1 <- (tmp1 == 1)
  tmp1.idx2 <- (tmp1 == 2)
  ggg.tmp[tmp1.idx1] <- 1
  ggg.tmp[tmp1.idx2] <- 2
  tmp1 <- ggg[[4*i - 1]]
  tmp1.idx1 <- (tmp1 == 1)
  tmp1.idx2 <- (tmp1 == 2)
  ggg.tmp[tmp1.idx1] <- 1
  ggg.tmp[tmp1.idx2] <- 2
 
  tmp1 <- ggg[[4*i]]
  tmp1.idx1 <- (tmp1 == 1)
  tmp1.idx2 <- (tmp1 == 2)
  ggg.tmp[tmp1.idx1] <- 1
  ggg.tmp[tmp1.idx2] <- 2

  ggg2[[i-1]] <- ggg.tmp
}

inverse.rdpg <- function(g, dim = 2){

  A.list <- adjacency.list(g$adjacency, g$labels)

  x.red <- svd.extract(A.list$red, k = dim)
  x.green <- svd.extract(A.list$green, k = dim)
  x.config <- cbind(x.red,x.green)

  return(x.config)
}

adjacency.list <- function(A, labels.mat){

    n <- nrow(A)

    A.red <- matrix(seq(0,0,length.out = n*n), nrow = n)
    A.green <- matrix(seq(0,0,length.out = n*n), nrow = n)

    red.idx <- (labels.mat == "red")
    green.idx <- (labels.mat == "green")

    A.red[red.idx] <- A[red.idx]
    A.green[green.idx] <- A[green.idx]

    return(list(red=A.red,green=A.green))
}

nonpsd.laplacian <- function(A){

    n = nrow(A)
    s <- apply(A,1,sum)

    L <- diag(s)/(n-1) + A
    return(L)
}

svd.extract <- function(A, k = 1){

    L <- nonpsd.laplacian(A)
    L.svd <- svd(L)
    
    L.svd.values <- L.svd$d[1:k]
    L.svd.vectors <- L.svd$v[,1:k]

    if(k == 1)
      L.coords <- L.svd.values * L.svd.vectors
    else
      L.coords <- L.svd.vectors %*% diag(L.svd.values)

    return(L.coords)
}




  



