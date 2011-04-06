kidney.egg.attributed <- function(n,m,l,pi0,piA){
  
  pi00 <- pi0*pi0
  pi0A <- pi0*piA; piA0 <- pi0A
  piAA <- piA*piA
  
  runif.mat <- matrix(0,nrow = n, ncol = n)
  runif.mat[col(runif.mat) > row(runif.mat)] <- runif(n*(n-1)/2)
  runif.mat <- (runif.mat + t(runif.mat))

  edge.mat <- matrix(0, nrow = n, ncol =n)
  edge.mat[1:m,1:m] <- sum(piAA)
  edge.mat[1:m,(m+1):n] <- sum(pi0A) 
  edge.mat[(m+1):n,1:m] <- sum(pi0A)
  edge.mat[(m+1):n,(m+1):n] <- sum(pi00)

  red.mat <- matrix(0, nrow = n, ncol = n)
  red.mat[1:m,1:m] <- piAA[1]
  red.mat[1:m,(m+1):n] <- pi0A[1]
  red.mat[(m+1):n,1:m] <- pi0A[1]
  red.mat[(m+1):n,(m+1):n] <- pi00[1]

  adjacency.mat <- matrix(0,nrow = n, ncol = n)
  red.idx <- (runif.mat <= red.mat)
  green.idx <- (runif.mat > red.mat & runif.mat <= edge.mat)

  adjacency.mat[red.idx] <- 1
  adjacency.mat[green.idx] <- 2
  diag(adjacency.mat) <- 0

  vertex.colors <- c(rep("red",length.out=l),
                     rep("blue", length.out = m-l),
                     rep("green",length.out=(n-m)))

  return(list(adjacency=adjacency.mat,
              v.colors = vertex.colors))
}
dirichlet.rdpg <- function(n,m,l, x0, xA){

  pi00.1 <- x0[,1] %o% x0[,1]
  pi00.2 <- x0[,2] %o% x0[,2]
  pi0A.1 <- x0[,1] %o% xA[,1]
  pi0A.2 <- x0[,2] %o% xA[,2]
  piA0.1 <- t(pi0A.1)
  piA0.2 <- t(pi0A.2)
  piAA.1 <- xA[,1] %o% xA[,1]
  piAA.2 <- xA[,2] %o% xA[,2]

  pi00 <- pi00.1 + pi00.2
  pi0A <- pi0A.1 + pi0A.2
  piA0 <- pi0A
  piAA <- piAA.1 + piAA.2

  edge.mat <- matrix(0, nrow = n, ncol = n)
  edge.mat[1:m,1:m] <- piAA 
  edge.mat[1:m,(m+1):n] <- piA0
  edge.mat[(m+1):n,1:m] <- pi0A
  edge.mat[(m+1):n,(m+1):n] <- pi00
  
  runif.mat <- matrix(0,nrow = n, ncol = n)
  runif.mat[col(runif.mat) > row(runif.mat)] <- runif(n*(n-1)/2)
  runif.mat <- (runif.mat + t(runif.mat))
  
  red.mat <- matrix(0, nrow = n, ncol = n)
  red.mat[1:m,1:m] <- piAA.1
  red.mat[1:m,(m+1):n] <- piA0.1
  red.mat[(m+1):n,1:m] <- pi0A.1
  red.mat[(m+1):n,(m+1):n] <- pi00.1
  
  adjacency.mat <- matrix(0,nrow = n, ncol = n)
  red.idx <- (runif.mat <= red.mat)
  green.idx <- (runif.mat > red.mat & runif.mat <= edge.mat)

  adjacency.mat[red.idx] <- 1
  adjacency.mat[green.idx] <- 2
  diag(adjacency.mat) <- 0

  vertex.colors <- c(rep("red",length.out=l),
                     rep("blue", length.out = m-l),
                     rep("green",length.out=(n-m)))

  return(list(adjacency=adjacency.mat,
              v.colors = vertex.colors))
}

adjacency.list <- function(A){

    A.red <- (A == 1) + 0
    A.green <- (A == 2) + 0

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

count.red.nn <- function(dist.mat, colors, k){

    idx <- which((colors == "blue") | (colors == "green"))
    num.red.nn <- numeric(length(idx))

    for(i in 1:length(idx)){
        dist.mat.i <- dist.mat[idx[i],]
        dist.mat.i.sorted <- sort(dist.mat.i, decreasing=FALSE,
                                  index.return = TRUE)

        k.nn.idx <- dist.mat.i.sorted$ix[2:(k+1)]
        num.red.nn[i] <- sum(colors[k.nn.idx] == "red")
    }

    sorted.red.nn <- sort(num.red.nn, decreasing = TRUE, index.return = TRUE)
    
    return(list(vector = sorted.red.nn$x, order = idx[sorted.red.nn$ix]))
}

dave.lda.nominate <- function(g, dim = 2){
  A.list <- adjacency.list(g$adjacency)
  x.red <- svd.extract(A.list$red, k = dim)
  x.green <- svd.extract(A.list$green, k = dim)
  x.separate <- cbind(x.red,x.green)

  gr <- c(rep(1, length(which((g$v.colors == "red")))),
          rep(2, length(which((g$v.colors == "blue") | (g$v.colors == "green")))))
  
  abc <- lda(x.separate, grouping = gr, CV = TRUE)
  idx <- which(g$v.colors != "red")
  
  uvw <- sort(abc$posterior[idx,1], decreasing = TRUE, index.return = TRUE)

  return(list( order = idx[uvw$ix], value = uvw$x))
}

carey.lda.nominate <- function(g, dim = 2){
  A.list <- adjacency.list(g$adjacency)
  x.red <- svd.extract(A.list$red, k = dim)
  x.green <- svd.extract(A.list$green, k = dim)
  x.separate <- cbind(x.red,x.green)

  gr <- c(rep(1, length(which((g$v.colors == "red") | (g$v.colors == "blue")))),
            rep(2, length(which((g$v.colors == "green")))))
  
  abc <- lda(x.separate, grouping = gr, CV = TRUE)
  idx <- which(g$v.colors != "red")
  
  uvw <- sort(abc$posterior[idx,1], decreasing = TRUE, index.return = TRUE)

  return(list( order = idx[uvw$ix], value = uvw$x))
}

p.mrr.nsrr <- function(xyz, idx){
  maxval <- max(xyz$value)
  maxidx <- which(xyz$value == maxval)

  k1 <- length(which(xyz$order[maxidx] <= m))
  p <- k1/length(maxidx)
  
  rk <- which(xyz$order <= m)[1]
  rk.val <- xyz$value[rk]
  rk.val.idx <- xyz$order[which(xyz$value == rk.val)]
  k1 <- length(which(rk.val.idx > m))
  k2 <- length(rk.val.idx) - k1

  mrr <- 1/(rk + k1/(1 + k2))

  return(list(p = p, mrr = mrr))
}


count.red.nn.nominate <- function(g, dim = 2){

    A.list <- adjacency.list(g$adjacency)

    x.red <- svd.extract(A.list$red, k = dim)
    x.green <- svd.extract(A.list$green, k = dim)
    x.separate <- cbind(x.red,x.green)
    
    d.mat <- as.matrix(dist(x.separate))
    xyz <- count.red.nn(d.mat, g$v.colors,
                        k = sum(g$v.colors == "red"))

    return(list(order=xyz$order, value=xyz$vector))
}

kidney.egg.rdpg.driver <- function(n,m,l,pi0, piA, mc, method = "count.red.nn"){

  prob.vec <- numeric(mc)
  mrr.vec <- numeric(mc)

  for(i in 1:mc){
    g.i <- kidney.egg.attributed(n,m,l,pi0,piA)
    xyz.i <- do.call(paste(method,".nominate", sep = ""), list(g.i,2))
##    xyz.i <- inverse.rdpg.nominate(g.i)

    ttt <- p.mrr.nsrr(xyz.i, m)
    mrr.vec[i] <- ttt$mrr
    prob.vec[i] <- ttt$p
  }
  return(list(prob.vec = prob.vec, p = mean(prob.vec), mrr.vec = mrr.vec, mrr = mean(mrr.vec)))
}

dirichlet.rdpg.driver <- function(n,m,l,alpha0,alphaA,mc, method = "count.red.nn"){

  prob.vec <- numeric(mc)
  rk.vec <- numeric(mc)

  for(i in 1:mc){
    x0 <- rdirichlet(n - m, r*alpha0 + 1)
    xA <- rdirichlet(m, r*alphaA + 1)

    x0 <- x0[,2:3]
    xA <- xA[,2:3]
    
    g.i <- dirichlet.rdpg(n,m,l,x0,xA)
    xyz.i <- do.call(paste(method,".nominate", sep = ""), list(g.i,2))
    
    ## xyz.i <- inverse.rdpg.nominate(g.i)

    rk <- which(xyz.i$order <= m)[1]
    rk.val <- xyz.i$value[rk]
    rk.val.idx <- xyz.i$order[which(xyz.i$value == rk.val)]
    k1 <- length(which(rk.val.idx > m))
    k2 <- length(rk.val.idx) - k1

    rk.vec[i] <- rk + k1/(1 + k2)
    
    maxval <- max(xyz.i$value)
    maxidx <- which(xyz.i$value == maxval)
    idx <- xyz.i$order[maxidx]

    k1 <- length(which(idx <= m))
    prob.vec[i] <- k1/length(idx)
    if(is.nan(prob.vec[i]))
      browser()
  }
  return(list(prob.vec = prob.vec, p = sum(prob.vec)/mc, rk.vec = rk.vec, mrr = mean(1/rk.vec)))
}

cep1 <- function(n, method = "count.red.nn", seed = 1729){
  set.seed(seed)
  m = 10
  l = 5
  r = 100
  mc <- 1000
  
  t <- seq(0.2,0.8,by=0.1)
  pvec <- numeric(length(t))
  mrr.vec <- numeric(length(t))

  alpha0 <- c(0.6,0.2,0.2)

  for(i in 1:length(t)){
    alphaA <- c(0.8 - t[i], t[i], 0.2)

    abc <- dirichlet.rdpg.driver(n,m,l,alpha0,alphaA,mc, method)
    pvec[i] <- abc$p
    mrr.vec[i] <- abc$mrr
  }
  return(list(p = pvec, mrr = mrr.vec))
}

