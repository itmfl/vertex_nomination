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

procrustes <- function(A.mat, B.mat){

    AB.mat <- t(A.mat) %*% B.mat
    AB.mat.svd <- svd(AB.mat)
    T <- t(AB.mat.svd$v)%*%AB.mat.svd$u
    
    return(B.mat %*% T)
}

count.red.nn <- function(dist.mat, colors, k){

    idx <- which((colors == "blue") | (colors == "green"))
    num.red.nn <-seq(0,0,length.out = length(idx))

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

inverse.rdpg.nominate <- function(g, dim = 2, do.procrustes = FALSE){

    A.list <- adjacency.list(g$adjacency, g$labels)

    x.red <- svd.extract(A.list$red, k = dim)
    x.green <- svd.extract(A.list$green, k = dim)
    x.separate <- cbind(x.red,x.green)
    
    if( do.procrustes == FALSE){
        d.mat <- as.matrix(dist(x.separate))
    }
    else {
        x.common <- svd.extract(g$adjacency, k = 2*dim)
        x.separate.procrustes <- procrustes(x.common,x.separate)
        d.mat <- as.matrix(dist(x.seperate.procrustes))
    }
    
    xyz <- count.red.nn(d.mat, g$v.colors,
                        k = sum(g$v.colors == "red"))

    return(list(order=xyz$order, value=xyz$vector))
}

g.statistics <- function(g){

    red.idx <- which(g$v.colors == "red")
    t1.vec <- apply(g$adjacency[,red.idx],1,sum)

    A.list <- adjacency.list(g$adjacency, g$labels)
    t2.vec <- apply(A.list$red, 1, sum)
    
    return(list(t1 = t1.vec, t2 = t2.vec))
}

tx.statistics.nominate <- function(g){

    g.stat <- g.statistics(g)
    blue.green.idx <- which(g$v.colors != "red") 

    tx <- g.stat$t1[blue.green.idx] 
    tx.sorted <- sort(tx, decreasing = TRUE, index.return = TRUE)
    
    return(list(value = tx.sorted$x, order = blue.green.idx[tx.sorted$ix]))
}

tn.statistics.nominate <- function(g){
    
    g.stat <- g.statistics(g)
    blue.green.idx <- which(g$v.colors != "red") 
    
    tn <- g.stat$t2[blue.green.idx]
    tn.sorted <- sort(tn, decreasing = TRUE, index.return = TRUE)
    
    return(list(value = tn.sorted$x, order = blue.green.idx[tn.sorted$ix]))
}
               
tf.statistics.nominate <- function(g){
    
    g.stat <- g.statistics(g)
    blue.green.idx <- which(g$v.colors != "red") 
    
    tf <- g.stat$t1[blue.green.idx] + g.stat$t2[blue.green.idx]
    tf.sorted <- sort(tf, decreasing = TRUE, index.return = TRUE)
    
    return(list(value = tf.sorted$x, order = blue.green.idx[tf.sorted$ix]))
}

graph.distance.nominate <- function(g){

    P <- transition.matrix(g$adjacency)
    D <- diffusion.distance(P, t = 5)

    red.idx <- which(g$v.colors == "red")
    blue.green.idx <- which(g$v.colors != "red")
    l <- length(red.idx)
    
    geometric.dist <- exp(1/l*apply(log(D[blue.green.idx,red.idx]),1,mean))

    geometric.dist.sorted <- sort(geometric.dist, decreasing = FALSE,
                                  index.return = TRUE)
    
    return(list(order = blue.green.idx[geometric.dist.sorted$ix],
                value = geometric.dist.sorted$x))
}

nominate.vertex <- function(g, method){
    
    if(method == "inverse.rdpg")
      return(inverse.rdpg.nominate(g))
    if(method == "diffusion.distance")
      return(graph.distance.nominate(g))
    if(method == "tx.statistics")
      return(tx.statistics.nominate(g))
    if(method == "tn.statistics")
      return(tn.statistics.nominate(g))
    if(method == "tf.statistics")
      return(tf.statistics.nominate(g))
}



