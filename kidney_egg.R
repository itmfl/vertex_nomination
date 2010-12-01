library("igraph")

#l = 10 
#m = 20 
#n = 100
#g1 <- graph.adjacency(g.tmp$adjacency, mode="undirected")

kidney.egg <- function(n,m,l,p.vec,s.vec){

    p1 <- p.vec[1]
    s1 <- s.vec[1]
    p2 <- p.vec[2]
    s2 <- s.vec[2]

    runif.mat <- matrix(0,nrow = n, ncol = n)
    runif.mat[col(runif.mat) > row(runif.mat)] <- runif(n*(n-1)/2)
    runif.mat <- (runif.mat + t(runif.mat))

    edge.mat <- matrix(p1 + p2, nrow = n, ncol =n)
    edge.mat[1:m,1:m] <- (s1 + s2)

    red.mat <- matrix(p1, nrow = n, ncol = n)
    red.mat[1:m,1:m] <- s1

    adjacency.mat <- matrix(0,nrow = n, ncol = n)
    edge.idx <- (runif.mat < edge.mat)
    adjacency.mat[edge.idx] <- 1

    labels.mat <- matrix("none", nrow = n, ncol = n)
    red.idx <- (runif.mat <= red.mat)
    green.idx <- (runif.mat > red.mat & runif.mat <= edge.mat)

    labels.mat[red.idx] <- "red"
    labels.mat[green.idx] <- "green"

    vertex.colors <- c(seq("red","red",length.out=l),
                     seq("blue","blue",length.out=m-l),
                     seq("green","green",length.out=(n-m)))

    return(list(adjacency=adjacency.mat,
              labels=labels.mat,
              v.colors = vertex.colors))
}
    
kidney.egg.rdpg <- function(n,m,l, kidney.red, kidney.green,
                            egg.red, egg.green){
  
  egg2egg.red <- sum(egg.red*egg.red)
  egg2egg.green <- sum(egg.green*egg.green)
  egg2egg.edge <- egg2egg.red + egg2egg.green
  egg2egg.none <- 1 - egg2egg.edge

  kidney2kidney.red <- sum(kidney.red*kidney.red)
  kidney2kidney.green <- sum(kidney.green*kidney.green)
  kidney2kidney.edge <- kidney2kidney.red + kidney2kidney.green
  kidney2kidney.none <- 1 - kidney2kidney.edge

  kidney2egg.red <- sum(kidney.red*egg.red)
  kidney2egg.green <- sum(kidney.green*egg.green)
  kidney2egg.edge <- kidney2egg.red + kidney2egg.green
  kidney2egg.none <- 1 - kidney2egg.edge

  labels.mat <- matrix("", nrow = n, ncol = n)
  labels.generation.mat <- matrix(0,nrow = n, ncol = n)
  labels.generation.mat[col(labels.generation.mat) > row(labels.generation.mat)] <- runif(n*(n-1)/2)
  labels.generation.mat <- (labels.generation.mat + t(labels.generation.mat))

  egg2egg.mat <- labels.generation.mat[1:m,1:m]
  kidney2kidney.mat <- labels.generation.mat[(m+1):n,(m+1):n]
  kidney2egg.mat1 <- labels.generation.mat[c(1:m,(m+1):n)]

  egg2egg.adjacency <- matrix(0,nrow=m,ncol=m)
  edge.idx <- (egg2egg.mat <= egg2egg.edge)
  egg2egg.adjacency[edge.idx] <- 1

  egg2egg.labels <- matrix("none", nrow=m, ncol = m)
  red.idx <- (egg2egg.mat <= egg2egg.red)
  green.idx <- ((egg2egg.mat > egg2egg.red) &
                (egg2egg.mat <= egg2egg.edge))
  egg2egg.labels[red.idx] <- "red"
  egg2egg.labels[green.idx] <- "green"
  

  kidney2kidney.adjacency <- matrix(0, nrow = (n-m), ncol = (n-m))
  edge.idx <- (kidney2kidney.mat <= kidney2kidney.edge)
  kidney2kidney.adjacency[edge.idx] <- 1
  
  kidney2kidney.labels <- matrix("none", nrow=(n-m), ncol = (n-m))
  red.idx <- (kidney2kidney.mat <= kidney2kidney.red)
  kidney2kidney.labels[red.idx] <- "red"
  green.idx <- ((kidney2kidney.mat > kidney2kidney.red)
                & (kidney2kidney.mat <= kidney2kidney.edge))
  
  kidney2kidney.labels[green.idx] <- "green"
 
  kidney2egg.adjacency1 <- matrix(0, nrow = m, ncol=(n-m))
  edge.idx <- (kidney2egg.mat1 <= kidney2egg.edge)
  kidney2egg.adjacency1[edge.idx] <- 1

  kidney2egg.labels1 <- matrix("none", nrow = m, ncol =(n-m))
  red.idx <- (kidney2egg.mat1 <= kidney2egg.red)
  kidney2egg.labels1[red.idx] <- "red"
  green.idx <- ((kidney2egg.mat1 > kidney2egg.red)
                & (kidney2egg.mat1 <= kidney2egg.edge))
  
  kidney2egg.labels1[green.idx] <- "green"

  labels.mat1 <- cbind(egg2egg.labels, kidney2egg.labels1)
  labels.mat2 <- cbind(t(kidney2egg.labels1), kidney2kidney.labels)
  labels.mat <- rbind(labels.mat1, labels.mat2)

  adjacency.mat1 <- cbind(egg2egg.adjacency, kidney2egg.adjacency1)
  adjacency.mat2 <- cbind(t(kidney2egg.adjacency1), kidney2kidney.adjacency)

  adjacency.mat <- rbind(adjacency.mat1, adjacency.mat2)
  diag(adjacency.mat) <- 0
  diag(labels.mat) <- "none"

  vertex.colors <- c(seq("red","red",length.out=l),
                     seq("blue","blue",length.out=m-l),
                     seq("green","green",length.out=(n-m)))

  return(list(adjacency=adjacency.mat,
              labels=labels.mat,
              v.colors = vertex.colors))
}

kidney.egg.plot <- function(kidney.egg.g){

  g <- graph.adjacency(kidney.egg.g$adjacency,mode="undirected")
  v.colors <- kidney.egg.g$v.colors
  
  t = get.edgelist(g)
  e.color <- seq("","",length.out = nrow(t))
  edge.labels = kidney.egg.g$labels

  for(i in 1:nrow(t)){
    e.color[i] = edge.labels[t[i,1]+1,t[i,2]+1] 
  }
  plot.igraph(g, layout=layout.circle, vertex.size=5,
              vertex.label=NA, vertex.color = v.colors,
              edge.color = e.color)  
}
