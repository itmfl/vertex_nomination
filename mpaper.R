library(igraph)
parse.tvks <- function(filename, l = 2, m = 5, n = 12){

  xyz <- read.table(filename, header = FALSE)
  A1 <- matrix(0, nrow = ncol(xyz) - 2, ncol = ncol(xyz) - 2)
  A2 <- matrix(0, nrow = ncol(xyz) - 2, ncol = ncol(xyz) - 2)

  for(i in 1:floor(nrow(xyz)/4)){
    tmp <- xyz[i,]
    tmp.idx <- which(tmp[2:(ncol(xyz) - 1)] == 1)
    tmp.topic <- tmp[ncol(xyz)]
    if(tmp.topic == 1){
      A1[tmp.idx[1],tmp.idx[2]] <- A1[tmp.idx[1], tmp.idx[2]] + 1
      A1[tmp.idx[2],tmp.idx[1]] <- A1[tmp.idx[2], tmp.idx[1]] + 1
    }
    else{
      A2[tmp.idx[1],tmp.idx[2]] <- A2[tmp.idx[1], tmp.idx[2]] + 1
      A2[tmp.idx[2],tmp.idx[1]] <- A2[tmp.idx[2], tmp.idx[1]] + 1
    }
  }

  A <- (A1 > 0 | A2 > 0)
  A <- A + 0

  red.idx <- (A1 > A2)
  green.idx <- (A2 >= A1) & (A2 != 0)

  E <- matrix("none", nrow = ncol(xyz) - 2, ncol = ncol(xyz) - 2)
  
  E[red.idx] <- "red"
  E[green.idx] <- "green"

  v.colors <- c(seq("red", "red", length.out = l),
                seq("blue", "blue", length.out = m - l),
                seq("green", "green", length.out = n - m))
  
  return(list(adjacency = A, labels = E, v.colors = v.colors))
}


  
      
    
    
