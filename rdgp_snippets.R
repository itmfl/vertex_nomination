pvpm2 <- function(p,s,lambda){

  x <- list(length(p))
  y <- list(length(p))

  for(i in 1:length(p)){
    a <- p[i]
    b <- s[i]

    if(a == b){
      x[[i]] <- sqrt(a)
      y[[i]] <- sqrt(b)
    }
    else{
    
      y1 <- sqrt((1 - lambda)*b)
      y2 <- sqrt(lambda*b)

      x1 <- sqrt(1 - lambda)*a/sqrt(b) - sqrt(lambda)*sqrt(a/b*(b-a))
      x2 <- sqrt(lambda)*a/sqrt(b) + sqrt(1 - lambda)*sqrt(a/b*(b-a))

      x[[i]] <- c(x1,x2)
      y[[i]] <- c(y1,y2)
    }
  }

  return(list(x = x, y = y))
  
}
