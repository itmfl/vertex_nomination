max0 <- function(tmp){

  if(tmp >= 0)
    return(tmp)
  else
    return(0)
}
binomial.sum <- function(m,n){
  s <- 0
  for(i in 0:(n-1)){
    s <- s + choose(m+i,i)
  }
  return(s)
}
