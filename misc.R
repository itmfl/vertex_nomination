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

test4 <- function(mc){

  GT = seq(0,0,length.out = mc)
  LT = seq(0,0,length.out = mc)
  EQ = seq(0,0,length.out = mc)
  for( i in 1:mc){
    a1 = runif(10)
    a2 = runif(10)
    if(max(a1) > max(a2))
      GT[i] = 1
    if(max(a2) > max(a1))
      LT[i] = 1
    if(max(a2) == max(a1))
      EQ[i] = 1
  }
  return(list(a=sum(GT),b=sum(EQ),c=sum(LT)))
}

test5 <- function(mc){
  GTa = seq(0,0,length.out = mc)
  LTa = seq(0,0,length.out = mc)
  EQa = seq(0,0,length.out = mc)

  GTb = seq(0,0,length.out = mc)
  LTb = seq(0,0,length.out = mc)
  EQb = seq(0,0,length.out = mc)
  
  GTc = seq(0,0,length.out = mc)
  LTc = seq(0,0,length.out = mc)
  EQc = seq(0,0,length.out = mc)
  
  for( i in 1:mc){
    a1 = rbinom(1,20,0.2)
    b1 = rbinom(1,10,0.3)
    a2 = rbinom(1,25,0.2)
    b2 = rbinom(1,15,0.3)
    c1 = a1 + b1 
    c2 = a2 + b2 

    GTa[i] <- (a1 > a2)
    LTa[i] <- (a1 < a2)
    EQa[i] <- (a1 == a2)
    
    GTb[i] <- (b1 > b2)
    LTb[i] <- (b1 < b2)
    EQb[i] <- (b1 == b2)

    GTc[i] <- (c1 > c2)
    LTc[i] <- (c1 < c2)
    EQc[i] <- (c1 == c2)
  }

  return(list(gta = sum(GTa), lta = sum(LTa), eqa = sum(EQa),
              gtb = sum(GTb), ltb = sum(LTb), eqb = sum(EQb),
              gtc = sum(GTc), ltc = sum(LTc), eqc = sum(EQc)))
}

pvpm <- function(a,b){

  s1 <- sqrt(b) + sqrt(a/b*(b-a)) - a/sqrt(b)
  s2 <- a/sqrt(b) + sqrt(a/b*(b-a)) - sqrt(b)

  lambda <- s2^2/(s1^2 + s2^2)

  y1 <- sqrt((1 - lambda)*b)
  y2 <- sqrt(lambda*b)

  x1 <- sqrt(1 - lambda)*a/sqrt(b) - sqrt(lambda)*sqrt(a/b*(b-a))
  x2 <- sqrt(lambda)*a/sqrt(b) + sqrt(1 - lambda)*sqrt(a/b*(b-a))

  return(list(x = c(x1,x2), y = c(y1,y2)))
}

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
  

