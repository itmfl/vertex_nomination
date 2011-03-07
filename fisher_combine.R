## Compute the Camp-Paulson approximation for the binomial distribution
camp.paulson <- function(n,p,k){

  a <- 1/(9*n - 9*k)
  b <- 1/(9*k + 9)
  r <- (k+1)*q/((n-k)*p)
  c <- (1 -b)*(r^(1/3))
  mu <- 1 - a
  sigma <- sqrt(b*(r^(2/3)) + a)
  l <- (c - mu)/sigma

  pnorm(l)
}

