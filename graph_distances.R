library("Matrix")
#library("MASS")

## Forest metrics
forest.metrics <- function(L, alpha){

    n = nrow(L)
    Q = solve(eye(n) + alpha*L)
    H = kappa(Q)

    return(H)
}

## Exponential distance
exp.dist <- function(P){
    w <- stationary.dist(P)
    Q <- outer(seq(1,1,length.out = nrow(P)), w)
    
    H = expm(P)*(1/Q)
    D = kappa(H)
}

## Expected commute time
ect <- function(P){
  
    n <- nrow(P)
    w <- power.method(t(P), tol = 1e-9)
    Q <- outer(seq(1,1,length.out=n), w)

    Z <- solve(diag(n) - P + Q)
    H <- kappa(Z * (1/Q))
    D <- sqrt(H)
    
    return(D)
}

diffusion.distance <- function(P, t, directed = FALSE){
  
  n <- nrow(P)
  w <- stationary.dist(P)
  Q <- outer(seq(1,1,length.out = n), w)

  if(directed == FALSE){
      P.2t <- mtx.exp(P,2*t)
      D.squared <- kappa(P.2t * (1/Q))
  }

  if(directed == TRUE){
      P.t <- mtx.exp(P,t)
      D.squared <- kappa(P.t %*% diag(1/w) %*% t(P.t))
  }
  
  D <- sqrt(D.squared)

  return(D)
}

laplacian.map <- function(W,k){

    L <- laplacian(W)
    decomp <- eigen(L)
    eigen.vals <- decomp$values[(nrow(L)-k):(nrow(L)-1)]
    eigen.vects <- decomp$vectors[,(nrow(L)-k):(nrow(L)-1)]
    eigen.vals.transformed <- 1/sqrt(eigen.vals)
    Psi <- eigen.vects * outer(seq(1,1,length.out = nrow(L)),
                               eigen.vals.transformed)
    
    return(Psi)
}

diffusion.map <- function(P, t, k=2){

    n <- nrow(P)
    w <- stationary.dist(P)
    
    Q <- outer(seq(1,1,length.out = n), w)
    Q.sqrt <- sqrt(Q)
    
    A <- t(Q.sqrt) * P * (1/Q.sqrt)
    
    decomp <- eigen(A)
    eigen.vals <- decomp$values[2:(k+1)]
    eigen.vects <- decomp$vectors[,2:(k+1)]
    eigen.vals.powered <- eigen.vals^t

    Psi <- outer(1/sqrt(w), seq(1,1,length.out = k)) * eigen.vects *
      outer(seq(1,1,length.out = n), eigen.vals.powered)
    return(Psi)
}

ect.map <- function(P,k=2){
    n <- nrow(P)
    w <- power.method(t(P), tol = 1e-9)
    Q <- outer(seq(1,1,length.out = n), w)
    Q.sqrt <- sqrt(Q)
    
    A <- t(Q.sqrt) * P * (1/Q.sqrt)

    decomp <- eigen(A)
    eigen.vals <- decomp$values[2:(k+1)]
    eigen.vects <- decomp$vectors[,2:(k+1)]
    eigen.vals.transformed <- 1/sqrt(1 - eigen.vals)

    Psi <- outer(1/sqrt(w), seq(1,1,length.out = k)) * eigen.vects *
      outer(seq(1,1,length.out = n), eigen.vals.transformed)

    return(Psi)
}

ect.pow <- function(P,t){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    Q <- seq(1,1,length.out=n) %*% t(w)

    Z <- solve(diag(n) - P + Q)
    Z.powered = mtx.exp(Z,t)
    H <- kappa(Z.powered * (1/Q))
    D <- sqrt(H)
    
    return(D)
}

out.of.sample.ect <- function(P,x){
    
    n <- nrow(P)
    w <- power.method(t(P), tol = 1e-9)
    Q <- outer(seq(1,1,length.out = n), w)

    A <- ginv(P - Q)

    Z <- solve(diag(n) - P + Q)
    z.star <- x%*%Z%*%A
    k.vect <- z.star * 1/w
}

gaussian.similarity <- function(dat, epsilon){
    n <- nrow(dat)
    W <- diag(n)
    
    W <- kappa(dat %*% t(dat))
    W = exp(-W/epsilon)

    diag(W) <- 0

    return(W)
}

gaussian.similarity.directed <- function(dat, epsilon, sparsity = 0.5){
    
    W <- gaussian.similarity(dat,epsilon)
    n <- nrow(W)

    for(i in 1:n){
        t <- runif(n)
        W[i,t < sparsity] <- 0
    }

    return(W)
}
        
laplacian <- function(W){
    S <- apply(W,1,sum)
    L <- diag(S) - W
    return(L)
}

sparsify.similarity <- function(W,k){
    n <- nrow(W)

    for(i in 1:n){
        tmp <- sort(W[i,], decreasing = TRUE, index.return = TRUE)
        index = tmp$ix[c(1:k)]
        W[i,index] <- 1
        W[i,-index] <- sum(tmp$x[index])/(100*(n-k))
    }

    W <- (W + t(W))/2

    return(W)
}

harmonic.sparsify <- function(W,k){
    n <- nrow(W)

    for(i in 1:n){
        tmp <- sort(W[i,], decreasing = FALSE, index.return = TRUE)
        W[i,] <- 0
        index = tmp$ix[1:k]
        W[i,index] <- 1/c(1:k)
    }

    return(W)
}
        
transition.matrix <- function(W){
    S <- apply(W,1,sum)
    P <- W/S
    return(P)
}

# Return the stationary distribution of P
stationary.dist <- function(P){
    
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    return(w)
}

transition.matrix.star <- function(P){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)

    Q <- outer(seq(1,1,length.out = n), w)
    P.star <- t(1/Q) * t(P) * Q

    return(P.star)
}
    
augmented.similarity <- function(W){
    n <- nrow(W)
    s = seq(1,1,length.out = n)
    W.new <- rbind(c(0,s),cbind(1,W))
    return(W.new)
}

mtx.exp <- function(X,n){
    if(n != round(n)){
        n <- round(n)
        warning("rounding exponent `n` to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X
    while(n > 0){
        if(n %% 2)
            phi <- phi %*% pot
        n <- n %/% 2
        pot <- pot %*% pot
    }
    return(phi)
}
            

power.method <- function(M,tol=1e-6){
    n = nrow(M)
    x = runif(n)
    x.new = x
    exit.status = FALSE
    while(exit.status != TRUE){
        x.new = M %*% x
        x.new.max = max(abs(x.new))
        x.new = x.new/x.new.max
        if(sum((x - x.new)^2) < tol)
          exit.status = TRUE
        x = x.new
    }
    x.norm = sum(abs(x.new))
    as.vector(abs(x.new)/x.norm)
}


kappa <- function(X){
    n <- nrow(X)
    D <- outer(diag(X),seq(1,1,length.out = n)) - X -
               t(X) + outer(seq(1,1,length.out = n), diag(X))
    return(D)
}

## Returns the tau transform of a square matrix D
tau <- function(D){
    n <- nrow(D)
    P <- projection.matrix(n)
    X <- (-1/2)*P%*%D%*%P
    
    return(X)
}

projection.matrix <- function(n){
    e <- seq(1,1,length.out=n)
    J <- (e%*%t(e))/n
    P <- diag(n) - J
}

