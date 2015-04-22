# Functions to calculate the mixing time for a MC
getMixingTime = function(P, eps){
  stateCount = dim(P)[1]
  N = mat.or.vec(nr=1, nc=stateCount)
  for(i in 1 : min(stateCount,100)){
    e_i = rep(0, stateCount)
    e_i[i] = 1
    N[i] = getMixingTimeStartingAt(P, eps, e_i)
  }
  return(max(N))
}

# Find the mixing time starting from initial distribution pi_0
getMixingTimeStartingAt<- function(P, eps, pi_0){
  pi_n = pi_0
  pi = rep(1/length(pi_0), length(pi_0))
  n = 0
  while(totalVariationDistance(pi, pi_n) >= eps){
    n = n + 1
    pi_n = P %*% pi_n
  }
  return(n)
}

# A distance measure for discrete probability distribution
totalVariationDistance <- function(p1, p2){
  d = 0.5*sum(abs(p1-p2))
  return(d)
}

# Matrix powers
pow <- function(A, n){
  if(n == 1)
    return(A)
  if(n == 2)
    return(A%*%A)
  if(n > 2)
    return(A%*%pow(A,n-1))
}


