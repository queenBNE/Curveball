# Functions for the Good-Shuffle algorithm 

# Run the goodshuffle algorithm for a matrix A, repeating step b-d N times 
goodshuffle.run <- function(A, N){
  R = dim(A)[1]
  C = dim(A)[2]
  hp = goodshuffle.step.a(A)
  for (rep in 1:steps){
    hp = goodshuffle.step.bcd(hp)
  }
  A = goodshuffle.step.f(hp, R, C)
  return(A)
}

# Perfrom one trade in the goodshuffle algorithm 
goodshuffle.trade <- function(A){
  R = dim(A)[1]
  C = dim(A)[2]
  hp = goodshuffle.step.a(A)
  hp = goodshuffle.step.bcd(hp)
  A = goodshuffle.step.f(hp, R, C)
  return(A)
}

goodshuffle.step.a <- function(A){
  hp=list()
  for (row in 1:dim(A)[1]) {hp[[row]]=(which(A[row,]==1))}
  return(hp)
}

goodshuffle.step.bcd <- function(hp){
  l_hp=length(hp)
  AB=sample(1:l_hp,2)  
  a=hp[[AB[1]]]
  b=hp[[AB[2]]]
  ab=intersect(a,b)
  l_ab=length(ab)
  l_a=length(a)
  l_b=length(b)
  if ((l_ab %in% c(l_a,l_b))==F){
    tot=setdiff(c(a,b),ab)
    l_tot=length(tot)
    L=l_a-l_ab
    new_a = a
    new_b = b
    while(setequal(new_a, a)){
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      new_a = c(ab,tot[1:L])
      new_b = c(ab,tot[(L+1):l_tot])
    }
    hp[[AB[1]]] = new_a
    hp[[AB[2]]] = new_b
  }
  return(hp)
}

goodshuffle.step.f <- function(hp, R, C){
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  return(rm)
}

goodshuffle.getTransitionMatrix <- function(states){
    stateCount = length(states)
    P = mat.or.vec(nr=stateCount, nc=stateCount)
    for(i in 1:stateCount)
      for(j in 1:stateCount)
        if(i != j)
          P[i,j] = goodshuffle.getTransitionProbability(states[[i]], states[[j]])
    diag(P) = 1 - colSums(P)
    return(P)
}

goodshuffle.getTransitionProbability <- function(A, B){
  r = dim(A)[1]
  # Find out which rows are different 
  t = A == B
  rowIdx = which(apply(t, 1, prod) == 0)
  # If more than two the transition probability equals 0 
  if(length(rowIdx) > 2)
    return(0)
  # If the same matrix we calculate the transition probability afterwards
  if(length(rowIdx) == 0)
    return(0)
  
  # The row indices of the two rows that are different
  i = rowIdx[1]
  j = rowIdx[2]
  s_i = length(which(A[i,] - A[j,] == 1))
  s_j = length(which(A[j,] - A[i,] == 1))
  p = 2/(r*(r-1)) * (factorial(s_i) * factorial(s_j)) / (factorial(s_i + s_j) - (factorial(s_i) * factorial(s_j)))
  return(p)
}

goodshuffle.getRuntime <- function(A, states, eps, steps){
  P = goodshuffle.getTransitionMatrix(states)
  N = getMixingTime(P, eps)
  
  # Run specified number of steps in the curveball algorithm and measure duration
  adjacencyLists = goodshuffle.step.a(A)
  t = system.time(
    for(i in 1:attempts){
      adjacencyLists = goodshuffle.step.bcd(adjacencyLists)
    })  
  return(c(N,t[1]))
}
