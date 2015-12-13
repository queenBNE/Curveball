# Functions for the Curveball algorithm 

# Run the curveball algorithm for a matrix A, repeating step b-d N times 
curveball.run <- function(A, N){
  R = dim(A)[1]
  C = dim(A)[2]
  hp = curveball.step.a(A)
  for (rep in 1:N){
    hp = curveball.step.bcd(hp)
  }
  A = curveball.step.f(hp, R, C)
  return(A)
}

# Run the curveball algorithm for incidence lists I, repeating step b-d N times 
curveball.run <- function(I, N){
  for (rep in 1:N){
    I = curveball.step.bcd(I)
  }
  return(I)
}

# Perfrom one trade in the curveball algorithm 
curveball.trade <- function(A){
  R = dim(A)[1]
  C = dim(A)[2]
  hp = curveball.step.a(A)
  hp = curveball.step.bcd(hp)
  A = curveball.step.f(hp, R, C)
  return(A)
}

curveball.step.a <- function(A){
  hp=list()
  for (row in 1:dim(A)[1]) {hp[[row]]=(which(A[row,]==1))}
  return(hp)
}

curveball.step.bcd <- function(hp){
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
    tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
    L=l_a-l_ab
    hp[[AB[1]]] = c(ab,tot[1:L])
    hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])
  }
  return(hp)
}

curveball.step.f <- function(hp, R, C){
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  return(rm)
}

curveball.getTransitionMatrix <- function(states){
  stateCount = length(states)
  P = mat.or.vec(nr=stateCount, nc=stateCount)
  for(i in 1:stateCount)
    for(j in 1:stateCount)
      if(i != j)
        P[i,j] = curveball.getTransitionProbability(states[[i]], states[[j]])
  diag(P) = 1 - colSums(P)
  return(P)
}

curveball.getTransitionProbability <- function(A, B){
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
  p = 2/(r*(r-1)) * (factorial(s_i) * factorial(s_j)) / factorial(s_i + s_j)
  return(p)
}  

curveball.getRuntime <- function(A, states, eps, steps){
  P = curveball.getTransitionMatrix(states)
  N = getMixingTime(P, eps)
  
  # Run specified number of steps in the curveball algorithm and measure duration
  adjacencyLists = curveball.step.a(A)
  t = system.time(
    for(i in 1:attempts){
      adjacencyLists = curveball.step.bcd(adjacencyLists)
    })  
  return(c(N,t[1]))
}
