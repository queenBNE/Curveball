# Perturbation scores 
perturbation.scores = function(M, N, FUN){
  swapped = list(M)
  d = dim(M)
  for(i in 2:N){
    swapped[[i]] = FUN(swapped[[i-1]])
  }  
  ps = sapply(swapped, function(A){return(perturbation.score(M,A))})
  return(ps)
}

perturbation.score = function(A, B){
  D=sum(A)
  idx = which(A == 1)
  same = length(which(B[idx] == 1))
  return(1 - same/D)
}