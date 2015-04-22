# Binary matrix randomisation functions

getAllBinaryMatricesWithRowAndColSums = function(M, max=100, prnt = F){
  matrices = list()
  # Set the initial matrix as the current state
  k <- 1
  matrices[[k]] <- M
  nextNew <- c(1)
  
  # Simple breadth-first algorithm - stop when no new matrices are found
  while(length(nextNew) > 0 && k < max){
    newMatrices <- nextNew 
    nextNew <- c()
    for(j in newMatrices){
      if(prnt)
        print(paste("Finding neigbors of state", j))
      checkerboards <- countCheckerboards(matrices[[j]], T)[[2]]
      if(length(checkerboards) == 4)
        checkerboards = matrix(checkerboards, nr=1, nc=4)
      for(c in 1:dim(checkerboards)[1]){
        chkr = checkerboards[c,]
        # Swap the edges to create a new network
        newState <- matrices[[j]]
        newState[chkr[1], chkr[2]] <- newState[chkr[1], chkr[2]] - 1
        newState[chkr[3], chkr[4]] <- newState[chkr[3], chkr[4]] - 1
        newState[chkr[1], chkr[4]] <- newState[chkr[1], chkr[4]] + 1
        newState[chkr[3], chkr[2]] <- newState[chkr[3], chkr[2]] + 1
        
        # Weights, todo
        w1=w2=0
        
        # Check if this network already existed in the stategraph
        isOld <- 0
        for(l in 1:length(matrices)){
          state <- matrices[[l]]
          areEqual <- prod(state == newState)
          isOld <- isOld + areEqual 
        }
        # If this is a new state, save it to the matrices and create an edge from current
        # state to the new state in the stategraph
        if(!isOld){
          k <- k + 1
          if(prnt)
            print(paste("New state", k))
          matrices[[k]] <- newState
          nextNew <- c(nextNew, k)
        }
      }
    }
  }
  return(matrices)
}

countCheckerboards <- function(A, saveList=F){
  count = 0
  idx = which(A == 1, arr.ind = T)
  l = dim(idx)[1]
  max=l*(l-1)/2
  checkerboards = mat.or.vec(nr=max, nc=4)
  for (i in 1:(l-1))
    for (j in (i+1):l){
      ri = idx[i,1]
      ci = idx[i,2]
      rj = idx[j,1]
      cj = idx[j,2]
      if(A[ri,cj] == 0 && A[rj, ci] == 0){
        count = count + 1
        if(saveList)
          checkerboards[count,] = c(ri,ci,rj,cj)
      }
    }
  if(saveList)
    return(list(count, checkerboards[1:count,]))
  else 
    return(list(count, NULL))
}


# Functions needed by speed comparision script 

# Generate a matrix Erdos-Renyi style 
generateRandomBinaryMatrix = function(r,c,m=0,p=0){
  result = mat.or.vec(nr=r, nc=c)
  if(m != 0)
    idx = sample(1:(r*c), m)
  else if(p != 0){
    samp = runif(r*c)
    idx = which(samp <= p)
  } 
  result[idx] = 1
  return(result)
}
