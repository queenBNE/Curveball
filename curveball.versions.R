curveball.randomise = function(I, N, directed=T, loops=F, return.all=F){
  if(directed)
    return(curveball.randomise.dir(I, N, loops, return.all))
  else if(loops){
    warning("The current implementation of the Curveball algorithm does not support randomization of undirected networks with self-loops")
    return(NULL)
  }
  return(curveball.randomise.undir(I, N, return.all))
}

curveball.randomise.dir =  function(I, N, loops, return.all){
  all = list()
  n=length(I) 
  for(run in 1:N){
    ij=sample(1:n,2)  
    i = ij[1]
    j = ij[2]    
    Ai=I[[i]]
    Aj=I[[j]] 
    AiAj.inter=intersect(Ai,Aj)   
    Ai_j = setdiff(Ai, AiAj.inter)
    if(!loops)
      Ai_j = setdiff(Ai_j, c(j)) 
    Aj_i = setdiff(Aj, AiAj.inter)
    if(!loops)
      Aj_i = setdiff(Aj_i, c(i))   
    if (length(Ai_j) > 0 & length(Aj_i) > 0){
        tot=union(Ai_j, Aj_i)
        tot.l = length(tot)
        tot=sample(tot, tot.l)
        L = length(Ai_j)
        I[[i]] = c(setdiff(Ai, Ai_j), tot[1:L])
        I[[j]] = c(setdiff(Aj, Aj_i), tot[(L+1):tot.l])
      }
    if(return.all)
      all[[run]] = I
    }
    if(return.all)
      return(all)
    return(I)
}

curveball.randomise.undir = function(I, N, return.all){
  all = list()
  n=length(I) 
  for(run in 1:N){
    ij=sample(1:n,2)  
    i = ij[1]
    j = ij[2]    
    Ai=I[[i]]
    Aj=I[[j]] 
    AiAj.inter=intersect(Ai,Aj)   
    Ai_j = setdiff(Ai, c(AiAj.inter, j))
    Aj_i = setdiff(Aj, c(AiAj.inter, i))
    
    if (length(Ai_j) > 0 & length(Aj_i) > 0){
      tot=union(Ai_j, Aj_i)
      tot.l = length(tot)
      tot=sample(tot, tot.l)
      L = length(Ai_j)
      I[[i]] = c(setdiff(Ai, Ai_j), tot[1:L])
      for(k in setdiff(I[[i]], Ai)){
        idx = which(I[[k]] == j)
        I[[k]] = union(i, I[[k]][-idx])
      }
      I[[j]] = c(setdiff(Aj, Aj_i), tot[(L+1):tot.l])
      for(l in setdiff(I[[j]],Aj)){
        idx = which(I[[l]] == i)
        I[[l]] = union(j, I[[l]][-idx])
      }
    }
    if(return.all)
      all[[run]] = I
  }
  if(return.all)
    return(all)
  return(I)
}



