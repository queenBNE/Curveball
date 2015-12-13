curveball.randomise.dir.multi =  function(I, N, loops=F, return.all){
  all = list()
  n=length(I) 
  for(run in 1:N){
    ij=sample(1:n,2)  
    i = ij[1]
    j = ij[2]    
    
    if(loops){
      if (length(I[[i]]) > 0 & length(I[[j]]) > 0){
        tot = c(I[[i]], I[[j]])
        l = length(tot)
        sample(tot, l)
        L = length(I[[i]])
        I[[i]] = tot[1:L]
        I[[j]] = tot[(L+1):l]
      }
    }else{
      Ai = I[[i]]
      fix_i = c()
      idx_j = which(Ai==j)
      if(length(idx_j) > 0){
        fix_i = rep(j, length(idx_j))
        Ai = Ai[-idx_j]
      }     
      Aj = I[[j]]
      fix_j = c()
      idx_i = which(Aj==i)
      if(length(idx_i) > 0){
        fix_j = rep(i, length(idx_i))
        Ai = Ai[-idx_i]
      } 
      L = length(Ai)
      if(L > 0 & length(Aj) > 0){
        tot = c(Ai, Aj)
        l = length(tot)
        sample(tot, l)
        I[[i]] = c(fix_i, tot[1:L])
        I[[j]] = c(fix_j, tot[(L+1):l])
      }
    }
    if(return.all)
      all[[run]] = I
  }
  if(return.all)
    return(all)
  return(I)
}