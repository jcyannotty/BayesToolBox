#------------------------------------------------
# Optimization helper functions
#------------------------------------------------
# Min: ||a - b||^2 where a is unconstrained & b is constrained to simplex
simplex_l2 = function(a,s=1){
  # Sort from max to min
  aorder = rev(order(a))
  astar = a[aorder]
  n = length(a)

  # Find the largest k such that a[k] - lam[k] > 0  
  lam = (cumsum(astar) - s)/(1:n)
  h = max(which((astar-lam)>0))
  
  # Return simplex projection
  b = sapply(astar,function(x) max(x-lam[h],0))
  b = b[order(aorder)]
  return(b)
}
