#--------------------------------------------------------------------
# Kernel Functions
# Notes:
# Kernels are in terms of the euclidean distance bw x1 and x2: sqrt(sum(x1-x2)^2)
# The product kerenel is the same as the single kerenl in terms of euclidena distance for sqr exponential
# This is not the case for the others
#--------------------------------------------------------------------
# Squared Exponential/Gaussian
# length scale = lambda, scale = sig
sqr_exp_kernel = function(x1,x2,sig,lam){
  # distance between points taken relative to euclidean
  if(length(lam) == 1){
    d = sum((x1-x2)^2)
    k = (sig^2)*exp(-d/(2*lam))
  }else{
    d = sum((x1-x2)^2/lam)
    k = (sig^2)*exp(-d/2)
  }
  return(k)
}


matern_kernel = function(x1,x2,sig,lam,nu){
  if(length(nu) == 1 & length(lam) == 1){
    # Kernel in terms of euclidean distance
    d = sqrt(sum(x1-x2)^2)
    a = (2)^(1-nu)/gamma(nu) 
    b = (sqrt(2*nu)*d/lam)^nu 
    if(d>0){k = (sig^2)*a*b*besselK(sqrt(2*nu)*d/lam,nu)}else(k = sig^2)
    
  }else{
    # Product kernel for matern 
    p = length(x1)
    if(length(lam) != p){lam = rep(lam[1],p)}
    if(length(nu) != p){lam = rep(nu[1],p)}
    d = abs(x1-x2)
    a = (2)^(1-nu)/gamma(nu) 
    b = (sqrt(2*nu)*d/lam)^nu 
    if(d>0){k = (sig^2)*exp(sum(log(a*b*besselK(sqrt(2*nu)*d/lam,nu))))}else(k = sig^2)
  }
  return(k)
}


power_exp_kernel = function(x1,x2,sig,lam,gam){
  if(length(gam) == 1 & length(lam) == 1){
    # Kernel in terms of euclidean distance
    d = sqrt(sum(x1-x2)^2)/lam
    k = (sig^2)*exp(-d^gam)
  }else{
    # Product kernel for matern 
    p = length(x1)
    if(length(lam) != p){lam = rep(lam[1],p)}
    if(length(gam) != p){gam = rep(gam[1],p)}
    
    # Get the component distance
    d = abs(x1-x2)/lam
    k = (sig^2)*exp(-sum(d^gam))
  }
  return(k)
}


rational_quad_kernel = function(x1,x2,sig,lam,alpha){
  if(length(lam) == 1 & length(alpha) == 1){
    # Kernel in terms of euclidean distance
    d = sqrt(sum(x1-x2)^2)
    k = (sig^2)*(1 + (d^2)/(2*alpha*lam^2))^(-alpha)
  }else{
    p = length(x1)
    if(length(lam) != p){lam = rep(lam[1],p)}
    if(length(alpha) != p){lam = rep(alpha[1],p)}
    
    # Get the component distance and kernel 
    d = abs(x1- x2)
    kc = (1 + (d^2)/(2*alpha*lam^2))^(-alpha)
    k = (sig^2)*exp(sum(log(kc)))
  }
  return(k)
}


wendland_kernel = function(x1,x2,sig,lam,D,q){
  d = sqrt(sum(x1-x2)^2)/lam
  j = floor(D/2) + q + 1
  if(q == 0){
    k = max((1-d),0)^j
  }
  if(q == 1){
    k = max((1-d),0)^(j+1)*((j+1)*d+1)
  }
  if(q == 2){
    k = max((1-d),0)^(j+2)*((j^2+4*j+3)*d^2+(3*j+6)*d + 3)/3
  }
  return(k)
}

# Laplacian kernel - based on L1 distance
laplace_kernel = function(x1,x2,sig,lam){
  if(length(lam) == 1){
    # Kernel in terms of euclidean distance
    d = sum(abs(x1-x2))/lam
    k = (sig^2)*exp(-d)
  }else{
    # Product kernel for matern 
    p = length(x1)
    if(length(lam) != p){lam = rep(lam[1],p)}
    if(length(gam) != p){lam = rep(gam[1],p)}
    
    # Get the component distance
    d = abs(x1-x2)/lam
    k = (sig^2)*exp(-sum(d))
  }
  return(k)
}

