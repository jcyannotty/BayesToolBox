#-----------------------------------------------------
# Samplers
#-----------------------------------------------------
# Gumbel sampling for categorical variables
rgumbel = function(n){
  u = runif(n,0,1)  
  rgum = -log(-log(u))
  return(rgum)
}

# Sample categorical variable with softmax pmf
rsoftmax = function(k,pvec,c=0,logprob = TRUE){
  # Draw gumbels 
  g = rgumbel(k)
  
  # define gamma
  if(logprob){
    # pvec = log unnormalized probability (already removes exp(..) with the log)
    gam = pvec + c
  }else{
    # pvec is unnormalized prob (usually includes the exponetial from the numerator) 
    gam = log(pvec) + c  
  }
  
  # Draw categories
  out = which.max(gam+g)
  return(out)
}

# Multivariate normal
rmvn = function(mvec, V){
  # Spectral decomposition to get V^1/2
  E = eigen(V)
  Ev = E$vectors
  Eval = E$values
  S = Ev%*%diag(sqrt(Eval))
  
  # Generate std normal
  N = nrow(V)
  z = rnorm(N, 0, 1)
  
  # Get the MVN random variable
  out = mvec + S%*%z
  return(out)
}

# Log metropolis hastings
# nprop: proposal of c -> n
# cprop: proposal of n -> c
logmh = function(log_nden, log_nprop, ntheta, log_cprop, log_cden, ctheta){
  # Acceptance ratio
  a = log_nden - log_cden + log_cprop - log_nprop
  a = min(0,a)
  u = runif(0,1)
  
  # Accept/reject
  if(log(u)< a){
    accept = TRUE
    theta = ctheta
  }else{
    accept = FALSE
    theta = ntheta
  }
  out = list(theta = theta, accept = accept)
  return(out)
}

# Grid 2-D Design Sampler 
# nc = number columns, nr = number rows, n = total points
# xmin, xmax = 2-d vector of min/max values for each dimension
grid_2d_design = function(n1,n2, xmin = c(-1,-1), xmax = c(1,1)){
  # Generate n uniform rvs
  n = n1*n2
  ux = runif(n,0,1)
  uy = runif(n,0,1)
  
  # Dimensions for each rectangle
  x1_len = (xmax[1] - xmin[1])/n1
  x2_len = (xmax[2] - xmin[2])/n2
  xgrid = expand.grid(1:n1,1:n2)
  
  # Get points
  x1 = ux*x1_len + x1_len*(xgrid[,1]-1) + xmin[1]
  x2 = uy*x2_len + x2_len*(xgrid[,2]-1) + xmin[2]
  
  # Join data
  xdata = cbind(x1,x2)
  return(xdata)
}

# Generate test grid around a given set of points
# d = length of square side, nbd = number of points in the neighborhood
grid_2d_testset = function(x_train, d = 1, nbd = 5){
  x_test = matrix(0, ncol = 2, nrow = nrow(x_train)*nbd)
  for(i in 1:nrow(x_train)){
    dx = d*runif(nbd,-0.5,0.5)
    dy = d*runif(nbd,-0.5,0.5)
    start = (i-1)*nbd + 1
    end = i*nbd
    x_test[start:end, ] = cbind(x_train[i,1]+dx,x_train[i,2]+dy)    
  }
  return(x_test)
}


# MOVE TO EXAMPLE FILE
# MALA Log Proposal
#mala_logprop = function(){}

# Test
# rgumbel(5)
# smax = 0
# nd = 50000
# pvec = c(0.1,0.2,0.3,0.2,0.2)
# for(i in 1:nd){
#   smax[i] = rsoftmax(5,pvec)
# }
# table(smax)/nd
# 
# 
# smax = 0
# pvec = exp(-c(4,2,1,2.3,0.4))
# const = sum(pvec)
# for(i in 1:nd){
#   smax[i] = rsoftmax(5,pvec,c = 0)
# }
# table(smax)/nd
# pvec/const

