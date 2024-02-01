#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GP Utility Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GP Sampling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample from a Gaussian process
sample_gp = function(mx, Rx){
  # Get length of mx
  n = length(mx)
  #Perform spectral decomposition on the new covariance matrix
  ev = eigen(Rx, symmetric = TRUE)
  psd = max(which(ev$values>0)) #Get the positive eigenvalues 
  evec = ev$vectors[,1:psd] #Get corresponding eigenvectors
  
  #Spectral decomposition of the partition of the covariance matrix that is full rank 
  Sp=evec[,1:psd]%*%diag(sqrt(ev$values[1:psd]))%*%t(evec[,1:psd])
  
  #Generate the predicted value from this distribution
  u=rnorm(n,0,1)
  out = mx+Sp%*%u
  return(out[,1])
}


# Predictive Distribution from GP
predict_dist_gp = function(y1,m1,m2,R11,R22,R12){
  #Set dimension parameters
  n1 = length(m1)
  n2 = length(m2)
  
  #Get inverse of R11 -- choleksy decomposition
  R11_chol = chol(R11)
  R11_inv = chol2inv(R11_chol)
    
  #Get mean and covaraince -- recall we assume mean 0 of the data
  mp = m2 + t(R12)%*%R11_inv%*%(y1-m1)
  Rp = R22 - t(R12)%*%R11_inv%*%R12
    
  out = list(mp = mp, Rp = Rp)
  return(out)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MH Sampling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Metropolis Hastings Functions
#MH Step with log density
mh_logstep = function(pdc, pdn, current_val, new_val){
  #Acceptance/Rejection -- LOG SCALE
  r = pdn - pdc
  if(is.infinite(r) | is.nan(r)){
    theta0 = current_val
    accept = 0
  }else{
    alpha = min(0, r)
    u = runif(1, min = 0, max = 1)
    if(alpha > log(u)){
      theta0 = new_val
      accept = 1
    }else{
      theta0 = current_val
      accept = 0
    }
  }
  
  mh_list = list(theta0 = theta0, accept = accept)
  return(mh_list)
}


#Adaptive MH function
adapt_delta = function(accept_num, last_accept, steps, current_delta, goal_rate = 0.44){
  a = accept_num - last_accept 
  r = a/steps
  
  delta_new = (r/goal_rate)*current_delta
  if(delta_new == 0){delta_new = current_delta}
  values = list(delta = delta_new, last_accept = accept_num)
  return(values)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Normalize Inputs on a [0,1] scale
norm_x = function(design_x, design_p){
  input_range = matrix(0, nrow = ncol(design_x), ncol = 2) 
  for(i in 1:ncol(design_x)){
    #Get range and scale the data -- this is using the prediction grid and the training points
    rng = range(c(design_p[,i], design_x[,i]))
    
    #Rescale the designs for the training x and testing x
    design_x[,i] = (design_x[,i] - rng[1])/(rng[2] - rng[1])
    design_p[,i] = (design_p[,i] - rng[1])/(rng[2] - rng[1])
    
    #store the ranges 
    input_range[i,] = rng
  }
  
  out = list(dx = design_x, dp = design_p, input_rng = input_range)
  return(out)
}

#Normalize Outputs - centered and scaled to have mean 0 and variance 1 (yt = ytrain, yf = ytest (true function))
norm_y = function(yt, yf){
  mean_yt=mean(yt)
  sd_yt=sd(yt)
  yt=(yt-mean_yt)/sd_yt
  yf=(yf-mean_yt)/sd_yt
  
  out = list(yt_norm = yt, yf_norm = yf, mean_yt = mean_yt, sd_yt = sd_yt)
  return(out)
}

#Get the pairwise differences between x vectors
pairwise_diff = function(design_mat){
  l_diff = list()
  for(i in 1:ncol(design_mat)){
    ltemp = list(abs(outer(design_mat[,i],design_mat[,i],"-")))
    #names(ltemp) = paste0('m',i)
    l_diff[[i]] = ltemp
  }
  names(l_diff) = paste0('l',1:ncol(design_mat))
  return(l_diff)
}
