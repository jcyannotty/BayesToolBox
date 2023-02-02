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

