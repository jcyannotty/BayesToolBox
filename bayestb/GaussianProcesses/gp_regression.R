#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GP Regression with custom kernel function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log prior for covariance parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_prior = function(theta,prtype,hyper){
  if(prtype == "invgamma"){
    alpha = hyper[1]
    beta = hyper[2]
    den = log(beta^alpha/gamma(alpha)) - (alpha + 1)*log(theta) - beta/theta  
  }else if(prtype == "beta"){
    den = log(dgamma(theta, shape = hyper[1], rate = hyper[2]))
  }else if(prtype == "uniform"){
    den = log(dunif(theta, a = hyper[1], rate = b[2]))
  }else{
    stop("Choose between: invgamma, beta, or uniform")
  }
  return(den)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Log Posterior: 
# Requires cov_function wrapper: cov_function(x1,x2,params) 
log_post = function(y, x, cov_function, pvec, hplist, priornames, nug = 1e-6){
  #Check conditions
  ny = length(y)
  for(j in 1:length(pvec)){
    if(priornames[j] == "beta"){
      if(pvec[j]<0 || pvec[j]>1){return(-Inf)}  
    }else if(priornames[j] == "invgamma"){
      if(pvec[j]<0){return(-Inf)}
    }else if(priornames[j] == "gamma"){
      if(pvec[j]<0){return(-Inf)}
    }else if(priornames[j] == "uniform"){
      if(pvec[j]<hplist[[j]][1] || pvec[j]>hplist[[j]][2]){return(-Inf)}
    }
  }

  # Get x pairs
  xpairs = expand.grid(1:ny,1:ny)
  
  #Get covaraince matrices for the GP modeling the function and the random error
  Rx = cov_function(x1 = x[xpairs[,1],], x2 = x[xpairs[,2],], params = pvec,
                    nrow = ny ,ncol = ny,sig2 = FALSE)
  #Rx = matrix(Rx, nrow = ny, ncol = ny, byrow = TRUE)

  #Combine these matrices to get one covaraince matrix
  V = Rx + diag(nug,ny)  
  
  #Get V inverse via Choleski Decomposition and log determinant
  #Check to see if all eigen values are positive to avoid an error
  if(any(eigen(V)$values < 0)){return(-Inf)}
  Vchol = chol(V)
  Vinv = chol2inv(Vchol)
  logdetV = 2*sum(log(diag(Vchol)))
  
  #Get Log-likelihood (ignoring constant of 2*pi)
  loglike = -0.5*logdetV - 0.5*t(y)%*%Vinv%*%y
  
  #Get the log priors for lambda 
  logpr = 0
  for(j in 1:length(pvec)){logpr = logpr + log_prior(pvec[j],priornames[j],hplist[[j]])}

  #Compute log posterior
  logpost = loglike + logpr
  return(logpost)
  }
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fitting Gaussian Process model
gp_train = function(y, x, cov_function, priornames,hp, mh_proposal, nd, nadapt, nburn, adaptevery = 100, nug = 1e-6){
  #Initialize posterior vectors
  N = nd + nadapt + nburn + 1
  ny = length(y)
  ntheta = length(priornames)
  postmat = matrix(0,nrow = N, ncol = ntheta)
  if(!is.matrix(x)){
    x = as.matrix(x, nrow = ny, ncol = 1)
  }
  
  #Initialize acceptance lists 
  accept_list = list()
  for(j in 1:ntheta){
    accept_list[[j]] = list(total = 0, last = 0)
    if(priornames[j] == "invgamma"){
      postmat[1,j] = hp[[j]][2]/(hp[[j]][1] - 1)
    }else if(priornames[j] == "beta"){
      postmat[1,j] = hp[[j]][1]/(hp[[j]][1] + hp[[j]][2])
    }else if(priornames[j] == "uniform"){
      postmat[1,j] = (hp[[j]][1] + hp[[j]][2])/2
    }
  } 

  #MCMC with Metropolis Hastings 
  adapt_iter = 0
  for(i in 2:N){
    # Sample each parameter
    for(j in 1:ntheta){
      cur_val = postmat[i-1,]
      if(j > 1){
        cur_val[1:(j-1)] = postmat[i,1:(j-1)]  
      }
      new_val = cur_val
      new_val[j] = runif(1, cur_val[j] - mh_proposal[[j]], cur_val[j] + mh_proposal[[j]])
      
      #Get densities 
      pdn = log_post(y, x, cov_function, new_val, hp, priornames, nug = nug)
      pdc = log_post(y, x, cov_function, cur_val, hp, priornames, nug = nug)
      
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
      postmat[i,j] = mh$theta0[j]
      if(i < nadapt + 1){
        accept_list[[j]]$total = accept_list[[j]]$total + mh$accept
      }
      
    }

    #Progress Tracker
    cat(i/N*100," percent complete\r")
    adapt_iter = adapt_iter + 1
  
    if(i < nadapt + 1 & adapt_iter == adaptevery){
      for(j in 1:ntheta){
        if((accept_list[[j]]$total-accept_list[[j]]$last)/adaptevery > 0.49 | 
           (accept_list[[j]]$total-accept_list[[j]]$last)/adaptevery < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_list[[j]]$total, accept_list[[j]]$last, adaptevery, mh_proposal[[j]])
          accept_list[[j]]$last = ar$last_accept
          mh_proposal[[j]] = ar$delta
        }else{
          accept_list[[j]]$last = accept_list[[j]]$total
        }  
      }
      
    }
  }
  
  cat("\n Complete.\n")
  
  # Get acceptance rates
  accept_rates = 0
  for(j in 1:ntheta){
    accept_rates[j] = accept_list[[j]]$total/nadapt
  }
  
  #Set the sample draws to save
  sdr = (nadapt+nburn+1):N
  out = list(y = y, x = x, mh_prop = mh_proposal, accept_rates = accept_rates,
             post = postmat[sdr,], priors = priornames)
  return(out)
}


# Predictions
gp_predict = function(fit, y, x_train, x_test, m_train, m_test, cov_function, nug = 1e-6, 
                      pred_y = FALSE){
  # Check if x is a matrix
  ny = length(y) #Number of observations
  if(!is.matrix(x_train)){
    x_train = as.matrix(x_train, nrow = ny, ncol = 1)
  }
  
  if(!is.matrix(x_test)){
    np = length(x_test)
    x_test = as.matrix(x_test, nrow = np, ncol = 1)
  }
  
  #Set dimension parameters
  p = ncol(x_train) #Number of predictors
  np = nrow(x_test) #Number of predictions
  Npost = nrow(fit$post) #Number of posterior draws
  nn = np + ny
  
  # Get x pairs
  x11_pairs = expand.grid(1:ny,1:ny)
  x12_pairs = expand.grid(1:ny,1:np)
  x22_pairs = expand.grid(1:np,1:np)
  
  fit_pred = matrix(NA, nrow = Npost, ncol = np)
  
  #Make predictions based on the Npost posterior draws
  for(i in 1:Npost){
    #Get covaraince matrices for the GP modeling the function and the random error
    pvec = fit$post[i,]
    R11 = cov_function(x1 = x_train[x11_pairs[,1],], x2 = x_train[x11_pairs[,2],], params = pvec,
                       nrow = ny, ncol = ny, sig2 = TRUE)
    R12 = cov_function(x1 = x_train[x12_pairs[,1],], x2 = x_test[x12_pairs[,2],], params = pvec, 
                       nrow = ny, ncol = np, sig2 = pred_y)
    R22 = cov_function(x1 = x_test[x22_pairs[,1],], x2 = x_test[x22_pairs[,2],], params = pvec, 
                       nrow = np, ncol = np,sig2 = pred_y)

    # Reshape
    #R11 = matrix(R11, nrow = ny, ncol = ny, byrow = TRUE)
    #R12 = matrix(R12, nrow = ny, ncol = np, byrow = FALSE)
    #R22 = matrix(R22, nrow = np, ncol = np, byrow = TRUE)
    
    # Add the nugget for stability
    R11 = R11 + diag(nug,ny)

    # Get predictive distribution
    out = predict_dist_gp(y,m_train,m_test,R11,R22,R12)    
    
    #Perform spectral decomposition on the new covariance matrix
    ev = eigen(out$Rp, symmetric = TRUE)
    psd = max(which(ev$values>0)) #Get the positive eigenvalues 
    evec = ev$vectors[,1:psd] #Get corresponding eigenvectors
    
    #Spectral decomposition of the partition of the covariance matrix that is full rank 
    Sp=evec[,1:psd]%*%diag(sqrt(ev$values[1:psd]))%*%t(evec[,1:psd])
    
    #Generate the predicted value from this distribution
    u=rnorm(np,0,1)
    fit_pred[i,]=out$mp+Sp%*%u
    
    #Progress Tracker
    cat(i/Npost*100," percent complete\r")
  }
  
  #Get summary statistics for the predictions
  pred_mean = apply(fit_pred, 2, mean)
  pred_sd = apply(fit_pred, 2, sd)
  
  out = list(fit_pred = fit_pred, pred_mean = pred_mean, pred_sd = pred_sd)
  return(out)
}


