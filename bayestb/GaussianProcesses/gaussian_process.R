#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Gaussian Processes:

#Description: Functions for Gaussian Process regression/emulation. File includes:  
#--Utility functions for normalizing data and computing pairwise differences
#--Covariance function(s) (add more functions)
#--Log priors and log posterior (both will depend on the covariance function)
#--MH step to accept or reject a new sample
#--Adaptive MH to help with the mixing during MCMC
#--GP function to fit a model
#--GP prediction function to compute predictions

#Status: Working. Need to add more covariance functions, create the necessary priors,
#--and add a feature in the logpost and/or GP functions to allow the user to specify
#--which covariance function to use.

#References:

#Notation:
#Lambda_y = random error precision
#Lambda_x = precision parameter multiplied to the covariance function 
#rho = GP correlation parameters using rho_covmat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Utility Functions:
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Covariance functions:
#rho correlation function
rho_covmat = function(diff_list, rho_vec, alpha = 2){
  R = matrix(1, nrow = nrow(data.frame(diff_list[[1]])), 
             ncol = ncol(data.frame(diff_list[[1]])))
  for(i in 1:length(rho_vec)){
    R = R*rho_vec[i]^(data.frame(diff_list[[i]])^alpha)
  }
  return(R)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Log Priors:
#Precision parameters -- gamma prior
prior_lam = function(lam, alpha, beta){
  den = log(dgamma(lam, shape = alpha, rate = beta))
  return(den)
}

#Correlation parameters -- beta prior (rho_covmat covariance function)
prior_rho = function(rho, alpha, beta){
  den = log(dbeta(rho, shape1 = alpha, shape2 = beta))
  return(den)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Log Posterior: (using normalized data & using rho_covmat)
log_post = function(y, ninfo, hp, lam_y, lam_x, rho, x_diff){
  #Check conditions
  if(any(rho<=0) || any(rho>=1) || lam_y<=0 || lam_x <= 0){
    return(-Inf)  
  }
  ny = ninfo$ny
  
  #Get covaraince matrices for the GP modeling the function and the random error
  Rx = rho_covmat(diff_list = x_diff, rho_vec = rho)
  epsilon = diag(1/lam_y, ny)
  
  #Combine these matrices to get one covaraince matrix
  V = (1/lam_x)*Rx + epsilon  
  
  #Get V inverse via Choleski Decomposition and log determinant
  #Check to see if all eigen values are positive to avoid an error
  if(any(eigen(V)$values < 0)){return(-Inf)}
  Vchol = chol(V)
  Vinv = chol2inv(Vchol)
  logdetV = 2*sum(log(diag(Vchol)))
  
  #Get Log-likelihood (ignoring constant of 2*pi)
  loglike = -0.5*logdetV - 0.5*t(y)%*%Vinv%*%y
  
  #Get the log priors for lambda 
  loglam_x = prior_lam(lam_x, hp$ax, hp$bx)
  loglam_y  = prior_lam(lam_y, hp$ay, hp$by)
  
  #Get log priors for rho 
  logrho = 0
  for(i in 1:ninfo$nx){
    logrho = logrho + prior_rho(rho[i], hp$rhoa[i], hp$rhob[i])
  }

  #Compute log posterior
  logpost = loglike + loglam_x + loglam_y + logrho
  return(logpost)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fitting Gaussian Process model
fit_gp = function(y, x_diff, mh_proposal, ninfo, hp, N, N_adapt){
  #Initialize posterior vectors
  post_lamy = 1
  post_lamx = 5
  post_rho = matrix(NA, nrow = N, ncol = ninfo$nx)
  
  post_rho[1,] = rep(0.5, ninfo$nx)
  
  #Initialize acceptance lists 
  accept_lamy = list(total = 0, last = 0)
  accept_lamx = list(total = 0, last = 0)
  accept_rho = list(total = rep(0, ninfo$nx), last = rep(0, ninfo$nx))
  
  #Simplify the names for the mh_proposal and ninfo
  rrho = mh_proposal$rrho; rx = mh_proposal$rx; ry = mh_proposal$ry; 
  nx = ninfo$nx
  
  #MCMC with Metropolis Hastings 
  for(i in 2:N){
    #Sample lambda_y -- get proposed move
    cur_val = post_lamy[i-1]
    new_val = runif(1, cur_val - ry, cur_val + ry)
    
    #Get densities 
    pdn = log_post(y, ninfo, hp, lam_x = post_lamx[i-1], lam_y = new_val, 
                   rho = post_rho[i-1,], x_diff)
    pdc = log_post(y, ninfo, hp, lam_x = post_lamx[i-1], lam_y = cur_val, 
                   rho = post_rho[i-1,], x_diff)
    
    #MH Step and update acceptance info
    mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
    post_lamy[i] = mh$theta0
    accept_lamy$total = accept_lamy$total + mh$accept
    
    #______________________________________________________
    #Sample lambda_x -- get proposed move
    cur_val = post_lamx[i-1]
    new_val = runif(1, cur_val - rx, cur_val + rx)
    
    #Get densities 
    pdn = log_post(y, ninfo, hp, lam_x = new_val, lam_y = post_lamy[i], 
                   rho = post_rho[i-1,], x_diff)
    pdc = log_post(y, ninfo, hp, lam_x = cur_val, lam_y = post_lamy[i], 
                   rho = post_rho[i-1,], x_diff)
    
    #MH Step and update acceptance info
    mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
    post_lamx[i] = mh$theta0
    accept_lamx$total = accept_lamx$total + mh$accept
    
    #______________________________________________________
    #Sample rho -- get proposed move
    cv = post_rho[i-1,] #Get the most recent rho values
    for(j in 1:nx){
      nv = cv
      nv[j] = runif(1, cv[j] - rrho[j], cv[j] + rrho[j])
      
      #Get densities 
      pdn = log_post(y, ninfo, hp, lam_x = post_lamx[i], lam_y = post_lamy[i], 
                     rho = nv, x_diff)
      
      pdc = log_post(y, ninfo, hp, lam_x = post_lamx[i], lam_y = post_lamy[i], 
                     rho = cv, x_diff)
      
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cv[j], new_val = nv[j])
      cv[j] = mh$theta0 #Update the jth element of cv
      accept_rho$total[j] = accept_rho$total[j] + mh$accept
    }
    post_rho[i,] = cv
    
    #______________________________________________________
    #Progress Tracker
    cat(i/N*100," percent complete\r")
    
    #Adaptive MCMC phase
    ss = N/20 #Set the step size
    if(i%%ss == 0 & i < (N_adapt+1)){
      #Adapt lambda_y
      if((accept_lamy$total-accept_lamy$last)/ss > 0.49 | (accept_lamy$total-accept_lamy$last)/ss < 0.39){
        #Apply adaptation and apply the results
        ar = adapt_delta(accept_lamy$total, accept_lamy$last, ss, ry)
        accept_lamy$last = ar$last_accept
        ry = ar$delta
      }else{accept_lamy$last = accept_lamy$total}
      
      #Adapt lambda_x
      if((accept_lamx$total-accept_lamx$last)/ss > 0.49 | (accept_lamx$total-accept_lamx$last)/ss < 0.39){
        #Apply adaptation and apply the results
        ar = adapt_delta(accept_lamx$total, accept_lamx$last, ss, rx)
        accept_lamx$last = ar$last_accept
        rx = ar$delta
      }else{accept_lamx$last = accept_lamx$total}
      
      #Adapt rho
      for(j in 1:nx){
        if((accept_rho$total[j]-accept_rho$last[j])/ss > 0.49 | (accept_rho$total[j]-accept_rho$last[j])/ss < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_rho$total[j], accept_rho$last[j], ss, rrho[j])
          accept_rho$last[j] = ar$last_accept
          rrho[j] = ar$delta
        }else{accept_rho$last[j] = accept_rho$total[j]}  
      }
      
    }
  }
  
  cat("\n Complete.\n")
  rate.rho=(accept_rho$total - accept_rho$last)/(N-N_adapt)
  rate.lambdax=(accept_lamx$total - accept_lamx$last)/(N-N_adapt)
  rate.lambday=(accept_lamy$total - accept_lamy$last)/(N-N_adapt)
  cat("[rate.rho =",rate.rho,"]\n")
  cat("[rate.lambdax =",rate.lambdax,"]\n")
  cat("[rate.lambday =",rate.lambday,"]\n")
  cat("[rrho =", rrho, "]\n")
  cat("[ry =", ry, "]\n")
  cat("[rx =", rx, "]\n")
  
  #Set the sample draws to save
  sdr = (N_adapt+1):N
  out = list(y = y, lambda_y = post_lamy[sdr], lambda_x = post_lamx[sdr],
             rho = post_rho[sdr,])
  return(out)
}

#Predictions with the GP 
#--pred_diff is np+ny X np+ny mattrix (order: test points, train points) 
predict_gp = function(fit, pred_diff, ninfo){
  #Set dimension parameters
  nx = ninfo$nx #Number of predictors
  ny = ninfo$ny #Number of observations
  np = ninfo$np #Number of predictions
  Npost = ninfo$Npost #Number of posterior draws
  nn = np + ny
  
  fit_pred = matrix(NA, nrow = Npost, ncol = np)
  #Make predictions based on the Npost posterior draws
  for(i in 1:Npost){
    #Now get the covariance matrices for eta and delta
    rho_vec = fit$rho[i,]
    R = rho_covmat(pred_diff, rho_vec) #pred_diff includes differences for training and testing points
    epsilon = 1/fit$lambda_y[i]*diag(c(rep(0,np),rep(1,ny))) #Predicting the real process, so no random error for predictions
    
    #Get the likelihood covariance matrix
    V = 1/fit$lambda_x[i]*R + epsilon
    V = as.matrix(V)
    #Get predictions using conditional normal identities
    V22 = V[(np+1):nn,(np+1):nn] #Covariance of training data
    V21 = V[(np+1):nn,1:np] #Covariance between train and pred
    V12 = V[1:np,(np+1):nn] #Covariance between pred and train
    V11 = V[1:np,1:np] #Covaraince of pred
    
    #Get inverse of V22 -- choleksi decomposition
    V22_chol = chol(V22)
    V22_inv = chol2inv(V22_chol)
    
    #Get mean and covaraince -- recall we assume mean 0 of the data
    mp = V12%*%V22_inv%*%(fit$y)
    Vp = V11 - V12%*%V22_inv%*%V21
    
    #Perform spectral decomposition on the new covariance matrix
    ev = eigen(Vp, symmetric = TRUE)
    psd = max(which(ev$values>0)) #Get the positive eigenvalues 
    evec = ev$vectors[,1:psd] #Get corresponding eigenvectors
    
    #Spectral decomposition of the partition of the covariance matrix that is full rank 
    Sp=evec[,1:psd]%*%diag(sqrt(ev$values[1:psd]))%*%t(evec[,1:psd])
    
    #Generate the predicted value from this distribution
    u=rnorm(np,0,1)
    fit_pred[i,]=mp+Sp%*%u
    
    #Progress Tracker
    cat(i/Npost*100," percent complete\r")
  }
  
  #Get summary statistics for the predictions
  pred_mean = apply(fit_pred, 2, mean)
  pred_sd = apply(fit_pred, 2, sd)
  
  out = list(y = fit$y, lambda_y = fit$lambda_y, lambda_x = fit$lambda_x,
             rho = fit$rho, fit_pred = fit_pred, pred_mean = pred_mean, pred_sd = pred_sd)
  return(out)
}
