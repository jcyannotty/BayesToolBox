#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Bayesian Calibration
#Description: Functions to perform Bayesian Calibration. Working as of 8/25/21

#References:

#Notation:
#z = (yf, yc) = vector of field data and computer output 
#rho_z = correlation parameters used in the covariance matrix of z (associated with eta)
#rho_del = correlation parameters used in covaraince matrix of discrepancy delta
#psi = correlation parameters used in the covariance of calibration parameters (associated with eta)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Data functions and parameters:
#Normalize Inputs and Output:
norm_x = function(design_c, design_f, num_cal){
  input_range = matrix(0, nrow = ncol(design_c), ncol = 2) 
  num_control = ncol(design_c) - num_cal #Get number of control inputs
  for(i in 1:ncol(design_c)){
    #Get range and scale the data
    rng = range(design_c[,i])
    if(i <= num_control){
      #Rescale both computer and field inputs
      design_c[,i] = (design_c[,i] - rng[1])/(rng[2] - rng[1])
      design_f[,i] = (design_f[,i] - rng[1])/(rng[2] - rng[1])
    }else{
      #Rescale just computer inputs
      design_c[,i] = (design_c[,i] - rng[1])/(rng[2] - rng[1])
    }
    #store the ranges 
    input_range[i,] = rng
  }
  
  out = list(dc = design_c, df = design_f, input_rng = input_range)
  return(out)
}

norm_y = function(yc, yf){
  mean_yc=mean(yc)
  sd_yc=sd(yc)
  yc=(yc-mean_yc)/sd_yc
  yf=(yf-mean_yc)/sd_yc
  
  out = list(yc_norm = yc, yf_norm = yf, mean_yc = mean_yc, sd_yc = sd_yc)
  return(out)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Covariance Matrix for GP:
#Get the pairwise differences between x vectors
pairwise_diff = function(design_mat){
  l_diff = list()
  for(i in 1:ncol(design_mat)){
    ltemp = list(outer(design_mat[,i],design_mat[,i],"-"))
    #names(ltemp) = paste0('m',i)
    l_diff[[i]] = ltemp
  }
  names(l_diff) = paste0('l',1:ncol(design_mat))
  return(l_diff)
}

#Covaraince matrix computation with the correlation parameter
rho_covmat = function(diff_list, rho_vec, alpha = 2){
  R = matrix(1, nrow = nrow(diff_list[[1]]), ncol = ncol(diff_list[[1]]))
  for(i in 1:length(rho_vec)){
    R = R*rho_vec[i]^(diff_list[[i]]^alpha)
  }
  return(R)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Log Priors:
#Precision parameters -- lambda_f, lambda_z, lambda_delta -- gamma prior
prior_lam = function(lam, alpha, beta){
  den = log(dgamma(lam, shape = alpha, rate = beta))
  return(den)
}

#Correlation parameters: rho_z, rho_delta, psi -- beta prior
prior_rho = function(rho, alpha, beta){
  den = log(dbeta(rho, shape1 = alpha, shape2 = beta))
  return(den)
}

#Calibration parameters theta -- uniform OR normal prior
prior_theta_unif = function(theta_vec){
  den = ifelse(any(theta_vec > 1) | any(theta_vec < 0), log(0), log(1))
  return(den)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Log Posterior: (using normalized data)
log_post = function(z, ninfo, hp, lam_f, lam_z, lam_del, rho_z, rho_del, psi, theta, z_diff, del_diff){
  #Check conditions
  if(any(rho_z<=0) || any(rho_z>=1) || any(psi<=0) || any(psi>=1) || any(rho_del<=0) || any(rho_del>=1) || lam_z<=0 || lam_f<=0 || lam_del<0){
    return(-Inf)  
  }
  nf = ninfo$nf
  
  #Get covaraince matrices for eta, delta, and random error components
  Rz = rho_covmat(diff_list = z_diff, rho_vec = c(rho_z, psi))
  Rdel = rho_covmat(diff_list = del_diff, rho_vec = rho_del)
  epsilon = diag(1/lam_f, nf)
  
  #Combine these matrices to get one covaraince matrix
  V = (1/lam_z)*Rz
  V[1:nf, 1:nf] = V[1:nf, 1:nf] + (1/lam_del)*Rdel + epsilon
  
  #Get V inverse via Choleski Decomposition and log determinant
  #Check to see if all eigen values are positive to avoid an error
  if(any(eigen(V)$values < 0)){return(-Inf)}
  Vchol = chol(V)
  Vinv = chol2inv(Vchol)
  logdetV = 2*sum(log(diag(Vchol)))
  
  #Get Log-likelihood (ignoring constant of 2*pi)
  loglike = -0.5*logdetV - 0.5*t(z)%*%Vinv%*%z
  
  #Get the log priors for lambda 
  loglam_f = prior_lam(lam_f, hp$af, hp$bf)
  loglam_z  = prior_lam(lam_z, hp$az, hp$bz)
  loglam_del = prior_lam(lam_del, hp$adel, hp$bdel)
  
  #Get log priors for rho and psi 
  logrho_z = logrho_del = logpsi = 0
  for(i in 1:ninfo$num_control){
    logrho_z = logrho_z + prior_rho(rho_z[i], hp$rhoa[i], hp$rhob[i])
    logrho_del = logrho_del + prior_rho(rho_del[i], hp$rdla[i], hp$rdlb[i])
  }
  for(i in 1:ninfo$num_cal){
    logpsi = logpsi + prior_rho(psi[i], hp$psia[i], hp$psib[i])
  }
  
  #Get log prior for theta
  logtheta = prior_theta_unif(theta_vec = theta)
  
  #Compute log posterior
  logpost = loglike + loglam_f + loglam_z + loglam_del + logrho_z + logrho_del + logpsi + logtheta
  return(logpost)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Metropolis Hastings functions:
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
#Calibration Function:
calibrate = function(z, z_diff, mh_proposal, ninfo, hp, N, N_adapt){
  #Initialize posterior vectors
  post_lamz = 1
  post_lamf = 50
  post_lamdel = 1
  post_rho = matrix(NA, nrow = N, ncol = ninfo$num_control)
  post_rhodel = matrix(NA, nrow = N, ncol = ninfo$num_control)
  post_psi = matrix(NA, nrow = N, ncol = ninfo$num_cal)
  post_theta = matrix(NA, nrow = N, ncol = ninfo$num_cal)
  
  post_rho[1,] = post_rhodel[1,] = rep(0.5, ninfo$num_control)
  post_psi[1,] = post_theta[1,] = rep(0.5, ninfo$num_cal)
  
  #Initialize acceptance lists 
  accept_lamz = list(total = 0, last = 0)
  accept_lamf = list(total = 0, last = 0)
  accept_lamdel = list(total = 0, last = 0)
  accept_rho = list(total = rep(0, ninfo$num_control), last = rep(0, ninfo$num_control))
  accept_rhodel = list(total = rep(0, ninfo$num_control), last = rep(0, ninfo$num_control))
  accept_psi = list(total = rep(0, ninfo$num_cal), last = rep(0, ninfo$num_cal))
  accept_theta = list(total = rep(0, ninfo$num_cal), last = rep(0, ninfo$num_cal))
  
  #Simplify the names for the mh_proposal and ninfo
  rrho = mh_proposal$rrho; rpsi = mh_proposal$rpsi; rrdel = mh_proposal$rrdel
  rf = mh_proposal$rf; rz = mh_proposal$rz; rdel = mh_proposal$rdel; rth = mh_proposal$rth
  ncal = ninfo$num_cal
  ncon = ninfo$num_control
  
  #Get discrepancy design pairwise differences
  del_diff = list()
  for(k in 1:ninfo$num_control){
    del_diff[[k]] = z_diff[[k]][1:ninfo$nf, 1:ninfo$nf]  
  }
  
  
  #MCMC with Metropolis Hastings 
  for(i in 2:N){
    #Sample lambda_z -- get proposed move
    cur_val = post_lamz[i-1]
    new_val = runif(1, cur_val - rz, cur_val + rz)
    
    #Get densities 
    pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i-1], lam_z = new_val, 
                   lam_del = post_lamdel[i-1], rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i-1], lam_z = cur_val, 
                   lam_del = post_lamdel[i-1], rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    #MH Step and update acceptance info
    mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
    post_lamz[i] = mh$theta0
    accept_lamz$total = accept_lamz$total + mh$accept
    
    #______________________________________________________
    #Sample lambda_f -- get proposed move
    cur_val = post_lamf[i-1]
    new_val = runif(1, cur_val - rf, cur_val + rf)
    
    #Get densities 
    pdn = log_post(z, ninfo, hp, lam_f = new_val, lam_z = post_lamz[i], 
                   lam_del = post_lamdel[i-1], rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    pdc = log_post(z, ninfo, hp, lam_f = cur_val, lam_z = post_lamz[i], 
                   lam_del = post_lamdel[i-1], rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    #MH Step and update acceptance info
    mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
    post_lamf[i] = mh$theta0
    accept_lamf$total = accept_lamf$total + mh$accept
    
    #______________________________________________________
    #Sample lambda_del -- get proposed move
    cur_val = post_lamdel[i-1]
    new_val = runif(1, cur_val - rdel, cur_val + rdel)
    
    #Get densities 
    pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                   lam_del = new_val, rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                   lam_del = cur_val, rho_z = post_rho[i-1,], rho_del = post_rhodel[i-1,],
                   psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
    #MH Step and update acceptance info
    mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cur_val, new_val = new_val)
    post_lamdel[i] = mh$theta0
    accept_lamdel$total = accept_lamdel$total + mh$accept
    
    #______________________________________________________
    #Sample rho_z -- get proposed move
    cv = post_rho[i-1,] #Get the most recent rho values
    for(j in 1:ncon){
      nv = cv
      nv[j] = runif(1, cv[j] - rrho[j], cv[j] + rrho[j])
      
      #Get densities 
      pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del = post_lamdel[i], rho_z = nv, rho_del = post_rhodel[i-1,],
                     psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
      pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del =  post_lamdel[i], rho_z = cv, rho_del = post_rhodel[i-1,],
                     psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cv[j], new_val = nv[j])
      cv[j] = mh$theta0 #Update the jth element of cv
      accept_rho$total[j] = accept_rho$total[j] + mh$accept
    }
    post_rho[i,] = cv
    
    #______________________________________________________
    #Sample rho_del -- get proposed move
    cv = post_rhodel[i-1,] #Get the most recent rho values
    for(j in 1:ncon){
      nv = cv
      nv[j] = runif(1, cv[j] - rrdel[j], cv[j] + rrdel[j])
      
      #Get densities 
      pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del =  post_lamdel[i], rho_z = post_rho[i,], rho_del = nv,
                     psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
      pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del =  post_lamdel[i], rho_z = post_rho[i,], rho_del = cv,
                     psi = post_psi[i-1,], theta = post_theta[i-1,], z_diff, del_diff)
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cv[j], new_val = nv[j])
      cv[j] = mh$theta0 #Update the jth element of cv
      accept_rhodel$total[j] = accept_rhodel$total[j] + mh$accept
    }
    post_rhodel[i,] = cv
    #______________________________________________________
    #Sample psi -- get proposed move
    cv = post_psi[i-1,] #Get the most recent rho values
    for(j in 1:ncal){
      nv = cv
      nv[j] = runif(1, cv[j] - rpsi[j], cv[j] + rpsi[j])
      
      #Get densities 
      pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del = post_lamdel[i], rho_z = post_rho[i,], rho_del = post_rhodel[i,],
                     psi = nv, theta = post_theta[i-1,], z_diff, del_diff)
      pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del = post_lamdel[i], rho_z = post_rho[i,], rho_del = post_rhodel[i,],
                     psi = cv, theta = post_theta[i-1,], z_diff, del_diff)
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cv[j], new_val = nv[j])
      cv[j] = mh$theta0 #Update the jth element of cv
      accept_psi$total[j] = accept_psi$total[j] + mh$accept
    }
    post_psi[i,] = cv
    #______________________________________________________
    #Sample theta -- get proposed move
    cv = post_theta[i-1,] #Get the most recent rho values
    for(j in 1:ncal){
      nv = cv
      nv[j] = runif(1, cv[j] - rth[j], cv[j] + rth[j])
      
      #Get densities 
      pdn = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del = post_lamdel[i], rho_z = post_rho[i,], rho_del = post_rhodel[i,],
                     psi = post_psi[i,], theta = nv, z_diff, del_diff)
      pdc = log_post(z, ninfo, hp, lam_f = post_lamf[i], lam_z = post_lamz[i], 
                     lam_del = post_lamdel[i], rho_z = post_rho[i,], rho_del = post_rhodel[i,],
                     psi = post_psi[i,], theta = cv, z_diff, del_diff)
      #MH Step and update acceptance info
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = cv[j], new_val = nv[j])
      cv[j] = mh$theta0 #Update the jth element of cv
      accept_theta$total[j] = accept_theta$total[j] + mh$accept
    }
    post_theta[i,] = cv
    
    #Progress Tracker
    cat(i/N*100," percent complete\r")
    
    #Adaptive MCMC phase
    ss = N/20 #Set the step size
    if(i%%ss == 0 & i < (N_adapt+1)){
      #Adapt lambda_z
      if((accept_lamz$total-accept_lamz$last)/ss > 0.49 | (accept_lamz$total-accept_lamz$last)/ss < 0.39){
        #Apply adaptation and apply the results
        ar = adapt_delta(accept_lamz$total, accept_lamz$last, ss, rz)
        accept_lamz$last = ar$last_accept
        rz = ar$delta
      }else{accept_lamz$last = accept_lamz$total}
      
      #Adapt lambda_f
      if((accept_lamf$total-accept_lamf$last)/ss > 0.49 | (accept_lamf$total-accept_lamf$last)/ss < 0.39){
        #Apply adaptation and apply the results
        ar = adapt_delta(accept_lamf$total, accept_lamf$last, ss, rf)
        accept_lamf$last = ar$last_accept
        rf = ar$delta
      }else{accept_lamf$last = accept_lamf$total}
      
      #Adapt lambda_del
      if((accept_lamdel$total-accept_lamdel$last)/ss > 0.49 | (accept_lamdel$total-accept_lamdel$last)/ss < 0.39){
        #Apply adaptation and apply the results
        ar = adapt_delta(accept_lamdel$total, accept_lamdel$last, ss, rdel)
        accept_lamdel$last = ar$last_accept
        rdel = ar$delta
      }else{accept_lamdel$last = accept_lamdel$total}
      
      #Adapt rho_z
      for(j in 1:ncon){
        if((accept_rho$total[j]-accept_rho$last[j])/ss > 0.49 | (accept_rho$total[j]-accept_rho$last[j])/ss < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_rho$total[j], accept_rho$last[j], ss, rrho[j])
          accept_rho$last[j] = ar$last_accept
          rrho[j] = ar$delta
        }else{accept_rho$last[j] = accept_rho$total[j]}  
      }
      
      #Adapt rho_del
      for(j in 1:ncon){
        if((accept_rhodel$total[j]-accept_rhodel$last[j])/ss > 0.49 | (accept_rhodel$total[j]-accept_rhodel$last[j])/ss < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_rhodel$total[j], accept_rhodel$last[j], ss, rrdel[j])
          accept_rhodel$last[j] = ar$last_accept
          rrdel[j] = ar$delta
        }else{accept_rhodel$last[j] = accept_rhodel$total[j]}  
      }
      
      #Adapt psi
      for(j in 1:ncal){
        if((accept_psi$total[j]-accept_psi$last[j])/ss > 0.49 | (accept_psi$total[j]-accept_psi$last[j])/ss < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_psi$total[j], accept_psi$last[j], ss, rpsi[j])
          accept_psi$last[j] = ar$last_accept
          rpsi[j] = ar$delta
        }else{accept_psi$last[j] = accept_psi$total[j]}  
      }
      
      #Adapt theta
      for(j in 1:ncal){
        if((accept_theta$total[j]-accept_theta$last[j])/ss > 0.49 | (accept_theta$total[j]-accept_theta$last[j])/ss < 0.39){
          #Apply adaptation and apply the results
          ar = adapt_delta(accept_theta$total[j], accept_theta$last[j], ss, rth[j])
          accept_theta$last[j] = ar$last_accept
          rth[j] = ar$delta
        }else{accept_theta$last[j] = accept_theta$total[j]}  
      }
      
    }
  }
  
  cat("\n Complete.\n")
  rate.rhoz=(accept_rho$total - accept_rho$last)/(N-N_adapt)
  rate.psi=(accept_psi$total - accept_psi$last)/(N-N_adapt)
  rate.rhodel=(accept_rhodel$total - accept_rhodel$last)/(N-N_adapt)
  rate.lambdaf=(accept_lamf$total - accept_lamf$last)/(N-N_adapt)
  rate.lambdaz=(accept_lamz$total - accept_lamz$last)/(N-N_adapt)
  rate.lambdadel=(accept_lamdel$total - accept_lamdel$last)/(N-N_adapt)
  rate.theta=(accept_theta$total - accept_theta$last)/(N-N_adapt)
  cat("[rate.rhoz=",rate.rhoz,"]\n")
  cat("[rate.psi=",rate.psi,"]\n")
  cat("[rate.rhodel=",rate.rhodel,"]\n")
  cat("[rate.lambdaf=",rate.lambdaf,"]\n")
  cat("[rate.lambdaz=",rate.lambdaz,"]\n")
  cat("[rate.lambdadel=",rate.lambdadel,"]\n")
  cat("[rate.theta=",rate.theta,"]\n\n")
  
  # cat("Accept Total lamdaf:",accept_lamf$total, '\n')
  # cat("Accept Total lamdaz:",accept_lamz$total, '\n')
  # cat("Accept Total lamdadel:",accept_lamdel$total, '\n\n')
  # 
  # cat("Last Accept Total lamdaf:",accept_lamf$last, '\n')
  # cat("Last Accept Total lamdaz:",accept_lamz$last, '\n')
  # cat("Last Accept Total lamdadel:",accept_lamdel$last, '\n\n')
  # 
  # cat("Proposal delta lamdaf:",rf, '\n')
  # cat("Proposal delta lamdaz:",rz, '\n')
  # cat("Proposal delta lamdadel:",rdel, '\n')
  
  #Set the sample draws to save
  sdr = (N_adapt+1):N
  out = list(z = z, lambda_z = post_lamz[sdr], lambda_f = post_lamf[sdr], lambda_del = post_lamdel[sdr],
             rho_z = post_rho[sdr,], rho_del = post_rhodel[sdr,], psi = post_psi[sdr,],
             theta = post_theta[sdr,])
  return(out)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prediction functions:
#pred_diff must be in the order of Predictions, Field Obs, and Computer Output (Xp, Xf, Xc)
predict_calibrate = function(fit, pred_diff, ninfo){
  #Set dimension parameters
  nf = ninfo$nf #Number of field obs
  nc = ninfo$nc #Number of computer outputs
  nz = nf+nc #Size of training data
  np = ninfo$np #Size of prediction set
  nn = np + nz #Size of training and prediction data
  ncon = ninfo$num_control #Number of control inputs
  ncal = ninfo$num_cal #Number of calibration inputs
  Npost = length(fit$lambda_f) #Number of posterior draws
  
  fit_pred = fit_eta = matrix(NA, nrow = Npost, ncol = np)
  #Make predictions based on the Npost posterior draws
  for(i in 1:Npost){
    #Update the pairwise difference matrix: 
    #--field obs and predictions will have 0 for the calibration differences
    #--Only the diff between the calibration parameters for the computer output settings will be effected
    for(k in 1:ncal){
      pred_diff[[k+ncon]][1:(nf+np),] = pred_diff[[k+ncon]][1:(nf+np),] + fit$theta[i,k]
      pred_diff[[k+ncon]][,1:(nf+np)] = pred_diff[[k+ncon]][,1:(nf+np)] - fit$theta[i,k]
    }
    
    #Get the delta matrix of differences for field obs and predictions
    del_diff = list()
    for(k in 1:ncon){
      del_diff[[k]] = z_diff[[k]][1:(nf+np), 1:(nf+np)]  
    }
    
    #Now get the covariance matrices for eta and delta
    rho_vec = c(fit$rho_z[i,], fit$psi[i,])
    R = rho_covmat(pred_diff, rho_vec)
    Rdel = rho_covmat(del_diff, fit$rho_del[i,])
    epsilon = 1/fit$lambda_f[i]*diag(c(rep(0,np),rep(1,nf),rep(0,nc))) #Predicting the real process, so no random error for predictions
    
    #Get the likelihood covariance matrix
    V = 1/fit$lambda_z[i]*R + epsilon
    V[1:(nf+np),1:(nf+np)] = V[1:(nf+np),1:(nf+np)] + 1/fit$lambda_del[i]*Rdel #Add discrepancy cov mat
    
    #Get predictions using conditional normal identities
    V22 = V[(np+1):nn,(np+1):nn] #Covariance of training data
    V21 = V[(np+1):nn,1:np] #Covariance between train and pred
    V12 = V[1:np,(np+1):nn] #Covariance between pred and train
    V11 = V[1:np,1:np] #Covaraince of pred
    
    #Get inverse of V22 -- choleksi decomposition
    V22_chol = chol(V22)
    V22_inv = chol2inv(V22_chol)
    
    #Get mean and covaraince -- recall we assume mean 0 of the data
    mp = V12%*%V22_inv%*%(fit$z)
    Vp = V11 - V12%*%V22_inv%*%V21
    
    #Perform spectral decomposition on the new covariance matrix
    ev = eigen(Vp, symmetric = TRUE)
    psd=max(which(ev$values>0)) #Get the positive eigenvalues 
    evec = ev$vectors[,1:psd] #Get corresponding eigenvectors
    
    #Spectral decomposition of the partition of the covariance matrix that is full rank 
    Sp=evec[,1:psd]%*%diag(sqrt(ev$values[1:psd]))%*%t(evec[,1:psd])
    
    #Generate the predicted value from this distribution
    u=rnorm(np,0,1)
    fit_pred[i,]=mp+Sp%*%u
    
    #Generate predictions for eta the computer model -- predictions no longer have discrepancy term  
    #Get the delta matrix of differences for field obs
    del_diff = list()
    for(k in 1:ncon){
      del_diff[[k]] = z_diff[[k]][(np+1):(nf+np), (np+1):(nf+np)]  
    }
    Rdel = rho_covmat(del_diff, fit$rho_del[i,])
    
    #Get the likelihood covariance matrix
    V = 1/fit$lambda_z[i]*R + epsilon
    V[(np+1):(nf+np),(np+1):(nf+np)] = V[(np+1):(nf+np),(np+1):(nf+np)] + 1/fit$lambda_del[i]*Rdel #Add discrepancy cov mat
    
    #Get predictions using conditional normal identities
    V22 = V[(np+1):nn,(np+1):nn] #Covariance of training data
    V21 = V[(np+1):nn,1:np] #Covariance between train and pred
    V12 = V[1:np,(np+1):nn] #Covariance between pred and train
    V11 = V[1:np,1:np] #Covaraince of pred
    
    #Get inverse of V22 -- choleksi decomposition
    V22_chol = chol(V22)
    V22_inv = chol2inv(V22_chol)
    
    #Get mean and covaraince -- recall we assume mean 0 of the data
    mp = V12%*%V22_inv%*%(fit$z)
    Vp = V11 - V12%*%V22_inv%*%V21
    
    #Perform spectral decomposition on the new covariance matrix
    ev = eigen(Vp, symmetric = TRUE)
    psd=max(which(ev$values>0)) #Get the positive eigenvalues 
    evec = ev$vectors[,1:psd] #Get corresponding eigenvectors
    
    #Spectral decomposition of the partition of the covariance matrix that is full rank 
    Sp=evec[,1:psd]%*%diag(sqrt(ev$values[1:psd]))%*%t(evec[,1:psd])
    
    #Generate the predicted value from this distribution
    u=rnorm(np,0,1)
    fit_eta[i,]=mp+Sp%*%u
  }
  
  #Get summary statistics for the predictions
  pred_mean = apply(fit_pred, 2, mean)
  pred_sd = apply(fit_pred, 2, sd)
  eta_mean = apply(fit_eta, 2, mean)
  eta_sd = apply(fit_eta, 2, sd)
  
  out = list(z = fit$z, lambda_z = fit$lambda_z, lambda_f = fit$lambda_f, lambda_del = fit$lambda_del,
             rho_z = fit$rho_z, rho_del = fit$rho_del, psi = fit$psi,theta = fit$theta,
             fit_pred = fit_pred, fit_eta = fit_eta, pred_mean = pred_mean, pred_sd = pred_sd,
             eta_mean = eta_mean, eta_sd = eta_sd)
  return(out)
}

