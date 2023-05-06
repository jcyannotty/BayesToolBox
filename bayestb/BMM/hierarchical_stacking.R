#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hierarchical Stacking Functions:
# Implemented based on Yao et al. 2021

#Description: This script contains functions to perform hierarchical stacking. Main features include:
#--Pareto Smoothed Importance Sampling (PSIS)
#--Prior and log posterior distribution for HS parameters
#--Function to fit a HS model and return the posterior distribution of all HS parameters
#--Function to convert the HS parameter posteriors into a posterior distribution for the weights
#--Visualization tools (Still in progress)
#--Obtain Predictions using the HS weights -- i.e. predictions with model mixing
#--The functions are general, hence we should be able to use any design matrix for modeling the unconstrained weights wstar
#--An intercept must be included in the HS design. Pay attention to the required order of parameters (mu, sig2, alphas)

#Status: Works -- need to add an adaptive stage to the fit_hs function. 
#--Trace plots are sometimes loaded with autocorrelation and other times they are fine. 

#References:
#https://vasishth.github.io/bayescogsci/book/expected-log-predictive-density-of-a-model.html
#https://jrnold.github.io/bayesian_notes/model-comparison.html
#https://cran.r-project.org/web/packages/loo/loo.pdf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 1: Setup 
#setwd('C:/Users/johny/Documents/Ohio State/Research/Bayesian Hierarchical Stacking')
library(MASS)
library(car)
library(loo)
library(ggplot2)
library(reshape2)
library(fdrtool)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Block 2: Density functions
#Section 2.1: Priors for the stacking weight parameters and hyperparameters 
#prior of alpha_mk ~ normal with mean 0 and Variance sig2_k
prior_alpha = function(a_mk, sig2_k, a_mean = 0){
  a_den = dnorm(a_mk, mean = a_mean, sd = sqrt(sig2_k))
  return(a_den)
}

#Prior of mu_k ~ Normal with mean mu_0 and variance tau2_0
prior_mu_k = function(mu_k, mu_0 = 0, tau2_0 = 1){
  m_den = dnorm(mu_k, mean = mu_0, sd = sqrt(tau2_0))
  return(m_den)
}

#Prior of sigma2_k ~ Truncated normal on (0, infty) with mean 
#----Simplified to a gamma distribution
prior_sig2_k = function(sig2_k, shape_0 = 1, scale_0 = 0.5, HalfNorm = FALSE){
  if(HalfNorm){
    sig2 = shape_0*scale_0
    theta = sqrt(pi/(2*sig2))
    s_den = dhalfnorm(sig2_k, theta = theta)  
  }else{
    s_den = dgamma(sig2_k, shape = shape_0, scale = scale_0) 
  }
  return(s_den)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 2.2: Log Joint Posterior of mu sigma2, and alpha (order matters) 
#Requirements for the function: 
#--loo_pred: loo predictive density matrix using Block 3 (nxK), 
#--hsp_matrix: hierarchical stacking parameter list that contains the most recent values of parameters stored in hsp_array, an array of dimension (N_hs x m x K-1), where N_hs = #MCMC draws, m = #hs.parameters 
#--hsp_basis: basis function for the hs parameters -- requires an intercept in first column. This is a design matrix nxm
#--hsp_matrix must contain vectors with the 1st two parameters being mu and sigma
#--hp: hyper parameter list
#--Check the intercept condition for basis function sometime before starting the MCMC - waste to check every iteration
log_post = function(loo_pd, hsp_matrix, hsp_basis, hp){
  #Get dimensions
  u = ncol(hsp_matrix) - 1 #number of unconstrained weights 
  J = ncol(hsp_basis) - 1 #number of alphas in the weight model
  K = u + 1 #number of models
  
  #Get the prior density per model and compute the model weights
  log_prior = 0
  w_star = matrix(0, nrow = nrow(hsp_basis), ncol = K) #This fixes the constrained wstar to 0 
  for(r in 1:u){
    a_vec = hsp_matrix[-c(1:2),r]
    sigk = hsp_matrix[2,r]
    muk = hsp_matrix[1,r]
    
    #Get log priors for the individual parameters
    log_ap = sum(log(prior_alpha(a_vec, sigk, a_mean = hp$alpha_mean[,r])))
    log_sp = log(prior_sig2_k(sigk, shape_0 = hp$sig_shape[r], scale_0 = hp$sig_scale[r]))
    log_mk = log(prior_mu_k(muk, mu_0 = hp$mu_mean[r], tau2_0 = hp$mu_var[r]))
    
    #Get the log prior for the model
    log_prior[r] = log_ap + log_sp + log_mk
    
    #Compute the unconstrained weights with the linear model
    bk = matrix(c(muk,a_vec), ncol = 1)
    ws = hsp_basis %*% bk
    w_star[,r] = ws 
  }
  #Apply the softmax transformation to get the weights on [0,1] scale -- this is a matrix n X K
  wts = exp(w_star)/rowSums(exp(w_star))
  
  #Get the loo predictive density ("likelihood")
  log_lik = sum(log(rowSums(loo_pd*wts))) #using pointwise multiplication
  
  #Get the final density
  log_den = log_lik + sum(log_prior)
  return(log_den)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Block 3: LOO and PSIS
#Section 3.1: PSIS Weights
#Requires the densities of each training pt from the MCMC posterior draws for each model
#yt_den = Density of training points, N x n x K  (N = #posterior draws, n = sample size, K = #models) 
get_psis_wts = function(yt_den){
  K = dim(yt_den)[3]
  ratios = 1/yt_den 
  
  #Initialize the psis array N x n x K array (N = mcmc size, n = train sample size, k = # models) 
  psis_weights = array(NA, dim = dim(ratios))
  for(j in 1:K){
    reff = relative_eff(1/ratios[,,j], chain_id = rep(1,nrow(ratios)))
    psis_obj = psis(log(ratios[,,j]), r_eff = reff)
    
    if(sum(psis_obj$diagnostics$pareto_k > 0.7) != 0){
      print(paste("Warning, unstable k. Model", j, '-- Max k =', 
                  round(max(psis_obj$diagnostics$pareto_k),3)))
      #print("PSIS Weights are not used for and observation with k > 0.7. Use default weights of 1 for this given obs.")
      h = which(psis_obj$diagnostics$pareto_k > 0.7)
    }else{
      print(paste("PSIS weights for Model", j, '= Complete'))
    }
    #Converting from log weights to normal scale
    psis_weights[,,j] = weights(psis_obj, log = FALSE, normalize = FALSE)  
  }
  
  
  return(psis_weights)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 3.2: Get the LOO densities
#loo_pred_den is n_train X K -- row i k represents p(Yi | Y(-i), Mk)
#psis_weights -- output from get_psis_wts function
get_loo_den = function(yt_den, psis_weights){
  K = dim(yt_den)[3]
  n = ncol(yt_den)
  loo_pred_den = matrix(NA, nrow=n, ncol = K)
  for(j in 1:K){
    loo_num = apply(yt_den[,,j]*psis_weights[,,j],2,sum)
    sum_w = apply(psis_weights[,,j],2,sum)
    loo_pred_den[,j] = loo_num/sum_w
  }
  
  #Each column corresponds to the loo desnisities under a different model
  colnames(loo_pred_den) = paste0('model', 1:K) 
  return(loo_pred_den)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Block 4: Metropolis Hastings Functions
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
      #if(alpha < -4){print(paste0("alpha = ",alpha, " -- U = ",u))}
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
#Block 5: Hierarchical Stacking MCMC 
#Output of fit_hs is a K-dim array of posterior draws for mu, sig2, and alpha's 
fit_hs = function(loo_pd, mh_proposal, hsp_basis, hp, N, N_burn){
  #Get dimensions
  J = ncol(hsp_basis) - 1 #number of alphas \
  K = ncol(loo_pd) #number of models
  n = nrow(loo_pd) #sample size--number of training pts
  m = J + 2 #number of parameters (including mu and sig2)
  u = K-1 #number of unconstrained models
  
  #Initialize the hs posterior array for the parameters used to model the unconstrained weights
  hsp_array = array(1, dim = c(N,m,K))
  hsp_accept = matrix(0, nrow = m, ncol = K-1)
  
  #Initialize the 1st row elements
  for(j in 1:u){
    hsp_array[1,,j] = runif(m, min = 0.8, max = 1.2) #Random Initialization -- easier to debug
  }
  
  #Shorten names of proposal steps from mh_proposal 
  #--note, adaptive mcmc is not used at the moment. It should be an option, 
  #--but leaving it as the default option will could lead to serious run time issues
  rmu = mh_proposal$rmu; ra = mh_proposal$ra; rsig = mh_proposal$rsig
  
  #Construct the MCMC
  for(i in 2:N){
    #Get new and current parameter matrices (note, fixing a row of an array returns a matrix where the kth column of the matrix is the ith row in the kth layer of the array)
    c_mat = n_mat = hsp_array[i-1,,]
    for(j in 1:u){
      #Sample mu_j -- get proposed move
      cv = hsp_array[i-1,1,j] #Current Value
      nv = runif(1, cv - rmu, cv + rmu) #new value
      
      #Update the point in n_mat
      n_mat[1,j] = nv
      
      #get proposed density
      pdc = log_post(loo_pd, hsp_matrix = c_mat, hsp_basis, hp)
      pdn = log_post(loo_pd, hsp_matrix = n_mat, hsp_basis, hp)
      
      #MH Step 
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = c_mat, new_val = n_mat)
      c_mat = mh$theta0
      hsp_accept[1,j] = hsp_accept[1,j] + mh$accept
      
      #____________________________________________________
      #Sample sigma2j
      cv = c_mat[2,j]
      nv = runif(1, cv - rsig, cv + rsig)
      
      #Update the new matrix
      n_mat = c_mat
      n_mat[2,j] = nv
      
      #Get the densities (corrects for the potential of getting -variance term)
      pdc = log_post(loo_pd, hsp_matrix = c_mat, hsp_basis, hp)
      pdn = ifelse(nv>0, log_post(loo_pd, hsp_matrix = n_mat, hsp_basis, hp), -Inf)
      
      #MH Step 
      mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = c_mat, new_val = n_mat)
      c_mat = mh$theta0
      hsp_accept[2,j] = hsp_accept[2,j] + mh$accept
      
      #____________________________________________________
      #Sample alpha_rj
      for(r in 3:(J+2)){
        #Get new and current values of the alpha
        n_mat = c_mat
        cv = c_mat[r,j]
        lb = min(cv - ra, ra + cv)
        ub = max(cv - ra, ra + cv)
        nv = runif(1, lb, ub)
        n_mat[r,j] = nv
        
        #Get the densities
        pdc = log_post(loo_pd, hsp_matrix = c_mat, hsp_basis, hp)
        pdn = log_post(loo_pd, hsp_matrix = n_mat, hsp_basis, hp)
        
        #MH Step 
        mh = mh_logstep(pdc = pdc, pdn = pdn, current_val = c_mat, new_val = n_mat)
        c_mat = mh$theta0
        hsp_accept[r,j] = hsp_accept[r,j] + mh$accept
      }
     
    }
    #Track progress
    cat("Progress:", i/N*100, '% complete\r')
    
    #Update the parameter array 
    hsp_array[i,,] = c_mat
    
  }
  
  #Package Results:
  #Set sample draws to save
  #sdr = (N_adapt+1):N
  sdr = (N_burn + 1):N
  out_post = hsp_array[sdr,,]
  colnames(out_post) = c('mu_k', 'sig2_k', paste0('a',1:J,'_k'))
  rownames(hsp_accept) = colnames(out_post)
  out = list(post = out_post, accept = hsp_accept)
  return(out)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 5.2: Get model weights
#Get Weights function -- uses the output from the fit_hs function
#--Outputs a list of two matrices (1) the point estimate for the model weights  
#--and (2) the sd of the weights. Each matrix is n x K
#--This function gets the weights based on the posterior distribution of w
get_wts = function(fit_hs_output, hsp_basis, q.upper = 0.975, q.lower = 0.025){
  #Get the mu and alpha posterior draws
  b_array = fit_hs_output[,-2,]
  
  #Get number of models (constrained and unconstrained)
  K = dim(b_array)[[3]]
  #Set the wstar array and matrix
  wstar_array = array(0, dim = c(nrow(hsp_basis), nrow(b_array),K))
  for(j in 1:(K-1)){
    for(i in 1:nrow(b_array)){
      b_vec = matrix(b_array[i,,j], ncol = 1)
      wstar_array[,i,j] = hsp_basis%*%b_vec    
    }
  }
  
  #Prepare for the softmax transformation (these are n X Npost x K arrays)
  w_array_num = exp(wstar_array) #numerator 
  w_array_den = rowSums(w_array_num, dim = 2) #denominator
  
  #Perform softmax transformation for each dimension of the array. Then obtain the posterior info and store
  wts_x = matrix(0, nrow = nrow(hsp_basis), ncol = K)
  wts_x_sd = matrix(0, nrow = nrow(hsp_basis), ncol = K)
  wts_x_lb = matrix(0, nrow = nrow(hsp_basis), ncol = K)
  wts_x_ub = matrix(0, nrow = nrow(hsp_basis), ncol = K)
  wts_post = array(0, dim = dim(w_array_num))
  for(i in 1:K){
    wts_post[,,i] = w_array_num[,,i]/w_array_den
    wts_x[,i] = apply(w_array_num[,,i]/w_array_den, 1, mean)
    wts_x_sd[,i] = apply(w_array_num[,,i]/w_array_den, 1, sd)
    wts_x_lb[,i] = apply(w_array_num[,,i]/w_array_den, 1, quantile, q.lower)
    wts_x_ub[,i] = apply(w_array_num[,,i]/w_array_den, 1, quantile, q.upper)
  }
  
  out = list(wts_mean = wts_x, wts_sd = wts_x_sd, wts_lb = wts_x_lb, wts_ub = wts_x_ub, wts_post = wts_post)
  return(out)
}

#Get weights using the posterior mean of the alpha and mu parameters
#--This is alternative way to get the weights. This is analogous with how you would fit
#--a Bayesian logistic regression model
get_wts_lm = function(fit_hs_output, hsp_basis){
  #Get the number of models
  K = dim(fit_hs_output)[3]
  
  #Get the mu and alpha posterior draws
  b_array = fit_hs_output[,-2,]
  
  #Store the posterior means in the beta matrix (each column is a different model)
  b_matrix = matrix(0, nrow = ncol(b_array), ncol = K)
  for(j in 1:K){
    b_matrix[,j] = apply(b_array[,,j], 2, mean)
  }
  
  #Cast as a matrix and multiply with hsp_basis to get avg wstar
  b_matrix = as.matrix(b_matrix)
  wstar_mat = hsp_basis%*%b_matrix
  
  #Prepare for the softmax transformation (these are n X Npost x K arrays)
  w_mat_num = exp(wstar_mat) #numerator 
  w_mat_den = rowSums(w_mat_num) #denominator
  
  #Perform softmax transformation for each dimension of the array. Then obtain the posterior mean and store
  wts_x = matrix(0, nrow = nrow(hsp_basis), ncol = K)
  for(i in 1:K){
    wts_x[,i] = w_mat_num[,i]/w_mat_den
  }
  
  return(wts_x)
}

#Get Credible Interval for Weights
#--Input the posterior distribution of weights
get_wts_ci = function(wts_post, ci_level = 0.95){
  K = dim(wts_post)[3] #Number of models
  n = nrow(wts_post) #Number of obs
  lb = (1-ci_level)/2 #Lower percentile
  ub = 1 - lb #Upper percentile
  wts_lb = wts_ub = matrix(0, nrow = n, ncol = K) #Lower and Upper bounds
  for(i in 1:K){
    wts_lb[,i] = apply(wts_post[,,i], 1, function(x) quantile(x, lb))
    wts_ub[,i] = apply(wts_post[,,i], 1, function(x) quantile(x, ub))
  }
  wts_list = list(wts_lb = wts_lb, wts_ub = wts_ub)
  return(wts_list)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 5.3: Mixed predictive distribution
#Get predictive distributions with the HS weights.
#--pred_dist: Npost x n x K array of posterior draws
#--weights are applied to obtain the final mixture distribution
predict_hs = function(pred_dist, hs_wts){
  #use the property of arrays -- fix a row and the result is a matrix where each column k corresponds to the row vector of the kth layer of the array 
  y_mix = matrix(0, nrow = nrow(pred_dist), ncol = ncol(pred_dist))
  for(j in 1:dim(pred_dist)[[3]]){
    xx = apply(pred_dist[,,j], 1, function(x) x*hs_wts[,j])
    y_mix = t(xx) + y_mix
  }
  return(y_mix)
}

#Mixed function using Bayesian HS -- must be known functions
#--f_matrix = matrix of function output -- n x K
#--wts_mean = output of get_wts_bc - n x K matrix
#mix_functions = function(f_matrix, wts_mean){
#  #Pointwise multiplication
#  f_mix = rowSums(f_matrix*wts_mean)
#  return(f_mix)
#}

#Mixed Function credible interval
#--f_matrix = matrix of function output
#--wts_post = posterior dist of weights n x Npost x K
mix_functions = function(f_matrix, wts_post, ci_level = 0.95){
  K = dim(wts_post)[3] #Number of models
  n = nrow(wts_post) #Number of obs
  lb = (1-ci_level)/2 #Lower percentile
  ub = 1 - lb #Upper percentile
  f_post = array(0, dim = dim(wts_post)) #posterior dist of f
  
  #Get the posterior dist of f
  for(i in 1:K){
    f_post[,,i] = apply(wts_post[,,i], 2, function(x) x*f_matrix[,i])
  }
  
  #Posterior Distribution of the mixed functions
  f_mix_post = rowSums(f_post, dim = 2)
  f_lb = apply(f_mix_post,1, quantile, lb)
  f_ub = apply(f_mix_post,1, quantile, ub)
  f_mean = rowMeans(f_mix_post)
  
  #Return results
  f_mix_list = list(f_lb = f_lb, f_ub = f_ub, f_mean = f_mean, f_post = f_mix_post) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Section 5.4: Get Expected Log Pointwise Density Estimates
get_elpd = function(pred_den, hs_wts){
  #use the property of arrays -- fix a row and the result is a matrix where each column k corresponds to the row vector of the kth layer of the array 
  y_mix = matrix(0, nrow = nrow(pred_den), ncol = ncol(pred_den))
  for(j in 1:dim(pred_den)[[3]]){
    xx = apply(log(pred_den[,,j]), 1, function(x) x*hs_wts[,j])
    y_mix = t(xx) + y_mix
  }
  return(y_mix)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Block 6: Visualizing the weights
#Section 6.1: Plot weights vs. an input
plot_wts_x = function(x_vec, wts, x_name){
  data = data.frame(rep(x_vec, ncol(wts)), melt(wts))
  colnames(data) = c('x', 'x_row', 'Model', 'Wt_x')
  data$Model = as.factor(data$Model)
  y_lab = 'W(X)'
  g_title = paste('Model Weights vs.', x_name)
  col_list = c('red', 'blue', 'green', 'purple', 'orange', 'pink')
  
  p = ggplot(data = data, aes(x = x, y = Wt_x, color = Model)) +
    geom_line() +
    geom_point() + 
    theme_bw() + 
    labs(title = g_title, 
         x = x_name, y = y_lab) + 
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.title.align = 0.5 , 
          axis.text.x = element_text(hjust = 1)) +
    scale_color_manual(values=col_list[1:K])
  return(p)
}

#Plot a 2D grid of weights
#--hs_wts_grid = wts obtained by applying the get_wts function to the hs_basis_grid
#--v1_ind and v2_ind = index of the predictor in the hsp_basis (note the intercept)
#--k = the model number, i.e. the column number of the wts df 
plot_wts_2d = function(hs_wts_grid, hs_basis_grid, v1_ind, v2_ind, k){
  
  #Get plot data
  data = data.frame(hs_wts_grid[,k], hs_basis_grid[,c(v1_ind, v2_ind)])
  colnames(data) = c('Weights', 'v1', 'v2')
  
  #Get names
  x_name = colnames(hs_basis_grid)[v1_ind]
  y_name = colnames(hs_basis_grid)[v2_ind]
  ptitle = paste('Model',k, 'by', x_name, '&', y_name)
  
  #Plot the data
  p = data %>% ggplot(aes(x = v1, y = v2)) + 
    geom_point(aes(color = Weights), size = 5) + 
    labs(title = ptitle, x = x_name, y = y_name) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    scale_color_gradientn(colors = c("navyblue", "darkmagenta", "darkorange1"),
                          breaks = c(0.25, 0.5, 0.75), 
                          limits = c(0,1))
    
    #scale_color_gradient2(low="midnightblue", mid = 'orchid', high="purple4", midpoint = 0.5, limits = c(0,1))
  return(p)
}

  
