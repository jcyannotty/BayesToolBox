#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Bayesian Regression Code
library(MASS)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set priors
#Sigma2 Inverse Gamma Sampler
sig2_sample = function(X, y, a_shape, b_scale, beta_vec, mu_vec, V){
  p = ncol(X) #Number of predictors
  n = nrow(X) #Sample size in design
  y = as.matrix(y, ncol = 1) #cast to a matrix
  
  #Get shape and scale (scale in inv gamma is the rate in gamma)
  V_inv = chol2inv(chol(V))
  sig2_shape = 0.5*(n + p + 2*a_shape)
  sig2_scale = 0.5*(t(y - X%*%beta_vec) %*% (y - X%*%beta_vec) + 
                      t(beta_vec - mu_vec)%*%V_inv%*%(beta_vec - mu_vec) + 2*b_scale)
  sigma2_new = 1/rgamma(n = 1, shape = sig2_shape, rate = sig2_scale)
  return(sigma2_new)
}

#Beta Multivariate Normal Sampler
beta_sample = function(X, y, sig2, mu_vec, V){
  #Get inverses
  V_inv = chol2inv(chol(V))
  Sig_post = chol2inv(chol(t(X)%*%X  + V_inv))
  
  #Get sum of precision weighted means
  b = t(X)%*%y + V_inv%*%mu_vec
  
  #Get posterior mean and variance
  beta_mean = Sig_post%*%b
  beta_var = Sig_post*sig2
  
  beta_new = mvrnorm(n = 1, mu = beta_mean, Sigma = beta_var)  
  return(beta_new)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Construct function for Gibb's Sampling 
bayes_reg = function(X, y, mu_vec, V, a_shape, b_scale, N, N_burn){
  #Initialize posterior draws
  beta_post = matrix(0, nrow = N, ncol = ncol(X))
  sig2_post = 1
  for(i in 2:N){
    #Draw sig2
    sig2_post[i] = sig2_sample(X, y, a_shape, b_scale, beta_vec = beta_post[i-1,], mu_vec, V)
    
    #Draw beta
    beta_post[i,] = beta_sample(X, y, sig2 = sig2_post[i], mu_vec, V)
    cat("Fitting Bayesian Regression Model: ", i/N*100, "Percent Complete\r")
  }
  cat("\nReturning Model...\n")
  #Compile results and remove burn in stage
  post_draws = data.frame(beta_post, sig2_post)
  colnames(post_draws) = c(paste0('b',0:(ncol(X)-1)), 'sig2')
  
  post_draws = post_draws[-c(1:N_burn),]
  return(post_draws)
}

#Get predictions from Bayesian Regression (also return posterior densities, useful for LOO)
#--Read in output from bayes_reg function and prediction design matrix
predict_bayes = function(post_draws, X, y=NULL, return.den = FALSE){
  #Get parameters
  N = nrow(post_draws)
  n = nrow(X)
  beta_mat = as.matrix(post_draws[,1:ncol(X)])
  X = as.matrix(X)
  #Get the posterior predictive distribution and density using the Gibb's Sampler Results
  if(return.den == TRUE){
    post_dist = post_den = matrix(0, nrow = N, ncol = n)
    
    #Get the distribution
    for(i in 1:N){
      #Get beta and mean vectors 
      beta_vec = beta_mat[i,]
      pred_mean = X%*%beta_vec
      
      #Get the posterior predictions and densities
      post_dist[i,] = rnorm(n, mean = pred_mean, sd = post_draws$sig2[i])
      post_den[i,] = dnorm(y, mean = pred_mean, sd = post_draws$sig2[i])
    }
    
    #Get the output list
    out = list(post_dist = post_dist, post_den = post_den)
    return(out)
  }else{
    post_dist = matrix(0, nrow = N, ncol = n)
    #Get the distribution
    for(i in 1:N){
      #Get beta and mean vectors 
      beta_vec = beta_mat[i,]
      pred_mean = X%*%beta_vec
      
      #Get the posterior predictions and densities
      post_dist[i,] = rnorm(n, mean = pred_mean, sd = post_draws$sig2[i])
    }
    out = list(post_dist = post_dist, post_den = NULL)
    return(out)
  }
}

#Posterior Mean function
predict_mean = function(post_draws, X){
  #Get parameters
  N = nrow(post_draws)
  n = nrow(X)
  beta_mat = as.matrix(post_draws[,1:ncol(X)])
  X = as.matrix(X) 
  
  #Get mean function
  post_mean = matrix(0, nrow = N, ncol = n)
  for(i in 1:N){
    #Get beta and mean vectors 
    beta_vec = beta_mat[i,]
    post_mean[i,] = X%*%beta_vec
  }
  return(post_mean)
}

#Get Marginal Likelihood
marginal_lhood = function(y_train, x_design,n_samples,beta_mean, beta_cov, a_shape, b_scale){
  #Compute Marginal Likelihoods numerically
  lhood = 0
  for(i in 1:n_samples){
    sig2 = 1/rgamma(n = 1, shape = a_shape, rate = b_scale)
    beta = mvrnorm(n = 1, mu = beta_mean, Sigma = beta_cov)
    lhood[i] = prod(dnorm(y_train, mean = x_design%*%beta, sd=sqrt(sig2)))  
  }
  mean_lhood = mean(lhood)
  return(mean_lhood)
}

#Get marginal likelihood for models with only error variance and deterministic mean
marginal_lhood_sig2 = function(y_train,f_mean,n_samples,a_shape, b_scale){
  #Compute Marginal Likelihoods numerically
  lhood = 0
  for(i in 1:n_samples){
    sig2 = 1/rgamma(n = 1, shape = a_shape, rate = b_scale)
    lhood[i] = prod(dnorm(y_train, mean = f_mean, sd=sqrt(sig2)))  
  }
  mean_lhood = mean(lhood)
  return(mean_lhood)
}


#Get Model Evidence
model_evidence = function(lhood_mean_vec, prior_model_vec){
  num = lhood_mean_vec*prior_model_vec
  model_ev = num/sum(num)
  return(model_ev)
}
