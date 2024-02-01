#-----------------------------------------------------
# HS Stan Wrappers
#-----------------------------------------------------
# Prediction wrapper
predict_hs_stan = function(fit_hs,x_basis,f_test,tau_mu, tau_dis, tau_sig, qlower = 0.025, 
                           qupper = 0.975){
  # Get parameters
  mu = extract(fiths)$mu
  mu0 = extract(fiths)$mu_0
  sigma = extract(fiths)$sigma
  beta_con = extract(fiths)$beta_con
  sigma_con = extract(fiths)$sigma_con
  beta_dis = extract(fiths)$beta_dis
  sigma_dis = extract(fiths)$sigma_dis
  
  # Get weight draws
  K = ncol(f_test)
  N = nrow(mu)
  n = nrow(f_test)
  wxstar = array(1,dim = c(N,n,K))
  wx = array(1,dim = c(N,n,K))
  px = matrix(0,nrow = N, ncol = n)
  for(i in 1:N){
    for(l in 1:(K-1)){
      b = as.matrix(beta_con[i,l,])*sigma_con[i,l]
      m = (mu0[i] + mu[i,l])*tau_mu
      wxstar[i,,l] = x_basis%*%b + m
    }
    wx[i,,] = exp(wxstar[i,,])/rowSums(exp(wxstar[i,,]))
    px[i,] = rowSums(wx[i,,]*f_test)
    cat("Progress: ",round(100*i/N,4), " Percent \r")
  }
  
  # Summary Statistics
  pred_mean = apply(px, 2, mean)
  pred_ub = apply(px, 2, quantile, qupper)
  pred_lb = apply(px, 2, quantile, qlower)
  
  wts_mean = wts_ub = wts_lb = matrix(0,nrow = n,ncol = K)
  for(l in 1:K){
    wts_mean[,l] = apply(wx[,,l],2,mean)
    wts_ub[,l] = apply(wx[,,l],2,quantile,qupper)
    wts_lb[,l] = apply(wx[,,l],2,quantile,qlower)
  }
  
  out = list(pred_draws = px, pred_mean = pred_mean, pred_ub = pred_ub, pred_lb = pred_lb,
             wts_draws = wx, wts_mean = wts_mean, wts_ub = wts_ub, wts_lb = wts_lb)
  return(out)
}


# Wrapper for fits
fitted_hs_stan = function(wxpost,f,qlower = 0.025,qupper = 0.975){
  N = nrow(wxpost)
  n = nrow(f)
  K = ncol(f)
  px = matrix(0,ncol = n, nrow = N)
  wts_mean = wts_lb = wts_ub = matrix(0,ncol = K, nrow = n)
  for(j in 1:K){
    px = px + t(t(wxpost[,,j])*f[,j]) 
    wts_mean[,j] = apply(wxpost[,,j],2,mean)
    wts_lb[,j] = apply(wxpost[,,j],2,quantile,qlower)
    wts_ub[,j] = apply(wxpost[,,j],2,quantile,qupper)
  }
  pred_mean = apply(px,2,mean)
  pred_lb = apply(px,2,quantile,qlower)
  pred_ub = apply(px,2,quantile,qupper)
  
  out = list(pred_mean = pred_mean, pred_lb = pred_lb, pred_ub = pred_ub,
             wts_mean = wts_mean, wts_lb = wts_lb, wts_ub = wts_ub)
  return(out)
}

# Normalization on inputs
normx = function(x){
  x = as.matrix(x)
  p = ncol(x)
  xt = x
  for(i in 1:p){
    xt[,i] = (x[,i]- mean(x[,i]))/sd(x[,i]) 
  }
  return(xt)
}
