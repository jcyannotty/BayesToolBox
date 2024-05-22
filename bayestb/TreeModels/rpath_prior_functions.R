#------------------------------------------------
# Random Path - Tree Functions
# Desc: Works with Object-Oriented like tree functions
#------------------------------------------------
# rpath phix 
rpath_phix = function(tree, path, xvec, xminvec, xmaxvec, q, gam=0.5){
  log_phi = 0
  if(length(xvec)==1){xvec = as.vector(xvec)}
  path = rev(path)
  if(length(path) > 1){
    for(i in 1:(length(path)-1)){
      nname = paste0("nid",path[i])
      v0 = tree[[nname]]$v
      c0 = tree[[nname]]$c
      lb = xminvec[i]
      ub = xmaxvec[i]
      if(path[i+1]%%2==0){
        r0 = FALSE
      }else{
        r0 = TRUE
      } 
      prb = rpath_split(xvec[v0],c0,xmin=lb,xmax=ub,r=r0,gam=gam,q=q)
      if(is.na(prb)){print(c0)}
      log_phi = log_phi + log(prb)
    }
    phi = exp(log_phi)
  }else{
    phi = 1  
  }
  return(phi)   
}


# Random path splits
rpath_split = function(x,c,xmin,xmax,r=TRUE,gam=0.5,q=4){
  # Get half points and distances
  a = c - gam*(c-xmin)
  b = c + gam*(xmax-c)
  d1 = (c-x)/(c-a)
  d2 = (x-c)/(b-c)
  
  # Get psi
  if(x<c){
    psi = 0.5*max(1-d1,0)^q
  }else{
    psi = 1 - 0.5*max(1-d2,0)^q
  }  
  # Get the probability of the right/left move (specified by r = TRUE/FALSE)
  if(r){
    prob = psi   
  }else{
    prob = 1 - psi
  }
  return(prob)
}


# Random Assignments
rpath_drawz = function(phix){
  B = ncol(phix)
  n = nrow(phix)
  z = apply(phix, 1, function(x) sample(1:B, size = 1, prob = x))
  z_ind = matrix(0,ncol = B, nrow = n)
  for(i in 1:n){
    z_ind[i,z[i]] = 1
  }
  return(z_ind)
}


# Disagreement Prob - prob of falling in a diff partition than the deterministic split
# bn = botton node umber 1,2,...,B (not the tree node id!!)
disagree_prob = function(bn, phix){
  prob = 1 - phix[bn]
  return(prob)
}

# Covariance function, conditional on tree and gamma
rpath_cov = function(phix1, phix2, tau2){
  out = tau2*sum(phix1*phix2)
  return(out)
}


# Plot phix
plot_phix = function(x, phix, cvec, title = "Phi(x)", x_lab = "x", y_lim = c(0,1)){
  cols = c("red","blue3","green4","orchid3","orange2","lightblue3", "pink2", "black", "gold")
  plot(x, phix[,1], type = 'l', col = 'red', panel.first = {grid(col = 'lightgrey')}, 
       xlab = x_lab, ylab = 'Phi(x)', main = title, lwd = 2, ylim = y_lim)
  for(i in 2:ncol(phix)){
    lines(x, phix[,i], type = 'l', col = cols[i], lwd = 2)
  }
  
  for(i in 1:length(cvec)){
    abline(v = cvec[i], col = 'grey30', lty = 'dashed')  
  }
}



# Plot rpath covariance
plot_rpath_cov = function(x, covm, title = "Cov(x,x')",y_lim = NULL,
                          legendpos = "topright", llabs = NULL, lncol = 2, cols = NULL){
  
  if(is.null(cols)) cols = c("red","blue3","green4","orchid3","orange2","lightblue3", "pink2", "black", "gold")
  if(is.null(y_lim)){y_lim = c(min(covm),max(covm))}
  plot(x, covm[,1], type = 'l', col = cols[1], panel.first = {grid(col = 'lightgrey')}, 
       xlab = 'x', ylab = 'Phi(x)', main = title, lwd = 2, ylim = y_lim)
  for(i in 2:ncol(covm)){
    lines(x, covm[,i], type = 'l', col = cols[i], lwd = 2)
  }
  
  if(is.null(llabs)){llabs = paste0("x'_",1:ncol(covm))}
  legend(legendpos,legend=llabs,col = cols[1:ncol(covm)],ncol = lncol,pch = 16)
}

# Rpath psi(x) derivative
rpath_ddx_psix = function(x,c,xmin,xmax,r=TRUE,gam=0.5,q=4){
  a = c - gam*(c-xmin)
  b = c + gam*(xmax-c)
  d1 = (c-x)/(c-a)
  d2 = (x-c)/(b-c)
  
  # Get psi
  if(x<c & x>a){
    const = q/(gam*(c-xmin))
    func = 0.5*max(1-d1,0)^(q-1)
  }else if(x>c & x<b){
    const = q/(gam*(xmax-c))
    func = 0.5*max(1-d2,0)^(q-1)
  }else if(x==c){
    func = NA; const = 1
  }else{
    func = const = 0
  }
  
  # If left move, derivative is negative
  if(!r){
    const = -const    
  }
  
  ddx = const*func
  return(ddx)  
}


# Rpath phi(x) derivative
rpath_ddx_phix = function(tree, path, xvec, xminvec, xmaxvec, q, gam=0.5){
  log_phi = 0
  ddx_sum = 0
  if(length(xvec)==1){xvec = as.vector(xvec)}
  path = rev(path)
  if(length(path) > 1){
    for(i in 1:(length(path)-1)){
      nname = paste0("nid",path[i])
      v0 = tree[[nname]]$v
      c0 = tree[[nname]]$c
      lb = xminvec[i]
      ub = xmaxvec[i]
      if(path[i+1]%%2==0){
        r0 = FALSE
      }else{
        r0 = TRUE
      } 
      prb = rpath_split(xvec[v0],c0,xmin=lb,xmax=ub,r=r0,gam=gam,q=q)
      ddx = rpath_ddx_psix(xvec[v0],c0,xmin=lb,xmax=ub,r=r0,gam=gam,q=q)
      
      if(prb>0){ddx_scale = ddx/prb}else{ddx_scale = 0}
      if(is.na(prb)){print(c0)}
      
      log_phi = log_phi + log(prb)
      ddx_sum = ddx_sum + ddx_scale
    }
    phi = exp(log_phi)
  }else{
    phi = 1
    ddx_sum = 0
  }
  out = list(phix = phi, ddx = phi*ddx_sum)
  return(out)   
}


# Plot phix
plot_ddx_phix = function(x, phix, cvec, title = "Phi(x)", x_lab = "x", y_lim = c(0,1)){
  cols = c("red","blue3","green4","orchid3","orange2","lightblue3", "pink2", "black", "gold")
  plot(x, phix[,1], col = 'red', panel.first = {grid(col = 'lightgrey')}, 
       xlab = x_lab, ylab = 'd/dx Phi(x)', main = title, pch = 16, cex = 0.7, ylim = y_lim)
  for(i in 2:ncol(phix)){
    points(x, phix[,i], col = cols[i], pch = 16, cex = 0.7)
  }
  
  for(i in 1:length(cvec)){
    abline(v = cvec[i], col = 'grey30', lty = 'dashed')  
  }
}


# Sbart calculations
sbart_phix = function(tree, path, xvec, gam=0.1){
  log_phi = 0
  if(length(xvec)==1){xvec = as.vector(xvec)}
  path = rev(path)
  if(length(path) > 1){
    for(i in 1:(length(path)-1)){
      nname = paste0("nid",path[i])
      v0 = tree[[nname]]$v
      c0 = tree[[nname]]$c
      if(path[i+1]%%2==0){
        r0 = FALSE
      }else{
        r0 = TRUE
      } 
      prb = sbart_split(xvec[v0],c0,r=r0,gam=gam)
      if(is.na(prb)){print(c0)}
      log_phi = log_phi + log(prb)
    }
    phi = exp(log_phi)
  }else{
    phi = 1  
  }
  return(phi)   
}


# Random path splits
sbart_split = function(x,c,r=TRUE,gam){
  x0 = (x-c)/gam
  psi = 1/(1 + exp(-x0))
  if(r){
    prob = psi   
  }else{
    prob = 1 - psi
  }
  return(prob)
}


#------------------------------------------------
# Old functions for vg
#------------------------------------------------
# variance function, conditional on tree and gamma
rpath_var = function(phix, beta, tau2){
  if(length(beta) == 1){beta = rep(beta, length(phix))}
  out = sum((beta^2 + tau2)*phix) - sum(beta%*%t(beta)*phix%*%t(phix))
  return(out)
}


# Variogram - generate x's in p-dimensional space
get_sphere_pts = function(x, nx, h){
  p = length(x)
  m = matrix(0,nrow = nx, ncol = p)
  for(i in 1:nx){
    z = rnorm(p, 0, 1)
    lam = sum(z^2)
    m[i,] = h*z/sqrt(lam)
  }
  return(m)
}


# Variogram (xmat = nxp)
get_variogram_pts = function(xmat,h,nx){
  n = nrow(xmat)
  p = nrow(xmat)
  
  # Get the x points which are h away
  vp = matrix(0,nrow = nx*n, ncol = p)
  for(i in 1:n){
    idx = ((i-1)*nx + 1):(i*nx)
    vp[idx,] = get_sphere_pts(xmat[i,],nx,h)  
  }
  return(vp)
}


# Rpath variogram for response
rpath_variogram_fx = function(phix,phixh,fx,fxh,beta,tau2){
  c12_sum = 0
  c11_sum = 0
  c22_sum = 0
  n = nrow(phix)
  nx = nrow(phixh)/n
  p = ncol(phix)
  # Covariance
  if(nx>1){
    print("WORKING")
    #  for(i in 1:n){
    #idx = ((i-1)*nx+1):(nx*i)
    #c12_temp = matrix(phix[i,],nrow = nx, ncol = p, byrow = TRUE)*phixh[idx,]
    #c12_sum = c12_sum + tau2*sum()
  }else{
    c12_sum = tau2*sum(rowSums(phix*phixh)*rowSums(fx*fxh))
  }
  # Variances
  c11_tot = apply(phix,1,function(x) rpath_var(x,beta,tau2))
  c22_tot = apply(phixh,1,function(x) rpath_var(x,beta,tau2))
  
  # Rescale variances by the fx^2 values
  c11_tot = c11_tot*rowSums(fx^2)
  c22_tot = c22_tot*rowSums(fxh^2)
  
  # Get averages
  c11_hat = mean(c11_tot)
  c22_hat = mean(c22_tot)
  c12_hat = c12_sum/(n*nx)
  
  # Get variogram
  v = c11_hat + c22_hat - 2*c12_hat
  return(v)
}

# Variogram for wts
rpath_variogram_wts = function(phix,phixh,beta,tau2){
  c12_sum = 0
  c11_sum = 0
  c22_sum = 0
  n = nrow(phix)
  nx = nrow(phixh)/n
  p = ncol(phix)
  # Covariance
  if(nx>1){
    for(i in 1:n){
      idx = ((i-1)*nx+1):(nx*i)
      c12_sum = c12_sum + tau2*sum(matrix(phix[i,],nrow = nx, ncol = p, byrow = TRUE)*phixh[idx,])
    }
  }else{
    c12_sum = tau2*sum(phix*phixh)
  }
  # Variances
  c11_tot = apply(phix,1,function(x) rpath_var(x,beta,tau2))
  c22_tot = apply(phixh,1,function(x) rpath_var(x,beta,tau2))
  
  # Get averages
  c11_hat = mean(c11_tot)
  c22_hat = mean(c22_tot)
  c12_hat = c12_sum/(n*nx)
  
  # Get variogram
  v = c11_hat + c22_hat - 2*c12_hat
  return(v)
}


