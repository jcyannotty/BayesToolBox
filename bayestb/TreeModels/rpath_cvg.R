#------------------------------------------------
# Random Path: Conditional  Variogram and Covariance 
#------------------------------------------------
# Load functions
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/TreeModels/tree_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/TreeModels/rpath_prior_functions.R")

#------------------------------------------------
# Covariance conditional on a tree and gamma
cov_tree = function(tree,x1,x2,xbnds,tau2,gam){
  # Convert to matrix if not a matrix
  if(!is.matrix(x1)){x1 = matrix(x1, nrow = 1)}
  if(!is.matrix(x2)){x2 = matrix(x2, nrow = 1)}
  
  # Get N
  N = nrow(x1)
  
  # Get bots
  bnv = getbots(tree)
  
  # Get path to root for each bn in bnv
  paths = list()
  for(i in 1:length(bnv)){
    paths[[i]] = pathtoroot(tree, bnv[i])  
  }
  
  # Get lower and upper bounds
  xl = xbnds[,1]; xu = xbnds[,2]
  bnds = list()
  for(i in 1:length(bnv)){
    bnds[[i]] = get_bounds(tree,paths[[i]],xl,xu)
  }
  
  # Get bottom nodes per x and map nid to bn number
  bart_tnode1 = c()
  bart_tnode2 = c()
  for(i in 1:N){
    b1 = get_xtnode(tree,x1[i,])
    b2 = get_xtnode(tree,x2[i,])
    bart_tnode1[i] = which(bnv == b1)
    bart_tnode2[i] = which(bnv == b2)
  }
  # Sample phix
  phix1 = matrix(0,nrow = N, ncol = length(bnv))
  phix2 = matrix(0,nrow = N, ncol = length(bnv))
  
  for(i in 1:N){
    for(b in 1:length(bnv)){
      phix1[i,b] = rpath_phix(tree, paths[[b]], x1[i,], 
                             bnds[[b]]$lbvec, bnds[[b]]$ubvec, gam=gam, q = 4)
      phix2[i,b] = rpath_phix(tree, paths[[b]], x2[i,], 
                              bnds[[b]]$lbvec, bnds[[b]]$ubvec, gam=gam, q = 4)
    }
  }

  # Compute covariance
  cv = tau2*rowSums(phix1*phix2)    
  return(cv)
}


# Semi-Variogram for mean function conditional on a single tree and gamma
svg_tree = function(tree, gam, tau2, xbnds, hgrid, N, sig2 = 0){
  # dims
  p = nrow(xbnds)
  nh = length(hgrid)

  # Generate x points
  x1 = matrix(0, nrow = N, ncol = p)
  for(i in 1:p){
    x1[,i] = runif(N, xbnds[i,1], xbnds[i,2])
  }
  
  # Compute covariance & variogram
  vgf = matrix(0, nrow = N, ncol = nh)
  x2 = x1*0
  for(j in 1:nh){
    # Get x+h points
    for(i in 1:p){
      x2[,i] = x1[,i] + hgrid[j]/sqrt(p) 
    }
    cv = cov_tree(tree, x1, x2, xbnds, tau2, gam)
    vgf[,j] = tau2 - cv + sig2
    cat("VG fit: ", round(j/nh,4),"\r")
  }
  
  # Compute semi-varioogram 
  out = list(vgf = vgf, vgfhat = apply(vgf,2,mean))
  return(out)
}

# Get phix for a tree
phix_tree = function(tree,x,xbnds,gam){
  # Cast to matrix 
  if(!is.matrix(x)){
    x = matrix(x, ncol = 1)
  }
  
  # Number of rows in x
  N = nrow(x)
  
  # Get bots
  bnv = getbots(tree)
  
  # Get path to root for each bn in bnv
  paths = list()
  for(i in 1:length(bnv)){
    paths[[i]] = pathtoroot(tree, bnv[i])  
  }
  
  # Get lower and upper bounds
  xl = xbnds[,1]; xu = xbnds[,2]
  bnds = list()
  for(i in 1:length(bnv)){
    bnds[[i]] = get_bounds(tree,paths[[i]],xl,xu)
  }
  
  # Get bottom nodes per x and map nid to bn number
  bart_tnode = c()
  for(i in 1:N){
    b = get_xtnode(tree,x[i,])
    bart_tnode[i] = which(bnv == b)
  }
  # Sample phix
  phix = matrix(0,nrow = N, ncol = length(bnv))
  for(i in 1:N){
    for(b in 1:length(bnv)){
      phix[i,b] = rpath_phix(tree, paths[[b]], x[i,], 
                              bnds[[b]]$lbvec, bnds[[b]]$ubvec, gam=gam, q = 4)
    }
  }
  return(phix)
}
