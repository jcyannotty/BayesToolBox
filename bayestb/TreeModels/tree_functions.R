#------------------------------------------------
# Bayesian Tree Generation
# Desc: Object-Oriented like tree functions

#------------------------------------------------
# Tree Functions
#------------------------------------------------
# Initialize a node
init_node = function(nid = 1){
  out = list(nid = nid, v = 0, c = 0, l = -1, r = -1, p = -1)
  return(out)
}


# Initialize a tree
init_tree = function(){
  root = init_node(1)
  tree = list("nid1" = root)
  return(tree)
}


# Birth step at the inode object
birth = function(inode,v0,c0){
  # Update internal node
  inode$v = v0
  inode$c = c0
  inode$l = inode$nid*2
  inode$r = inode$nid*2 + 1
  
  # Create new terminal nodes
  ltnode =  init_node(inode$l)
  rtnode =  init_node(inode$r)
  ltnode$p = inode$nid
  rtnode$p = inode$nid
  
  # Return everything as a big list
  out = list(inode = inode, ltnode = ltnode, rtnode = rtnode)
  return(out)
}


# Update tree after birth
append_tree = function(tree, new_nodes){
  for(i in 1:length(new_nodes)){
    new_id = paste0("nid",new_nodes[[i]]$nid)
    tree[[new_id]] = new_nodes[[i]]
  }
  return(tree)
}


# Prune the tree after death
prune_tree = function(tree, rm_nodes){
  for(i in 1:length(new_nodes)){
    rm_id = paste0("nid",rm_nodes[[i]]$nid)
    tree = tree[names(tree) %in% rm_id == FALSE]
  }
  return(tree)
}


# Get the terminal node for a given x
get_xtnode = function(tree,xvec){
  cnode = tree$nid1
  tnode = 0
  if(length(xvec)==1){xvec = as.vector(xvec)}
  while(tnode == 0){
    if(cnode$l == -1){
      tnode = cnode$nid
    }else{
      v0 = cnode$v; c0 = cnode$c
      if(xvec[v0] < c0){
        new_nid = paste0("nid",cnode$l) 
      }else{
        new_nid = paste0("nid",cnode$r)
      }
      cnode = tree[[new_nid]]
    }  
  }
  return(tnode)
} 


# Get the terminal nodes of the tree
getbots = function(tree){
  bots = c()
  for(i in 1:length(tree)){
    if(tree[[i]]$l == -1){
      bots = append(bots,tree[[i]]$nid)
    }
  }
  return(bots)
}


# Get the path from a terminal node to the root
pathtoroot = function(tree, bnid){
  path = c(bnid)
  cnid = bnid
  while(cnid > 1){
    nname = paste0("nid",cnid)
    cnid = tree[[nname]]$p
    path = append(path, cnid)
  } 
  return(path)
}


# Get the marginal bounds at each internal node on a path
get_bounds = function(tree, path, xmin, xmax){
  lbvec = 0; ubvec = 0
  if(length(xmax) == 1){xmax = as.vector(xmax)}
  if(length(xmin) == 1){xmin = as.vector(xmin)}
  path = rev(path)
  q = length(path)
  if(q>1){
    for(i in 1:(q-1)){
      # get node info
      nname = paste0("nid",path[i])
      v0 = tree[[nname]]$v
      c0 = tree[[nname]]$c
      
      # Update the lbvec and ubvec
      lbvec[i] = xmin[v0]
      ubvec[i] = xmax[v0]
      
      # Update bounds
      if(path[i+1]%%2==0){
        # Left move, update upper bound
        xmax[v0] = c0
      }else{
        # Left move, update lower bound
        xmin[v0] = c0
      }
    }
  }
  out = list(lbvec = lbvec, ubvec = ubvec) 
  return(out)
}


# Get the cuts used in the tree (1-dim)
get_cuts = function(tree){
  cuts = c()
  vars = c()
  for(i in 1:length(tree)){
    if(tree[[i]]$l != -1){
      cuts = append(cuts,tree[[i]]$c)
      vars = append(vars,tree[[i]]$v)
    }
  }
  out = list(cuts = cuts, vars = vars)
  return(out)
}


# Create cut points for tree building
create_cuts = function(xmin, xmax, ncut){
  if(length(xmax) == 1){xmax = as.vector(xmax)}
  if(length(xmin) == 1){xmin = as.vector(xmin)}
  xcuts = list()
  for(i in 1:length(xmax)){
    xcuts[[i]] = seq(xmin[i],xmax[i],length = ncut)
  }
  return(xcuts)
}

#------------------------------------------------
# Stochastic Tree Generating Process
#------------------------------------------------
# Split probability
psplit = function(a,b,d){
  pr = a*(1+d)^(-b)
  return(pr)
}


# Generate a new tree
generate_tree = function(a,b,xcuts,xmin,xmax){
  tree = init_tree()
  cnode = tree$nid1
  tnodes = c()
  depth_by_id = 2^(1:10)
  p = length(xcuts)
  queue = c(1)
  
  while(length(queue)>0){
    depth = min(which(depth_by_id > cnode$nid)) - 1
    p = psplit(a,b,depth)
    u = runif(1,0,1)
    # Split this node or make terminal
    if(u < p){
      # Get path and bounds along the path
      path = pathtoroot(tree,cnode$nid)
      bnds = get_bounds(tree,path,xmin,xmax)

      # Select a var and then a cut 
      v = sample(1:p,1)

      # Append to the bounds based on the parent's split var
      vmax = xmax[v]; vmin = xmin[v]
      q = length(path)
      if(q>1){
        for(i in 1:(length(path)-1)){
          nname2 = paste0("nid",path[i])
          if(tree[[nname2]]$v == v){
            vmin = bnds$lbvec[i]
            vmax = bnds$ubvec[i]
          }
        }
        # Now include the parent cutpoint c in the bounds
        #nname = paste0("nid",path[q])
        pname = paste0("nid",tree[[nname]]$p)
        # Left vs right node dictates the change of vmax or vmin
        if(tree[[pname]]$l == path[1]){
          vmax = tree[[pname]]$c
        }else{
          vmin = tree[[pname]]$c
        }
      }    
        
      if(q>1){
        cuts = xcuts[[v]][which(xcuts[[v]]>=vmin & xcuts[[v]]<=vmax)]
      }else{
        cuts = xcuts[[v]] # special case for root node
      }  
      c = sample(cuts,1)
      
      # Birth step and append
      out = birth(cnode,v,c)
      tree = append_tree(tree,out)
      
      # Remove this internal node and Add to front of queue (left, right)
      queue = queue[-1]
      queue = c(out$ltnode$nid,out$rtnode$nid,queue)
    }else{
      # Term node, add to the list and pop from the queue
      tnodes = append(tnodes,cnode$nid)
      queue = queue[-1]
    }
    
    # Update the cnode based on the queue
    if(length(queue)>0){
      nname = paste0("nid",queue[1])
      cnode = tree[[nname]]
    }
  }
  return(tree)
}


#tree1 = generate_tree(0.95,2,xcuts,0,8)
#tree2 = generate_tree(0.95,2,xcuts,0,8)
#tree3 = generate_tree(0.95,0.5,xcuts,0,8)
