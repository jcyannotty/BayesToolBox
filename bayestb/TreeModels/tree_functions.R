#------------------------------------------------
# Bayesian Tree Generation
# Desc: Object-Oriented like tree functions
#------------------------------------------------
init_node = function(nid = 1){
  out = list(nid = nid, v = 0, c = 0, l = -1, r = -1, p = -1)
  return(out)
}

init_tree = function(){
  root = init_node(1)
  tree = list("nid1" = root)
  return(tree)
}

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

append_tree = function(tree, new_nodes){
  for(i in 1:length(new_nodes)){
    new_id = paste0("nid",new_nodes[[i]]$nid)
    tree[[new_id]] = new_nodes[[i]]
  }
  return(tree)
}

prune_tree = function(tree, rm_nodes){
  for(i in 1:length(new_nodes)){
    rm_id = paste0("nid",rm_nodes[[i]]$nid)
    tree = tree[names(tree) %in% rm_id == FALSE]
  }
  return(tree)
}

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


getbots = function(tree){
  bots = c()
  for(i in 1:length(tree)){
    if(tree[[i]]$l == -1){
      bots = append(bots,tree[[i]]$nid)
    }
  }
  return(bots)
}


pathtoroot = function(tree, bnid){
  path = c(bnid)
  cnid = bnid
  while(cnid > 1){
    nname = paste0("nid",cnid)
    cnid = tree[[cnid]]$p
    path = append(path, cnid)
  } 
  return(path)
}


get_bounds = function(tree, path, xmax, xmin){
  lbvec = 0; ubvec = 0
  if(length(xmax) == 1){xmax = as.vector(xmax)}
  if(length(xmin) == 1){xmin = as.vector(xmin)}
  path = rev(path)
  for(i in 1:(length(path)-1)){
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
  out = list(lbvec = lbvec, ubvec = ubvec) 
  return(out)
}

