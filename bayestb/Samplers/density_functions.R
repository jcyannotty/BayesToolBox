#-----------------------------------------------------
# Densities
#-----------------------------------------------------
# scaled inverse chi-squared
dscinvchi2 = function(x,v,lam){
  # Use inv-gamma identity (invg scale = gamma rate)
  num = (lam*v/2)^(v/2)*exp(-v*lam/(2*x))
  den = gamma(v/2)*x^(1+v/2)
  out = num/den
  return(out)
}

# Scaled and shifted t-distribution
dtscaled = function(x,df,mean,scale){
  out = dt((x-mean)/scale,df)/scale
  return(out)
}
