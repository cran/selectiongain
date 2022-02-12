# mutistagecor.R author: xuefei mi, 11-03-2013, for selectiongain package v2.0.2

`multistagegain.each` <-
function(corr,Q, alg=GenzBretz(),Vg=1)
{ 
# pre-set interanl parameters
  Vg=Vg 
  lim.y=-200
  sum.dim=length(Q)+1
  k=c(lim.y,Q)
  gain.array = array (0,c(sum.dim-1))

# main function begins
  
# The output is given as
# (G1(y), G2(y)-G1(y),G3(y)-G2(y), ...), where Gi(y) refers to the total selection
# gain after the first i stages of selection.
   
  for (i in 2:c(sum.dim))   
  {
     gain.array[i-1]= multistagegain(Q=k[2:i],corr[1:i,1:i],alg= alg, Vg=Vg)
  }
  for (i in 2:c(sum.dim-1))
  { 
     gain.array[i]=gain.array[i]-sum(gain.array[(i-1):1])   
  } 
gain.array#*Vg^0.5
}

