`multistagegain.each` <-
function(k,corr,alphaofx,sum.dim,alg= GenzBretz())
{

   gain.array = array (0,c(sum.dim-1))
   
   for (i in 2:c(sum.dim))
   
   {
     alphaofx = pmvnorm(lower=k[1:i],corr=corr[1:i,1:i])
     gain.array[i-1]= multistagegain(k[1:i],corr[1:i,1:i],alphaofx,sum.dim=i,alg= GenzBretz())
   
   }
   
    
   for (i in 2:c(sum.dim-1))
   
   {
     
     gain.array[i]=gain.array[i]-gain.array[i-1] 
   
   }
   
   
gain.array

}

