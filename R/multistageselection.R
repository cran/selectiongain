`multistageselection` <-
function(k,corr,alphaofx,sumdimofxandy,alg= GenzBretz())
{


 
dim=sumdimofxandy  
# in one stage selection case dim = 2

for (i in 1:dim)
{
   if (is.infinite(k[i])==TRUE)
   {
     k[i]=-100
   }
}


if (dim<2)
{
stop("dimension of k must bigger than 1, otherwise this is not one-stage selection")

}else if(dim==2)
{
gainresult=corr[1,2]*dnorm(k[2])/alphaofx


}else
{

A=array(0,c(dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
          if(i!=j)
         {
              A[i,j]= (k[j]-corr[i,j]*k[i])/ (1-corr[i,j]^2)^0.5
          }
    }
}


part.corr=array(1,c(dim,dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(i!=j && q!=j && i!=q)

            { part.corr[i,j,q]= (corr[i,j]-corr[i,q]*corr[j,q])/ ((1-corr[i,q]^2)^0.5 * (1-corr[j,q]^2)^0.5)
             }
          }

      }
}





j3q<-function (q,A,part.corr,dim)

{      
    
          
    lower=A[q,-q]
    corr= part.corr[-q,-q,q]


    output=pmvnorm(lower = lower, upper = rep(Inf,c(dim-1)), mean = rep(0, length(lower)), 
    corr = corr, sigma = NULL, algorithm =  alg) 
    output
 }

calculatx1<-function(A,part.corr,dim,corr,k,alpha3)
{ 
   output=0
   i=1
   for (i in 1 : dim)
 {
   output= output+ corr[1,i]*dnorm(k[i])*j3q(i,A,part.corr,dim)/alpha3
 }
   


  output
}
  

gainresult<-calculatx1(A=A,part.corr=part.corr,dim=dim,corr=corr,k=k,alpha3=alphaofx)

}
gainresult
}

