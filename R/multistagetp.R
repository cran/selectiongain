`multistagetp` <-
function(alpha,corx,alg=GenzBretz())
{
   corr=corx
    dim=length(alpha)
    
   # if alpha =1, there will be a problem, the user should reduce alpha slightly
  
    for (i in 1:dim)
    {
        if (alpha[i]==1)
        {
           alpha[i]=alpha[i]*0.9999
        }
    
    }
    
 # this is for the root search function   
    qnormbar<-function(kwanted,koutput,corr,howmanyk,alpha, algorithm = alg)
   {   
        corr2dim=corr[1:howmanyk,1:howmanyk] 
        
        tempofk1=koutput[1:c(howmanyk-1)]
          
        alphatotal=1
        for (i in 1:howmanyk)
        {
           alphatotal=alphatotal*alpha[i]
        }

        pmvnorm(lower = c(tempofk1,kwanted), upper = rep(Inf,howmanyk), mean = rep(0,howmanyk),
                   corr = corr2dim, algorithm = alg) - alphatotal   
   }
     
 # repeat the root search for dim dimensions    
     
    koutput=0
    for (i in  1: dim ) 
   {
    
         if (i == 1)     
         {
           koutput[1]= qnorm(alpha[1],lower.tail = FALSE)[[1]][1]
  
         }else
         {   
              index1=qnormbar(kwanted=100,koutput,corr,i,alpha, algorithm = alg)
              index2=qnormbar(kwanted=-100,koutput,corr,i,alpha, algorithm = alg)
              
              #          kwanted should be -100<=  kwanted <= 100
              
             if (index1>0)
            {
                koutput[i]=100
            }else if(index2<0)
            {
               koutput[i]=-100
            }          
            else
            {
             koutput[i]= uniroot(qnormbar,interval=c(-100,100),koutput=koutput,corr,i,alpha)[[1]][1]  
            } 
         }
         if (is.infinite(koutput[i])==TRUE)
         {
           koutput[i]=-100
           # this is for the special case that koutput[1] is INF 
         }
     }

koutput
}

