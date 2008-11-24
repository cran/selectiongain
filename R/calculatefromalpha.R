`calculatefromalpha` <-
function(alpha,dim,corr)
{

    qnormbar<-function(kwanted,koutput,corr,howmanyk,alpha)
   {   
        corr2dim=corr[1:howmanyk,1:howmanyk] 
        
        tempofk1=koutput[1:c(howmanyk-1)]
          
        alphatotal=1
        for (i in 1:howmanyk)
        {
           alphatotal=alphatotal*alpha[i]
        }

        pmvnorm(lower = c(tempofk1,kwanted), upper = rep(Inf,howmanyk), mean = rep(0,howmanyk),
                   corr = corr2dim) - alphatotal   
   }
     
    koutput=0
    for (i in  1: dim )
       {
    
         if (i == 1)     
         {
           koutput[1]= qnorm(alpha[1],lower.tail = FALSE)[[1]][1]
  
         }else
         {

           koutput[i]= uniroot(qnormbar,interval=c(-100,100),koutput=koutput,corr,i,alpha)[[1]][1]  
         }
     }

koutput
}

