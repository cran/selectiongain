# mutistageoptimal.nlm.R author: xuefei mi, 28-03-2013, for selectiongain package v2.0.2


`multistageoptimal.nlm`<-function(N.upper, N.lower=rep(N.fs,length(N.upper)), corr, ini.value=c(c(length(N.upper):1)+Budget/(CostC+sum(CostTv)+1)), Budget, CostC=1, CostTv=rep(1,length(N.upper)), N.fs=1, iterlim=100, alg=GenzBretz())
{

# check, if initial value is greater than budget

   cost =  ini.value[1]*CostC+sum(CostTv*ini.value)
   if (cost>Budget)
   {
    stop( "initial cost is higher than the budget, try other values")
   }

# MultiStageOptim() calls multistagegain() to calculate the selection gain for multi-stage selection
# here xN is a variable, and all the other parameters are fixed
# this is just a shell function, which fits the requirement of the optimization function
      
   MultiStageOptim<-function(xN,N.upper, N.lower,corr,Budget,CostC,CostTv,N.fs,alg)       
   {       
     length.xN=length(xN)     
     cost =  xN[1]*CostC+sum(CostTv*xN)  
     if (cost <= Budget && all(xN<=N.upper)&& all(xN>N.lower) && cost)     
     {     
       xN.matrix=embed(c(xN,N.fs),2)

       if (all(xN.matrix[,1] < xN.matrix[,2]))
       {   
         alpha = xN.matrix[,1]/ xN.matrix[,2]
         Quantile= multistagetp(alpha, corx=corr[-1,-1], alg=alg)                               
         output<- -multistagegain(Q=Quantile,corr=corr, alg=alg)
     
         if (cost<=Budget*0.9)
         {
           output=output*0.8
         }             
       }else
         { 
          output=0
         }
     }else
     {
       xN.matrix=embed(c(xN,N.fs),2)

       if (all(xN.matrix[,1] < xN.matrix[,2]))
       {
         alpha = xN.matrix[,1]/ xN.matrix[,2]
         Quantile= multistagetp(alpha, corx=corr[-1,-1], alg=alg)                               
         output<- -multistagegain(Q=Quantile,corr=corr,alg=alg)
       }else
       { 
         output=0
       }
       output = 0
     }
     output
   }
   
     # nlm algorithm
      result<-nlm(MultiStageOptim,p=ini.value,N.upper=N.upper,N.lower=N.lower,corr=corr,Budget=Budget,CostC=CostC,CostTv=CostTv,N.fs=N.fs,alg=alg,hessian = TRUE,iterlim=iterlim)
  
     result.table = array(0,length(N.upper)+1)
     sample.size=result$estimate
     result.table=c(sample.size, -result$minimum)
     names(result.table)=c(rownames(sample.size, do.NULL = FALSE, prefix = "N"),"value")
     result.table 
}

