

`multistageoptimal.grid`<-function(N.upper, N.lower=rep(1,length(N.upper)), corr, num.grid=11, Budget, CostC=1, CostTv=rep(1,length(N.upper)),N.fs, detail=FALSE,alg=GenzBretz())

{
  

if (length(N.upper)>4 || length(N.upper)< 2)
   {
    stop( "dimension should be between 2 and 4")
   }

# for N.upper=4,3,2, we prepare for the loops
# comments only available for the case of N.upper=4


if (length(N.upper)==4)
{
   
 # grid for 1st dimension
  z=array(0,rep(num.grid,4))
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid)
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
          
     
     for (i in 1:loop.time)
     {
            # grid for 2nd dimension
            
           if ((N.upper[2]-N.lower[2]+1)<num.grid)
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {
           
            # grid for 3rd dimension 
            
                  if ((N.upper[3]-N.lower[3]+1)<num.grid)
                   {
                      loop.time=N.upper[3]-N.lower[3]+1
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }else{
                      loop.time=num.grid
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }
           
                for (k in 1:loop.time)
                {
                      # grid for 4th dimension
                    
                      if ((N.upper[4]-N.lower[4]+1)<num.grid)
                      {
                        loop.time=N.upper[4]-N.lower[2]+1
                        xNfour <- seq.int(N.lower[4], N.upper[4], length.out=loop.time)
                      }else{
                        loop.time=num.grid
                        xNfour <- seq.int(N.lower[4], N.upper[4], length.out=loop.time)
                       }
           
                
                
                    for (l in 1:loop.time)
                    {
                       
                      # for cost < budget, and xN within the optimization area, we calculate the gain
                       
                           xN=c(xNone[i],xNtwo[j],xNthree[k],xNfour[l])
      
                           cost =  xN[1]*CostC+sum(CostTv*xN)
  
                         if (cost <= Budget && all(xN<=N.upper)&& all(xN>N.lower))
                          {
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] < xN.matrix[,2]))
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corx=corr[-1,-1], alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                     # the gain is saved in a 4 dimensional table
                             z[i,j,k,l]=output    
                    }
                          
                      
                          
                          
                    }
     
                }
     }
     
                              
result=max(z)
  
result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

sample.size=c(xNone[location[1]],xNtwo[location[2]],xNthree[location[3]],xNfour[location[4]])

} else if (length(N.upper)==3)
{

 
  z=array(0,rep(num.grid,3))
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid)
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
          
     
     for (i in 1:loop.time)
     {
           
           if ((N.upper[2]-N.lower[2]+1)<num.grid)
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {
                  if ((N.upper[3]-N.lower[3]+1)<num.grid)
                   {
                      loop.time=N.upper[3]-N.lower[3]+1
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }else{
                      loop.time=num.grid
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }
           
                for (k in 1:loop.time)
                {
  
                           xN=c(xNone[i],xNtwo[j],xNthree[k])
      
                           cost =  xN[1]*CostC+sum(CostTv*xN)
  
                         if (cost <= Budget && all(xN<=N.upper)&& all(xN>N.lower))
                          {
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] < xN.matrix[,2]))
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corx=corr[-1,-1], alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                            
                             
                             z[i,j,k]=output    
                         
                          
                    }
     
                }
     }
result=max(z)
  
result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

sample.size=c(xNone[location[1]],xNtwo[location[2]],xNthree[location[3]])


} else if (length(N.upper)==2)
{


  z=array(0,rep(num.grid,3))
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid)
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
          
     
     for (i in 1:loop.time)
     {
           
           if ((N.upper[2]-N.lower[2]+1)<num.grid)
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {
  
                           xN=c(xNone[i],xNtwo[j])
      
                           cost =  xN[1]*CostC+sum(CostTv*xN)
  
                         if (cost <= Budget && all(xN<=N.upper)&& all(xN>N.lower))
                          {
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] < xN.matrix[,2]))
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corx=corr[-1,-1], alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                            
                             
                             z[i,j]=output    
                         
                          
                    
     
                }
     }
result=max(z)
  
result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

sample.size=c(xNone[location[1]],xNtwo[location[2]])


} else
{
}

# format the table with row and col names


result.table=c(sample.size, result)
 
names(result.table)=c(rownames(sample.size, do.NULL = FALSE, prefix = "N"),"value")

result.table
   
}
