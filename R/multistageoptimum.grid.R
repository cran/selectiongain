
# package selectiongain

# modified at 28-06-2013, for 1MAS+2PS

`multistageoptimum.grid`<-function(corr, Vg=1, num.grid=NA, width=NA,  Budget, CostProd=rep(0,length(N.upper)), CostTest=rep(1,length(N.upper)),Nf, alg=GenzBretz(),detail=FALSE,fig=FALSE, N.upper, N.lower)
{

# pre-define parameters 

   CostInitial=CostProd
   Vgca=Vg
   CostTv=CostTest  
   percent =0.0001
   N.fs=Nf

# main function begins
  if (length(CostTest)!= dim(corr)[1]-1)
  {
    stop( "dimension of CostTest has to be dim(corr)[1]-1")
  }  
  if (length(CostProd)!= dim(corr)[1]-1)
  {
    stop( "dimension of CostProd has to be dim(corr)[1]-1")
  }  


  if (Budget> sum( N.upper*(CostTest+CostProd)))
  {
    stop("N.upper is too small, try to set sum( N.upper*(CostTest+CostProd))>Budget")
  }
  	

  if (length(N.upper)>5 || length(N.upper)< 2)
  {
    stop( "dimension of N.upper should be between 2 and 5")
  }
   
  if (length(N.upper)!= length(N.lower))
  {
    stop( "dimension of upper and lower limit has to be equal")
  }
 
   if (length(N.upper)!= dim(corr)[1]-1)
  {
    stop( "dimension of upper has to be dim(corr)[1]-1")
  }

  if ( is.na(num.grid[1]) & is.na(width[1]))
  {
    stop( "one of the num.grid or width should not be NA")
  }
 
  if ( !is.na(num.grid[1]) & !is.na(width[1]))
  {
    stop( " only one of the num.grid and width can be used, not both of them")
  }
 
  if ( !is.na(width[1]) & is.na(num.grid[1]))
  {
       if ( length(width) == length(N.upper))
       {
         num.grid=ceiling((N.upper-N.lower)/width)[-1]+1
       }else if(length(width) != length(N.upper))
       {
           stop( "if width is not NA, then length of width must be identical to N.upper") 
       }
  }
     
  if (length(num.grid)==1)
  {
    num.grid=rep(num.grid, length(N.upper)-1)
  }
      

# for N.upper=5,4,3,2, we prepare for the loops
# comments only available for the case of N.upper=4

# the first stage constraint is removed, because it will be come a constraint factor by Budget

  Nup=N.upper[1]
  Nlo=N.lower[1]

  N.upper=N.upper[-1]
  N.lower=N.lower[-1]


  if (length(N.upper)==4)
  {
   
# grid for 1st dimension
     z=array(0,num.grid)
     if ((N.upper[1]-N.lower[1]+1)<num.grid[1])
     {
        loop.time=N.upper[1]-N.lower[1]+1
        xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
     }else
     {
        loop.time=num.grid[1]
        xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
     }
        maxz12=array(0,num.grid[1:2])       
     
     for (i in 1:loop.time)
     {
         # grid for 2nd dimension
            
           if ((N.upper[2]-N.lower[2]+1)<num.grid[2])
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid[2]
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {           
              # grid for 3rd dimension 
            
                  if ((N.upper[3]-N.lower[3]+1)<num.grid[3])
                   {
                      loop.time=N.upper[3]-N.lower[3]+1
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }else{
                      loop.time=num.grid[3]
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }
           
                for (k in 1:loop.time)
                {
                    # grid for 4th dimension
                    
                      if ((N.upper[4]-N.lower[4]+1)<num.grid[4])
                      {
                        loop.time=N.upper[4]-N.lower[2]+1
                        xNfour <- seq.int(N.lower[4], N.upper[4], length.out=loop.time)
                      }else{
                        loop.time=num.grid[4]
                        xNfour <- seq.int(N.lower[4], N.upper[4], length.out=loop.time)
                       }
           
                
                
                    for (l in 1:loop.time)
                    {
                       
                      # for cost < budget, and xN within the optimization area, we calculate the gain
                       
                           xN=c(xNone[i],xNone[i],xNtwo[j],xNthree[k],xNfour[l])
      
                           cost = round( sum(xN*CostInitial)+sum(CostTv*xN) ,3 )
  
                         if (cost <= Budget && all(xN[-1]<=N.upper)&& all(xN[-1]>=N.lower) & cost >= percent*Budget )
                          {
                               xN[1]= floor((Budget-cost)/(CostInitial[1]+CostTv[1])+xNone[i] )
                               
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] <= xN.matrix[,2])& xN[1]<=Nup & xN[1]>=Nlo )
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corr=corr, alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)*Vgca[1]^0.5
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
                             if (detail==TRUE )
                             { 
                               if (i==1 && j==1 && k==1 &&l==1)
                               {
                                 detail.table=c(xN[1],xNone[i],xNtwo[j],xNthree[k],xNfour[l],output)
                               }else
                               {
                                 detail.table=rbind(detail.table,c(xN[1],xNone[i],xNtwo[j],xNthree[k],xNfour[l],output))
                               }
                             }
                                            
                    }
     
                }

            maxz12[i,j]=max(z[i,j,,])
     }

      result=max(z)
  
  #    result.table = array(0,length(N.upper)+1)

      location = which(z==result,arr.ind =TRUE )

    #  sample.size=c(xNone[location[1]],xNtwo[location[2]],xNthree[location[3]],xNfour[location[4]],result)
           
             xNzero=floor((Budget-sum((CostInitial[2:5]+CostTv[2:5])*c(xNone[location[1]],xNtwo[location[2]],xNthree[location[3]],xNfour[location[4]])))/(CostInitial[1]+CostTv[1]))
           
           max.allocation=c(xNzero, xNone[location[1]],xNtwo[location[2]],xNthree[location[3]],xNfour[location[4]],result)

      if (detail==TRUE)
      {
      
       # save(detail.table,file="detail output")
			
        sample.size = rbind(detail.table,max.allocation)
      
      }
   
    }


  }else if (length(N.upper)==3)
{

 
  z=array(0,num.grid)
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid[1])
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid[1]
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
          
         maxz12=array(0,num.grid[1:2])     
  
     for (i in 1:loop.time)
     {
           
           if ((N.upper[2]-N.lower[2]+1)<num.grid[2])
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid[2]
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {
                  if ((N.upper[3]-N.lower[3]+1)<num.grid[3])
                   {
                      loop.time=N.upper[3]-N.lower[3]+1
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }else{
                      loop.time=num.grid[3]
                      xNthree <- seq.int(N.lower[3], N.upper[3], length.out=loop.time)
                    }
           
                for (k in 1:loop.time)
                {
  
                           xN=c(xNone[i],xNone[i],xNtwo[j],xNthree[k])
      
                          # cost =  sum(xN*CostInitial)+sum(CostTv*xN)+xNone[i]
                           
                            cost = round( sum(xN*CostInitial)+sum(CostTv*xN) ,3 )
  
                         if (cost <= Budget && all(xN[-1]<=N.upper)&& all(xN[-1]>=N.lower) & cost >= percent*Budget)
                          {
                                xN[1]= floor( (Budget-cost)/(CostInitial[1]+CostTv[1])+xNone[i])
                               
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] <= xN.matrix[,2])& xN[1]<=Nup & xN[1]>=Nlo)
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corr=corr, alg=alg)     
                           
                               output <- multistagegain(Q=Quantile,corr=corr,alg=alg)*Vgca[1]^0.5
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                            
                             
                             z[i,j,k]=output   
 
                             if (detail==TRUE )
                             { 
                               if (i==1 && j==1 && k==1)
                               {
                                 detail.table=c(xN[1],xNone[i],xNtwo[j],xNthree[k],output)
                               }else
                               {
                                 detail.table=rbind(detail.table,c(xN[1],xNone[i],xNtwo[j],xNthree[k],output))
                               }
                             }
                          
                          
                    }

                 maxz12[i,j]=max(z[i,j,])
     
                }
          
     }
result=max(z)
  
#result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

xNzero=floor((Budget-sum((CostInitial[2:4]+CostTv[2:4])*c(xNone[location[1]],xNtwo[location[2]],xNthree[location[3]])))/(CostInitial[1]+CostTv[1]))

max.allocation=c(xNzero, xNone[location[1]],xNtwo[location[2]],xNthree[location[3]],result)

if (detail==TRUE)
    {
       
      sample.size = rbind(detail.table,max.allocation)
      
    }
   
} else if (length(N.upper)==2)
{


  z=array(0,num.grid)
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid[1])
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid[1]
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
        maxz12=array(0,num.grid[1:2])              
     
     for (i in 1:loop.time)
     {
           
           if ((N.upper[2]-N.lower[2]+1)<num.grid[2])
           {
            loop.time=N.upper[2]-N.lower[2]+1
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }else{
            loop.time=num.grid[2]
            xNtwo <- seq.int(N.lower[2], N.upper[2], length.out=loop.time)
           }
           
           for (j in 1:loop.time)
           {
  
                           xN=c(xNone[i],xNone[i],xNtwo[j])
      
                         #  cost =  sum(xN*CostInitial)+sum(CostTv*xN)
                           
                            cost = round( sum(xN*CostInitial)+sum(CostTv*xN) ,3 )
  
                         if (cost <= Budget && all(xN[-1]<=N.upper)&& all(xN[-1]>=N.lower) & cost >= percent*Budget)
                          {
                               xN[1]= floor((Budget-cost)/(CostInitial[1]+CostTv[1])+xNone[i] )
                               
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] <= xN.matrix[,2])& xN[1]<=Nup & xN[1]>=Nlo)
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corr=corr, alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)*Vgca[1]^0.5
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                            
                             
                             z[i,j]=output  
                             if (detail==TRUE )
                             { 
                               if (i==1&& j==1)
                               {
                                 detail.table=c(xN[1],xNone[i],xNtwo[j],output)
                               }else
                               {
                                 detail.table=rbind(detail.table,c(xN[1],xNone[i],xNtwo[j],output))
                               }
                             }
                          
                              maxz12[i,j]=max(z[i,j])  
  
                         
                          
                    
     
                }
     }
result=max(z)
  
#result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

#sample.size=c(xNone[location[1]],xNtwo[location[2]],result)

xNzero=floor((Budget-sum((CostInitial[2:3]+CostTv[2:3])*c(xNone[location[1]],xNtwo[location[2]])))/(CostInitial[1]+CostTv[1]))

max.allocation=c(xNzero, xNone[location[1]],xNtwo[location[2]],result)


if (detail==TRUE)
    {
       
      sample.size = rbind(detail.table,max.allocation)
      
    }
   


} else if (length(N.upper)==1)
{


  z=array(0,num.grid)
  
  
           if ((N.upper[1]-N.lower[1]+1)<num.grid[1])
           {
            loop.time=N.upper[1]-N.lower[1]+1
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }else{
            loop.time=num.grid[1]
            xNone <- seq.int(N.lower[1], N.upper[1], length.out=loop.time)
           }
  # maxz12=array(c(0,0),rep(loop.time,2))            
     
     for (i in 1:loop.time)
     {
           
           
  
                           xN=c(xNone[i],xNone[i])
      
                          # cost =  sum(xN*CostInitial)+sum(CostTv*xN)
                           
                            cost = round( sum(xN*CostInitial)+sum(CostTv*xN) ,3 )
  
                         if (cost <= Budget && all(xN[-1]<=N.upper)&& all(xN[-1]>=N.lower) & cost >= percent*Budget)
                          {
                               xN[1]= floor((Budget-cost)/(CostInitial[1]+CostTv[1])+xNone[i] )
                               
                               xN.matrix=embed(c(xN,N.fs),2)

                               if (all(xN.matrix[,1] <= xN.matrix[,2])& xN[1]<=Nup & xN[1]>=Nlo)
                             {
                               alpha = xN.matrix[,1]/ xN.matrix[,2]
     
                               Quantile= multistagetp(alpha, corr=corr, alg=alg)     
                           
                               output<- multistagegain(Q=Quantile,corr=corr,alg=alg)*Vgca[1]^0.5
                              }else
                              { 
                               output=0
                              }
                             
                           }else
                           {
                               output=0
                           }
                            
                             
                              z[i]=output  
                             if (detail==TRUE )
                             { 
                               if (i==1)
                               {
                                 detail.table=c(xN[1],xNone[i],output)
                               }else
                               {
                                 detail.table=rbind(detail.table,c(xN[1],xNone[i],output))
                               }
                             }
                          
                             # maxz12[i,j]=max(z[i,j])  
            
     }
result=max(z)
  
#result.table = array(0,length(N.upper)+1)

location = which(z==result,arr.ind =TRUE )

#sample.size=c(xNone[location[1]],xNtwo[location[2]],result)

xNzero=floor((Budget-sum((CostInitial[2]+CostTv[2])*c(xNone[location[1]])))/(CostInitial[1]+CostTv[1]))

max.allocation=c(xNzero, xNone[location[1]],result)
#sample.size=max.allocation

if (detail==TRUE)
    {
      sample.size = rbind(detail.table, max.allocation)      
    }

} 

if (fig==TRUE)
{
    if ((length(N.upper)==1))
    {
     stop("for two dimension with budget constraint, only one degree of freedom, can not make 2-D figure")
    }
    
    
    require(grDevices)
    #op <- par(mfrow = c(2, 2))
    #points(x=greenzx,y=greenzy, col = "green", pch = 20)
    
    N1=xNone
    N2=xNtwo
  # image(N1, N2, maxz12, col=terrain.colors(100))
   
    image(N1, N2, maxz12, col=0,xlab=expression(N["2"],font=3),ylab=expression(N["3"],font=3),font.lab=3)
     
   z0=which(z==0)
    contour(N1, N2, maxz12,add=TRUE,col = "peru",labcex =1,levels =round( seq((min(z[-z0])+(max(z[-z0])-min(z[-z0]))/3), max(z[-z0]), by = (max(z[-z0])-min(z[-z0]))/10),3))
		
    points(x=xNone[location[1]],y=xNtwo[location[2]], col = "red", pch = 20)
#    abline(v=2000, lty = 2,col="black")
 #   abline(a=0,b=1, lty = 2,col="black")
 }



# format the table with col name 
    colword<-c("N1")
    
    for (i in 1:length(N.lower))
    {
         colword<-c(colword,paste("N",i+1,sep=""))
    }
    if (detail==FALSE)
    {
       names(max.allocation)<-c(colword,"gain")
       max.allocation
    }else
    {
       colnames(sample.size)<-c(colword,"gain")
       sample.size
    }
    
   
}