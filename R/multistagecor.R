# mutistagecor.R author: xuefei mi, first version 17-04-2013, for selectiongain package v2.0.6

# mutistagecor.R author: xuefei mi, second version 11-08-2013, for selectiongain package since v2.0.18

# the package are totally rewritten


`multistagecor` <-function(maseff=NA,VGCAandE=c(0,0,0,0,0),VSCA=c(0,0,0,0),T, L,Rep,index=FALSE,covscatype=c("LonginII"),detail=FALSE)
{

# we had two optition of VGCA,  either VGCA or VGCAandE  

# funtion for calculating cov without markers

covwithoutmas  <-function(V=c(0,0,0,0,0),VSCA=c(0,0,0,0),T, L,Rep,index=FALSE,covscatype=c("LonginII"))
{
    
  VarianceType=V[1]
  
  if (VarianceType=="VC2" )
  {
    Vg=1
    Vgl=0.5
    Vgy=0.5
    Vgly=1
    Ve=2
  }
  else if (VarianceType=="VC1")
  {
    Vg=1
    Vgl=0.25
    Vgy=0.25
    Vgly=0.5
    Ve=1
  }
  else if (VarianceType=="VC3")
  {
    Vg=1
    Vgl=1
    Vgy=1
    Vgly=2
    Ve=4
  }
 
  else  if (length(VGCAandE)==5& all(is.numeric(VGCAandE)))
  { 
    Vg=V[1]
    Vgl=V[2]
    Vgy=V[3]
    Vgly=V[4]
    Ve=V[5]
  }else
  {
   stop("the value of VGCAandE has to be numeric with length 5, or VC1,VC2,VC3 as defined in Longin's paper")
  }
 		
 
 
  dim = length(L)+1
  
  if (length(L) != length(Rep))
  {
    warning("Location and Replicates do not have same length", call. = FALSE)
  }
  
   if (length(L) != length(T))
  {
    warning("Location and testers do not have same length", call. = FALSE)
  }
   
  if (length(L) ==1 && index==TRUE)
  {
    warning("selection index need n >1", call. = FALSE)
  }
  

  # preparation for the VSCA
  
  
  VSCAType=VSCA[1]
  
  if (VSCAType=="VC2.2" )
  {
    Vs=0.5
    Vsl=0.25
    Vsy=0.25
    Vsly=0.5

  }
  else if (VSCAType=="VC2.1")
  {
    Vg=0.5
    Vgl=0.125
    Vgy=0.125
    Vgly=0.25

  }else if(all(is.numeric(VSCA))& length (VSCA)==4)
  {
    Vs=VSCA[1]
    Vsl=VSCA[2]
    Vsy=VSCA[3]
    Vsly=VSCA[4]

   }else
   {
     stop("the value of VSCA has to be numeric with length 4, or VC2.1,VC2.2, as defined in Longin's paper")
   }
  
 
# calculate the covariance  

if (covscatype=="LonginII")
 { 
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1]    + Ve/L[i-1]/Rep[i-1]/T[i-1]
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1],T[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1],T[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covscatype=="Heffner")
 {
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
	LT=L*T
	LTR=L*T*Rep
	
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy/(i-1)+(Vgl+Vgly)/sum(L[1:(i-1)])+  (Vs + Vsy)/sum(T[1:(i-1)]) + (Vsl+Vsly)/sum(LT[1:(i-1)])   + Ve/sum(LTR[1:(i-1)])
  
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {    
        cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+Vs/max(T[i-1],T[j-1])+
        Vsl/max(T[i-1],T[j-1])/max(L[i-1],L[j-1])
        cov[j,i]=cov[i,j]
      }
     }
  } 
 }
  
# calculate the optimal selection index = G^-1 /P


cov

}

# main function begins

        dim = length(L)+1 
        Vg=1

	if (all(is.numeric(VGCAandE)))
	{ 
	   if (length(VGCAandE)==5)
           { 
              Vg=VGCAandE[1] 
	   }else
           {
             stop("the value of VGCAandE has to be numeric with length 5, or VC1,VC2,VC3 as defined in Longin's paper")
           }
 		
        }
        
        if (length(L)!=length(Rep) | length(T)!=length(L))
        {
            stop("T, L and Rep must have the same length")   
        }
        
        
        if (is.na(maseff))
	 {
              tempb="empty"
	      cov=covwithoutmas(V=VGCAandE,VSCA=VSCA,T=T,L=L,Rep=Rep,index=index,
              covscatype=covscatype)
              if (index==TRUE)
             {                
                P=cov[-1,-1]
                G=matrix(rep(1,(dim-1)^2), nrow = dim-1, ncol=dim-1, byrow=TRUE) * Vg
                for (i in 2:c(dim-1))
                {
                   tempb=matrix(rep(1/dim,i),nrow=1,ncol=i,byrow=TRUE)
                   tempp=P[1:i,1:i]
                   tempg=G[1:i,1:i]
                   tempb= t((solve(tempp)%*% tempg )  %*% t( tempb))
                   tempb=tempb / sum(tempb)
                   P[i,i]=  tempb %*% tempp %*% t(tempb)
                   P[1:(i-1),i]= tempb %*% tempp[,1:(i-1)]
                   P[i,1:(i-1)]=t(P[1:(i-1),i])
                }

# calculate the covariance and give output
                   cov[2:dim,2:dim]=P
              }else if (index!=TRUE)
							{
							 
							}
	
                   output= list(cov2cor(cov),tempb,cov)
	   } else if(!is.na(maseff))
       {
             if (length(maseff)!=1)
             {
                 stop("if maseff is not NA the then the length of it has to be 1.")   
             }
        
             
             tempb="empty"
             if (index==TRUE)
             {
                warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
             }else if (index!=TRUE)
             {
                cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,T=T[-1], L=L[-1],Rep=Rep[-1],
								covscatype=covscatype)
               	dim=dim-1
                corp=cov
							
	              cormas= diag(dim+1)
		            cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
		            cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
		            cormas[3:c(dim+1),1]=corp[2:c(dim),1]
	            	cormas[1,2]=maseff^2*Vg
		            cormas[2,1]=maseff^2*Vg
	            	cormas[2,2]=maseff^2*Vg
	              cormas[1,1]=Vg
		
	          	for (i in 3:c(dim+1))
	      	   {
		           cormas[2,i]=maseff^2*Vg
	 	           cormas[i,2]=cormas[2,i]
	           } 
              


             }
             output= list(cov2cor(cormas),tempb,cormas)
         }
if (detail==TRUE)
{
output
}else (detail!=TRUE)
{
output[[1]]
}

}












