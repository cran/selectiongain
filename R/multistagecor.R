
# mutistagecor.R author: xuefei mi, 17-04-2013, for selectiongain package v2.0.6

# i am trying to add AccForGs=0 into the parameters

# i am trying to add Testers into the parmeters

`multistagecor` <-function(V=c(0,0,0,0,0),L=c(0,0),Rep=c(1,1),index=FALSE,coe)
{

# set variance components

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
 
  else
  {
    Vg=V[1]
    Vgl=V[2]
    Vgy=V[3]
    Vgly=V[4]
    Ve=V[5]
   }
 
  dim = length(L)+1
  
  if (length(L) != length(Rep))
  {
    warning("L and Rep do not have same length", call. = FALSE)
  }
   
  if (length(L) ==1 && index==TRUE)
  {
    warning("selection index need n >1", call. = FALSE)
  }
  
# calculate the covariance  
  
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1]+ Ve/L[i-1]/Rep[i-1]
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }

  tempb="empty"

# calculate the optimal selection index = G^-1 /P

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
  }
  list(cov2cor(cov),tempb)
}

