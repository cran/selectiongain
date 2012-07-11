
`multistageoptimal.pb`<-function(N.upper, N.lower, num.grid=11, Budget,  CostC, CostTv, V, L1, Rep, N.fs, detail=FALSE, fig=FALSE, alg)
{

length=num.grid
Location1=L1
NumSelected=N.fs
VarianceType=V
CostDH=CostC

# TwoStageDHoptim() calls TwoStageDH() to calculate the selection gain for two stage selection.
# here xN is a variable, and all the other parameters are fixed.
# this is just a shell function which fits the requirement of the optimization function.
     
  
   TwoStageDHoptim<-function(xN,Budget,CostDH=CostC,Location1,VarianceType,NumSelected,alg)       
{       
 Location2 = (Budget-(CostDH+Location1*CostTv[1])*xN[1])/xN[2]/CostTv[1]
                                      
 output<-  TwoStageDH(Budget=Budget, CostDH=CostDH,Location1=Location1,Location2=Location2, VarianceType=VarianceType,NumSelected=NumSelected,None=xN[1], Ntwo=xN[2],alg=alg)[[1]][1]
 output= 0-output
 output
}

# calculate the selection gain for two stage selection.
                                        
TwoStageDH<-function(Budget=200, CostDH=CostC,Location1=1,Location2=1, VarianceType="VC2",NumSelected=1,None=0, Ntwo=0,messageon=FALSE,index=TRUE,alg)
{

if (Ntwo<=None&&Location1>0&&Location2>0&&Ntwo>0)
{
if (VarianceType=="VC2" )

{
if (messageon==TRUE)
{
 message("VarianceType are 1:0.5:0.5:1:2 (VC2)")
}
 Vg=1
 Vgl=0.5
 Vgy=0.5
 Vgly=1
 Ve=2
 

}else if (VarianceType=="VC1")
{
if (messageon==TRUE)
{
message("VarianceType are 1:0.25:0.25:0.5:1 (VC2)")
}

 Vg=1
 Vgl=0.25
 Vgy=0.25
 Vgly=0.5
 Ve=1

}else if (VarianceType=="VC3")
{
if (messageon==TRUE)
{
message("VarianceType are 1:1:1:2:4 (VC3)")
}
 Vg=1
 Vgl=1
 Vgy=1
 Vgly=2
 Ve=4


}else if (VarianceType=="VC4")
{
if (messageon==TRUE)
{
message("VarianceType are 0.5:0.25:0.25:0.5:2 (VC4)")
}

 Vg=0.5
 Vgl=0.25
 Vgy=0.25
 Vgly=0.5
 Ve=2

}else if (FALSE)
{
 Vg=V[1]
 Vgl=V[2]
 Vgy=V[3]
 Vgly=V[4]
 Ve=V[5]

}

if(Ntwo==0)
{
Ntwo= floor( (Budget-(Location1+CostDH)*None)/(Location2))
Ntwo= max(1, Ntwo)
}


cov= diag(3)*Vg

cov[1,1]=Vg
cov[1,2]=Vg
cov[2,1]=Vg
cov[1,3]=Vg
cov[3,1]=Vg
cov[2,2]=Vg+Vgy+(Vgl+Vgly+Ve/Rep[1])/Location1
cov[3,3]=Vg+Vgy+(Vgl+Vgly+Ve/Rep[2])/Location2

cov[2,3]=Vg+Vgl/Location2
cov[3,2]=cov[2,3]


# for optimal selection index

if (index==TRUE)

{
    b1=0.5
    b2=0.5
    
    P=cov[c(2,3),c(2,3)]
    G=matrix(rep(1,4), nrow = 2, ncol=2, byrow=TRUE)

    G[1,1]=Vg
    G[2,2]=Vg
G[1,2]=min(G[1,1],G[2,2])
G[2,1]=G[1,2]



tempb= solve(P)%*% G %*% c(b1,b2)
b1=tempb[1]/sum(tempb)
b2=tempb[2]/sum(tempb)
    

cov[3,3]=b1^2*cov[2,2]+b2^2*cov[3,3]+2*b1*b2*cov[3,2]

cov[2,3]=b1*cov[2,2]+b2*cov[3,2]
cov[3,2]=cov[2,3]



}

correlation=cov2cor(cov)



alpha1<-Ntwo/None
alpha2<-NumSelected/Ntwo

       k<-multistagetp(alpha=c(alpha1,alpha2),corx=correlation[-1,-1])
       
       
    

outputgain<-multistagegain(Q=k,corr=cov2cor(cov),alg=alg)


t=rep(0,3)


outputgain=outputgain*Vg^0.5 

# Ntwo is not possible > None 


}else
{
outputgain=0

}


list(outputgain)
}

# main part of the function


# divide the interval [N.lower[1],N.upper[1]] and [N.lower[2],N.upper[2]] into n parts 


  xNone <- seq(N.lower[1], N.upper[1], length=length)

  yNtwo <- seq(N.lower[2], N.upper[2], length=length)
  
  z=array(0,c(length,length))
  
  output.table=array(0,c(length^2,5))
  
# calculate the z=gain, and store it,  this is a grid algorithm
# we need it for the contour plot and for the initial value of the nlm 

# the output.table will store all the possible grids, when detail == True, it will be send out.
  
for (i in 1:length)
{
 for (j in 1:length)
 {
  xN=c(xNone[i],yNtwo[j])
  z[i,j] <- - TwoStageDHoptim(c(xNone[i],yNtwo[j]),Budget=Budget,CostDH=CostC, Location1=Location1, VarianceType=VarianceType,NumSelected=NumSelected,alg)
  Location2 = floor((Budget-(CostDH+Location1*Rep[1])*xN[1])/xN[2]/Rep[2])
  output.table[(i-1)*length+j,]=c(Location1,Location2,xNone[i],yNtwo[j],z[i,j])
 
 }
}



  gainz=max(z)
  locationz=which(z==gainz,arr.ind =TRUE )
  
  greenzx=xNone[locationz[1]]
  greenzy=yNtwo[locationz[2]]
  
  if (fig==TRUE)
  {
  contour(xNone, yNtwo, z)
  points(x=greenzx,y=greenzy, col = "green", pch = 20)
  }
  
  output=array(0,c(3,9))
  colnames(output)<-c("NumSelected","Budget","Location1","Location2","N1","N2","Rep1","Rep2","gain")
  rownames(output)<-c("grid","nlm","round-nlm")
   
  Location2 = floor((Budget-(CostDH+Location1*Rep[1])*greenzx)/greenzy/Rep[2])
 
 # this is the result from the grid algorithm
   
  output[1,]=c(NumSelected,Budget,Location1,Location2,greenzx,greenzy,Rep[1],Rep[2],gainz)


# use the nlm algorithm to calculate the maximum
# the initial value is from the grid search algorithm described above

  nlm.f2 <- nlm(TwoStageDHoptim,p=c(greenzx,greenzy),Budget=Budget,CostDH=CostDH,Location1=Location1,VarianceType= VarianceType,NumSelected=NumSelected,alg=alg,hessian = TRUE)
  
  xN=c(0,0)
  
  xN[1]=nlm.f2$estim[1]
  xN[2]=nlm.f2$estim[2]
  Location2 =floor((Budget-(CostDH+Location1*Rep[1])*xN[1])/xN[2]/Rep[2])
 
# this is the result from the nlm algorithm  
  
  output[2,]=c(NumSelected,Budget,Location1,Location2,nlm.f2$estim[1],nlm.f2$estim[2],Rep[1],Rep[2],-nlm.f2$minimum)
  
  xN[1]=round(nlm.f2$estim[1])
  xN[2]=round(nlm.f2$estim[2])
  Location2 = floor((Budget-(CostDH+Location1*Rep[1])*xN[1])/xN[2]/Rep[2])
  
  gainz.round=-TwoStageDHoptim(c(round(nlm.f2$estim[1]),round(nlm.f2$estim[2])),Budget=Budget, CostDH=CostDH,Location1=Location1, VarianceType=VarianceType,NumSelected=NumSelected,alg=alg)

# this is the result from the nlm algorithm with location been rounded to the closest integer 
   
  output[3,]=c(NumSelected,Budget,Location1,Location2,round(nlm.f2$estim[1]),round(nlm.f2$estim[2]),Rep[1],Rep[2],gainz.round)

if (detail == FALSE)
 { 
   list (output)
 }else
 {
   list (output,output.table)
 }
  
  
}

