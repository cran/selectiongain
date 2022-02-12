# mutistagecor.R author: xuefei mi, first version 17-04-2013, for selectiongain package v2.0.6

# mutistagecor.R author: xuefei mi, second version 11-08-2013, for selectiongain package since v2.0.18

# the package are totally rewritten

# mutistagecor.R author: xuefei mi, 3th version 03-01-2014, for selectiongain package since v2.0.30

#M for type of testers are added

# mutistagecor.R author: xuefei mi, 4th version 11-01-2014, for selectiongain package since v2.0.30

# parental effect was considered

# v45, for longin 2015, genomic selection in wheat paper. 2015-03-01

# v46, add 3 method for parenet-bs1-gca etc. 2015-03-06

# mutistagecor.R author: xuefei mi, 2015-10-08, for selectiongain package since v2.0.35, 40, 47

# defualt values of Vline in covwithoutmas was added

# v48, integrated bug fixing of gsonly scheme in TAG 2015 paper. marked by ####PleaseCheck

# v51, integrated cov calculation for multi traits

# v52 Modifications introduced by JJM on 13.05.2016 are marked as #JJM_index:

# v54, integrated naive GS into the 2traits-1stages, done by xuefei mi

# v55. modified the naive GS for 2 traits one stage. Now it is called "2-traits-1-stage-GS-direct" approach
# according to Schulthess et al 2015 TAG, done by JM
# Implemented the "2-traits-1-stage-GS-reverse" approach for two trais one stage. done by JM

# v60. Dic 05 2016. Introducing last changes for index selection of two traits including GS. Removing the direct and reverse
# creating three diferent indexed and implementing the variations to different stages.

`multistagecor` <-function(maseff=NA,VGCAandE=c(0,0,0,0,0),VSCA=c(0,0,0,0),VLine=c(0,0,0,0,0),
              ecoweight=c(0,rep(1,length(T-1))/length(T-1)), rhop=0,  T, L,M=rep(1,length(T)),
              Rep,index=FALSE,indexTrait=c("Optimum"),covtype=c("LonginII"),detail=FALSE,
              VGCAandE2=c(0,0,0,0,0), VSCA2=c(0,0,0,0), COVgca=c(0,0,0,0,0), COVsca=c(0,0,0,0),
              maseff2=NA, q12=NA, q22=NA) #qx2= proportion of genetic variance explained by molecular markers
{

# we had two optition of VGCA,  either VGCA or VGCAandE
V=VGCAandE
V2=VGCAandE2
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

  if (length(L) != length(M) & covtype=="LonginII-M")
  {
    warning("Location and type of testers do not have same length", call. = FALSE)
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


  VLineType=VLine[1]

  if (VLineType=="VC4" )
  {
    VL=1
    VLl=0.15
    VLy=0.15
    VLly=0.50
    VLe=0.50

  }
  else if (VLineType=="VC5")
  {
     VL=1
    VLl=0.30
    VLy=0.30
    VLly=1
    VLe=1

  } else if (VLineType=="VC6")
  {
     VL=1
    VLl=0.60
    VLy=0.60
    VLly=2
    VLe=2

  }else if(all(is.numeric(VLine))& length (VLine)==5)
  {
    VL=VLine[1]
    VLl=VLine[2]
    VLy=VLine[3]
    VLly=VLine[4]
    VLe=VLine[5]

   }else
   {
     stop("the value of VLine has to be numeric with length 5, or VC4,VC5,VC6, as defined in Longin's paperII")
   }

  if (is.numeric(q12) & is.numeric(q22)) {
    rq1=maseff/sqrt(q12)
    rq2=maseff2/sqrt(q22)
    if (rq1<=1 & rq2<=1){
    }else{
      stop("the accuracy of trait 1 and 2 and the prortion of variance q12 and q22 produce a correlation rq
           larger than 1. Check Dekkers 2007 and recompute your values")
    }
    } else {
  }



# function for calculating cov without markers


covwithoutmas  <-function(V=c(0,0,0,0,0),VSCA=c(0,0,0,0),VLine=c(0,0,0,0,0), ecoweight=c(1,1),
                          rhop, T, L,M=length(T), Rep,index=FALSE,indexTrait=c("Optimum"),covtype=c("LonginII"),V2=c(0,0,0,0,0),
                          VSCA2=c(0,0,0,0), COVgca=c(0,0,0,0,0), COVsca=c(0,0,0,0),maseff=NA, maseff2=NA, q12=NA, q22=NA)
  ####PleaseCheck : I think in this line "VLineVLine" actually corresponds to "VLine"
  # modifyed by mi at 2015-11-05, yes this is a typo
{

# we had two optition of VGCA,  either VGCA or VGCAandE

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

  if (length(L) != length(M) & covtype=="LonginII-M")
  {
    warning("Location and type of testers do not have same length", call. = FALSE)
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


  VLineType=VLine[1]

  if (VLineType=="VC4" )
  {
    VL=1
    VLl=0.15
    VLy=0.15
    VLly=0.50
    VLe=0.50

  }
  else if (VLineType=="VC5")
  {
     VL=1
    VLl=0.30
    VLy=0.30
    VLly=1
    VLe=1

  } else if (VLineType=="VC6")
  {
     VL=1
    VLl=0.60
    VLy=0.60
    VLly=2
    VLe=2

  }else if(all(is.numeric(VLine))& length (VLine)==5)
  {
    VL=VLine[1]
    VLl=VLine[2]
    VLy=VLine[3]
    VLly=VLine[4]
    VLe=VLine[5]

   }else
   {
     stop("the value of VLine has to be numeric with length 5, or VC4,VC5,VC6, as defined in Longin's paperII")
   }



# calculate the covariance

if (covtype=="LonginII")
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
 }else if (covtype=="LonginII-M")
 {
  cov= diag(dim)*Vg
  cov[1,]=Vg
  cov[,1]=Vg
  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1] + Ve/L[i-1]/Rep[i-1]/T[i-1]
  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental")
 {

	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5

  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  cov[1,2]<-covlg*a2+a1*VL
  cov[2,1]=cov[1,2]

	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]

	cov[1,i]=a1*covlg+a2*Vg
	cov[i,1]=cov[1,i]

	# covVLl,Vgl are assumed to be 0, if not 0, add the code here

	}

  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]


  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-gca")
 {

	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5

  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)

  cov[1,1]<-Vg

  cov[1,2]=covlg
  cov[2,1]=cov[1,2]

	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]

	cov[1,i]=Vg
	cov[i,1]=cov[1,i]

	# covVLl,Vgl are assumed to be 0, if not 0, add the code here

	}

  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]


  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-perse")
 {

	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5

  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)

  cov[1,1]=VL

  cov[1,2]=VL
  cov[2,1]=cov[1,2]

	cov[2,2]=VL+VLy+ (VLl+VLly)/L[1]+VLe/L[1]/Rep[1]/T[1]
	 for (i in 3:dim)
  {
	cov[2,i]=covlg
	cov[i,2]=cov[2,i]

	cov[1,i]=covlg
	cov[i,1]=cov[1,i]

	# covVLl,Vgl are assumed to be 0, if not 0, add the code here

	}

  for (i in 3:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]


  }

  for (i in 3:dim)
  {
    for (j in 3:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-BS1")
 {

	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5

  cov= diag(dim)*(a1^2*VL+a2^2*Vg+2*a1*a2*covlg)
  cov[1,2:dim]=a1*covlg+a2*Vg
  cov[2:dim,1]=cov[1,2:dim]

	# if might need to be modifyed if dim>3



  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]


  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="LonginII-Parental-BS1-perse")
 {

	a1=ecoweight[1]
	a2=ecoweight[2]
	covlg=rhop*(Vg*VL)^0.5

  cov= diag(dim)*(VL)
  cov[1,2:dim]=a1*covlg
  cov[2:dim,1]=cov[1,2:dim]

	# if might need to be modifyed if dim>3



  for (i in 2:dim)
  {
    cov[i,i]=Vg+Vgy+(Vgl+Vgly)/L[i-1] + (Vs + Vsy)/T[i-1]/M[i-1] + (Vsl+Vsly)/L[i-1]/T[i-1] /M[i-1]   + Ve/L[i-1]/Rep[i-1]/T[i-1]


  }

  for (i in 2:dim)
  {
    for (j in 2:dim)
    {
      if (i!=j& i<j)
      {
         cov[i,j]=Vg+Vgl/max(L[i-1],L[j-1])+ Vs/max(T[i-1]*M[i-1],T[j-1]*M[j-1]) +Vsl/max(L[i-1],L[j-1])/max(T[i-1]*M[i-1],T[j-1]*M[j-1])
         cov[j,i]=cov[i,j]
      }
     }
  }
 }else if (covtype=="Heffner")
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
 }else if (covtype=="2traits_PS")  #SelIndexProject Scenario 1
 {
   a1=ecoweight[1]
   a2=ecoweight[2]

   covg= diag(3)*Vg
   covp= diag(3)*Vg

   LT=L*T
   LTR=L*T*Rep

   if (length(VGCAandE2)==5& all(is.numeric(VGCAandE2)))
   {
     Vg2=V2[1]
     Vg2l=V2[2]
     Vg2y=V2[3]
     Vg2ly=V2[4]
     Ve2=V2[5]
   }else
   {
     stop("the value of VGCAandE2 has to be numeric with length 5")
   }

   if (length(VSCA2)==4& all(is.numeric(VSCA2)))
   {
     Vs2=VSCA2[1]
     Vsl2=VSCA2[2]
     Vsy2=VSCA2[3]
     Vsly2=VSCA2[4]

   }else
   {
     stop("the value of VSCA2 has to be numeric with length 4")
   }

   if (length(COVgca)==5& all(is.numeric(COVgca)))
   {
     covg12=COVgca[1]
     covg12l=COVgca[2]
     covg12y=COVgca[3]
     covg12ly=COVgca[4]
     cove12=COVgca[5]
   }else
   {
     stop("the value of COVgca has to be numeric with length 5") #JJM_index: VGCAandE2 replaced by COVgca
   }

   if (length(COVsca)==4& all(is.numeric(COVsca)))
   {
     covs12=COVsca[1]
     covs12l=COVsca[2]
     covs12y=COVsca[3]
     covs12ly=COVsca[4]

   }else
   {
     stop("the value of COVsca has to be numeric with length 4") #JJM_index: VGCAandE2 replaced by COVgca
   }

    covg[2,2]=Vg
    covp[2,2]=Vg + Vgy + (Vgl+Vgly)/L[1] + (Vs + Vsy)/T[1] + (Vsl+Vsly)/L[1]/T[1] + Ve/L[1]/Rep[1]/T[1]

    covg[3,3]=Vg2
    covp[3,3]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[1] + (Vs2 + Vsy2)/T[1] + (Vsl2+Vsly2)/L[1]/T[1] + Ve2/L[1]/Rep[1]/T[1]

    covg[1,1]<- a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12   # JJM_index: this computation equals the V(H)=a'Ga see wricke pag 338
    covp[1,1]<- covg[1,1]+0
    # we will not need covp [1,1] temply. so i let the rests to be 0.

    covg[2,3]<-covg12
    covg[3,2]<-covg12

    covp[2,3]<-covg12+covg12y+(covg12ly+covg12l)/L[1]+ (covs12+covs12y)/T[1]+ (covs12ly+covs12l)/L[1]/T[1] + cove12/(L[1]*Rep[1]*T[1]) #JJM_index: added the term +cove12/(L[1]*Rep[1]*T[1])
    covp[3,2]<-covp[2,3]

    if (indexTrait==c("Optimum"))
    {
    vecB<- solve(covp[2:3,2:3]) %*%  covg[2:3,2:3] %*% c(a1,a2)  # JJM_index P-1Ga see wricke pag 339

    b1<-1
    b2<-vecB[2]/vecB[1]
    }
    else if(indexTrait==c("Base"))
    {
    b1<-1
    b2<-a2/a1
    }
    else if(indexTrait==c("Restricted"))
    {
    vecB<- (diag(2)- solve(covp[2:3,2:3]) %*% covg[2:3,3:3] %*%
            solve(t(covg[2:3,3:3]) %*% solve(covp[2:3,2:3]) %*% covg[2:3,3:3]) %*% t(covg[2:3,3:3])) %*%
            solve(covp[2:3,2:3]) %*%  covg[2:3,2:3] %*% c(a1,a2)  # JJM_index see wricke pag 343

    b1<-1
    b2<-vecB[2]/vecB[1]
    }


    cov= diag(dim)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) #JJM_index_2: removed dim-1 and set to dim. This value corresponds to the variance of the target trait or net merit a'*G*a
	  # need to pay attention for higher dimension it is not like that

    cov[2,2]= (b1^2*covp[2,2] + b2^2*covp[3,3] + 2*b1*b2*covp[2,3]) #JJM_index: this is fine and corresponds to V(I)

    cov[2,1]= (a1*b1*covg[2,2] + a2*b2*covg[3,3] + a1*b2*covg[2,3] + a2*b1*covg[2,3]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338

    cov[1,2]=cov[2,1]


    covGeno=covg #JJM_index: we will need this covariance matrix to estimate the gain for each trait separately
    covPheno=covp

    #JJM_index_2 producing correlation matrices for estimating the gain for one trait at a time
    covT1=cov

    covT1[1,1]<-covg[2,2]
    covT1[1,2]<-b1*covg[2,2]+b2*covg[2,3]
    covT1[2,1]<-covT1[1,2]


    corrT1<-cov2cor(covT1)


    covT2=cov

    covT2[1,1]<-covg[3,3]
    covT2[1,2]<-b1*covg[2,3]+b2*covg[3,3]
    covT2[2,1]<-covT2[1,2]

    corrT2<-cov2cor(covT2)






 }else if (covtype=="2traits_GS")   #SelIndexProject Scenario 2 assuming a prediction accuracy for the two traits, one stage selection
 {
   a1=ecoweight[1]
   a2=ecoweight[2]

   covg= diag(3)*Vg
   covp= diag(3)*Vg

   LT=L*T
   LTR=L*T*Rep

   if (length(VGCAandE2)==5& all(is.numeric(VGCAandE2)))
   {
     Vg2=V2[1]
     Vg2l=V2[2]
     Vg2y=V2[3]
     Vg2ly=V2[4]
     Ve2=V2[5]
   }else
   {
     stop("the value of VGCAandE2 has to be numeric with length 5")
   }

   if (length(VSCA2)==4& all(is.numeric(VSCA2)))
   {
     Vs2=VSCA2[1]
     Vsl2=VSCA2[2]
     Vsy2=VSCA2[3]
     Vsly2=VSCA2[4]

   }else
   {
     stop("the value of VSCA2 has to be numeric with length 4")
   }

   if (length(COVgca)==5& all(is.numeric(COVgca)))
   {
     covg12=COVgca[1]
     covg12l=COVgca[2]
     covg12y=COVgca[3]
     covg12ly=COVgca[4]
     cove12=COVgca[5]
   }else
   {
     stop("the value of COVgca has to be numeric with length 5") #JJM_index: VGCAandE2 replaced by COVgca
   }

   if (length(COVsca)==4& all(is.numeric(COVsca)))
   {
     covs12=COVsca[1]
     covs12l=COVsca[2]
     covs12y=COVsca[3]
     covs12ly=COVsca[4]

   }else
   {
     stop("the value of COVsca has to be numeric with length 4") #JJM_index: VGCAandE2 replaced by COVgca
   }

   covg[1,1]<- a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12   # JJM_index: this computation equals the V(H)=a'Ga see wricke pag 338
   covp[1,1]<- covg[1,1]+0
   # we will not need covp [1,1] temply. so i let the rests to be 0.


   #including the stage of genomic prediction
   covg[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
   covp[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
   #Dekkers, J. C. M. (2007). Prediction of response to marker-assisted and genomic selection using selection index theory. Journal of Animal Breeding and Genetics, 124(6), 331-341. doi:10.1111/j.1439-0388.2007.00701.x

   covg[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph
   covp[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph

   rq1=maseff/sqrt(q12)
   rq2=maseff2/sqrt(q22)

   covg[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
   covg[3,2]=covg[2,3]

   covp[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
   covp[3,2]=covp[2,3]


   if (indexTrait==c("Optimum"))
   {

     b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
     b2<- a2/a1   # According to Ceron-rojas 2015 G3 page 2157

   }
   else if(indexTrait==c("Base"))
   {
     b1<-1
     b2<-a2/a1

   }
   else if(indexTrait==c("Restricted"))
   {

     vecA<- (diag(2)- solve(covp[2:3,2:3]) %*% covg[2:3,c(3)] %*%
               solve(t(covg[2:3,c(3)]) %*% solve(covp[2:3,2:3]) %*% covg[2:3,c(3)]) %*% t(covg[2:3,c(3)])) %*%
       c(a1,a2) #new approach. idea of JM


     b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
     b2<- vecA[2]/vecA[1]   # According to Ceron-rojas 2015 G3 page 2157

   }

   #dim= length(L)+1, so in this case we will have dim = 3 because L = GS, PS so length = 2
   cov= diag(dim)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) #JJM_index_2: removed dim-1 and set to dim. It corresponds to the variance of the target trait or net merit a'*G*a
   # need to pay attention for higher dimension it is not like that


   #including the GS stage
   cov[2,2]= (b1^2*covp[2,2] + b2^2*covp[3,3] + 2*b1*b2*covp[2,3])  # as the usuall derivation of sel index
   cov[2,1]= (a1*b1*covg[2,2] + a2*b2*covg[3,3] + a1*b2*covg[2,3] + a2*b1*covg[2,3])   # as the usuall derivation of sel index
   cov[1,2]= cov[2,1]


   covGeno=covg #JJM_index: we will need this covariance matrix to estimate the gain for each trait separately
   covPheno=covp

   #JJM_index_2 producing correlation matrices for estimating the gain for one trait at a time
   covT1=cov

   covT1[1,1]<-Vg
   covT1[1,2]<-b1*covg[2,2]+b2*covg[2,3]
   covT1[2,1]<-covT1[1,2]


   corrT1<-cov2cor(covT1)


   covT2=cov

   covT2[1,1]<-Vg2
   covT2[1,2]<-b1*covg[2,3]+b2*covg[3,3]
   covT2[2,1]<-covT2[1,2]

   corrT2<-cov2cor(covT2)




}else if (covtype=="2traits_GS-PS")   #SelIndexProject Scenario 3 assuming a prediction accuracy for the two traits
 {
   a1=ecoweight[1]
   a2=ecoweight[2]

   covg= diag(5)*Vg # diag is set to 5 to include the accuracies fo both traits
   covp= diag(5)*Vg # diag is set to 5 to include the accuracies fo both traits

   LT=L*T
   LTR=L*T*Rep

   if (length(VGCAandE2)==5& all(is.numeric(VGCAandE2)))
   {
     Vg2=V2[1]
     Vg2l=V2[2]
     Vg2y=V2[3]
     Vg2ly=V2[4]
     Ve2=V2[5]
   }else
   {
     stop("the value of VGCAandE2 has to be numeric with length 5")
   }

   if (length(VSCA2)==4& all(is.numeric(VSCA2)))
   {
     Vs2=VSCA2[1]
     Vsl2=VSCA2[2]
     Vsy2=VSCA2[3]
     Vsly2=VSCA2[4]

   }else
   {
     stop("the value of VSCA2 has to be numeric with length 4")
   }

   if (length(COVgca)==5& all(is.numeric(COVgca)))
   {
     covg12=COVgca[1]
     covg12l=COVgca[2]
     covg12y=COVgca[3]
     covg12ly=COVgca[4]
     cove12=COVgca[5]
   }else
   {
     stop("the value of COVgca has to be numeric with length 5") #JJM_index: VGCAandE2 replaced by COVgca
   }

   if (length(COVsca)==4& all(is.numeric(COVsca)))
   {
     covs12=COVsca[1]
     covs12l=COVsca[2]
     covs12y=COVsca[3]
     covs12ly=COVsca[4]

   }else
   {
     stop("the value of COVsca has to be numeric with length 4") #JJM_index: VGCAandE2 replaced by COVgca
   }

   covg[4,4]=Vg
   covp[4,4]=Vg+Vgy+ (Vgl+Vgly)/L[2]+ (Vs + Vsy)/T[2] + (Vsl+Vsly)/L[2]/T[2] + Ve/L[2]/Rep[2]/T[2]

   covg[5,5]=Vg2
   covp[5,5]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[2] +  (Vs2 + Vsy2)/T[2] + (Vsl2+Vsly2)/L[2]/T[2]  + Ve2/L[2]/Rep[2]/T[2]

   covg[1,1]<- a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12   # JJM_index: this computation equals the V(H)=a'Ga see wricke pag 338
   covp[1,1]<- covg[1,1]+0
   # we will not need covp [1,1] temply. so i let the rests to be 0.

   covg[4,5]<-covg12
   covg[5,4]<-covg12

   covp[4,5]<-covg12+covg12y+(covg12ly+covg12l)/L[2]+ (covs12+covs12y)/T[2]+ (covs12ly+covs12l)/L[2]/T[2] + cove12/(L[2]*Rep[2]*T[2])#JJM_index: added the term +cove12/(L[1]*Rep[1]*T[1])
   covp[5,4]<-covp[4,5]


   #including the stage of genomic prediction
   covg[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
   covp[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
   #Dekkers, J. C. M. (2007). Prediction of response to marker-assisted and genomic selection using selection index theory. Journal of Animal Breeding and Genetics, 124(6), 331-341. doi:10.1111/j.1439-0388.2007.00701.x

   covg[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph
   covp[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph

   rq1=maseff/sqrt(q12)
   rq2=maseff2/sqrt(q22)

   covg[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
   covg[3,2]=covg[2,3]

   covp[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
   covp[3,2]=covp[2,3]


  if (indexTrait==c("Optimum"))
   {
     vecB<- solve(covp[4:5,4:5]) %*%  covg[4:5,4:5] %*%   c(a1,a2) # JJM_index P-1Ga see wricke pag 339. Economic value for markers is set to 0 as they do not have any intrinsic economic value (lande and thopmosn)

     b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
     b2<- a2/a1   # According to Ceron-rojas 2015 G3 page 2157
     #Ceron-Rojas, J. J., Crossa, J., Arief, V. N., Basford, K., Rutkoski, J., Jarqu?n, D., . DeLacy, I. (2015). A Genomic Selection Index Applied to Simulated and Real Data. G3 (Bethesda, Md.), 5(10), 2155-64. doi:10.1534/g3.115.019869
     b3<- 1    #JJM_index_2 According to Prof Utz
     b4<- vecB[2]/vecB[1]# JJM_index_2 according to Prof Utz
   }
   else if(indexTrait==c("Base"))
     {
       b1<-1
       b2<-a2/a1
       b3<-1
       b4<-a2/a1
   }
   else if(indexTrait==c("Restricted"))
   {

     vecA<- (diag(2)- solve(covp[2:3,2:3]) %*% covg[2:3,c(3)] %*%
             solve(t(covg[2:3,c(3)]) %*% solve(covp[2:3,2:3]) %*% covg[2:3,c(3)]) %*% t(covg[2:3,c(3)])) %*%
             c(a1,a2) #new approach. idea of JM
     vecB<- (diag(2)- solve(covp[4:5,4:5]) %*% covg[4:5,c(5)] %*%
             solve(t(covg[4:5,c(5)]) %*% solve(covp[4:5,4:5]) %*% covg[4:5,c(5)]) %*% t(covg[4:5,c(5)])) %*%
             solve(covp[4:5,4:5]) %*%  covg[4:5,4:5] %*% c(a1,a2)  # JJM_index see wricke pag 343

     b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
     b2<- vecA[2]/vecA[1]   # According to Ceron-rojas 2015 G3 page 2157
     b3<- 1
     b4<- vecB[2]/vecB[1]
   }

     #dim= length(L)+1, so in this case we will have dim = 3 because L = GS, PS so length = 2
     cov= diag(dim)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) #JJM_index_2: removed dim-1 and set to dim. It corresponds to the variance of the target trait or net merit a'*G*a
     # need to pay attention for higher dimension it is not like that

     cov[3,3]= (b3^2*covp[4,4] + b4^2*covp[5,5] + 2*b3*b4*covp[4,5]) #JJM_index: this is fine anc orrespond to V(I)

     cov[3,1]= (a1*b3*covg[4,4] + a2*b4*covg[5,5] + a1*b4*covg[4,5] + a2*b3*covg[4,5]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338

     cov[1,3]=cov[3,1]

     #including the GS stage
     cov[2,2]= (b1^2*covp[2,2] + b2^2*covp[3,3] + 2*b1*b2*covp[2,3])  # as the usuall derivation of sel index
     cov[2,1]= (a1*b1*covg[2,2] + a2*b2*covg[3,3] + a1*b2*covg[2,3] + a2*b1*covg[2,3])   # as the usuall derivation of sel index
     cov[1,2]= cov[2,1]


     covx1x1s <- maseff^2*Vg  #according to dekkers 2007 table 1
     covx1x2s <- maseff^2*covg12  #according to dekkers 2007 table 1
     covx2x1s <- maseff2^2*covg12  #according to dekkers 2007 table 1
     covx2x2s <- maseff2^2*Vg2  #according to dekkers 2007 table 1


     cov[3,2]= (b1*b3*covx1x1s + b2*b4*covx2x2s + b1*b4*covx1x2s + b2*b3*covx2x1s) # According to UTZ 1969
     cov[2,3]=cov[3,2]

     covGeno=covg #JJM_index: we will need this covariance matrix to estimate the gain for each trait separately
     covPheno=covp

     #JJM_index_2 producing correlation matrices for estimating the gain for one trait at a time
     covT1=cov

     covT1[1,1]<-covg[4,4]
     covT1[1,2]<-b1*covg[2,2]+b2*covg[2,3]
     covT1[2,1]<-covT1[1,2]
     covT1[1,3]<-b3*covg[4,4]+b4*covg[4,5]
     covT1[3,1]<-covT1[1,3]

     corrT1<-cov2cor(covT1)


     covT2=cov

     covT2[1,1]<-covg[5,5]
     covT2[1,2]<-b1*covg[2,3]+b2*covg[3,3]
     covT2[2,1]<-covT2[1,2]
     covT2[1,3]<-b3*covg[4,5]+b4*covg[5,5]
     covT2[3,1]<-covT2[1,3]

     corrT2<-cov2cor(covT2)

} else if (covtype=="2traits_PS-PS") #SelIndexProject Scenario 4 assuming only PS in two stages for the two traits
 {
   a1=ecoweight[1]
   a2=ecoweight[2]

   covg= diag(5)*Vg
   covp= diag(5)*Vg

   LT=L*T
   LTR=L*T*Rep

   if (length(VGCAandE2)==5& all(is.numeric(VGCAandE2)))
   {
     Vg2=V2[1]
     Vg2l=V2[2]
     Vg2y=V2[3]
     Vg2ly=V2[4]
     Ve2=V2[5]
   }else
   {
     stop("the value of VGCAandE2 has to be numeric with length 5")
   }

   if (length(VSCA2)==4& all(is.numeric(VSCA2)))
   {
     Vs2=VSCA2[1]
     Vsl2=VSCA2[2]
     Vsy2=VSCA2[3]
     Vsly2=VSCA2[4]

   }else
   {
     stop("the value of VSCA2 has to be numeric with length 4")
   }

   if (length(COVgca)==5& all(is.numeric(COVgca)))
   {
     covg12=COVgca[1]
     covg12l=COVgca[2]
     covg12y=COVgca[3]
     covg12ly=COVgca[4]
     cove12=COVgca[5]
   }else
   {
     stop("the value of COVgca has to be numeric with length 5") #JJM_index: VGCAandE2 replaced by COVgca
   }

   if (length(COVsca)==4& all(is.numeric(COVsca)))
   {
     covs12=COVsca[1]
     covs12l=COVsca[2]
     covs12y=COVsca[3]
     covs12ly=COVsca[4]

   }else
   {
     stop("the value of COVsca has to be numeric with length 4") #JJM_index: VGCAandE2 replaced by COVgca
   }

   covg[2,2]=Vg
   covp[2,2]=Vg + Vgy + (Vgl+Vgly)/L[1] + (Vs + Vsy)/T[1] + (Vsl+Vsly)/L[1]/T[1] + Ve/L[1]/Rep[1]/T[1]

   covg[3,3]=Vg2
   covp[3,3]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[1] + (Vs2 + Vsy2)/T[1] + (Vsl2+Vsly2)/L[1]/T[1] + Ve2/L[1]/Rep[1]/T[1]


   covg[4,4]=Vg
   covp[4,4]=Vg+Vgy+ (Vgl+Vgly)/L[2]+ (Vs + Vsy)/T[2] + (Vsl+Vsly)/L[2]/T[2] + Ve/L[2]/Rep[2]/T[2]

   covg[5,5]=Vg2
   covp[5,5]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[2] +  (Vs2 + Vsy2)/T[2] + (Vsl2+Vsly2)/L[2]/T[2]  + Ve2/L[2]/Rep[2]/T[2]


   covg[1,1]<- a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12   # JJM_index: this computation equals the V(H)=a'Ga see wricke pag 338
   covp[1,1]<- covg[1,1]+0
   # we will not need covp [1,1] temply. so i let the rests to be 0.

   covg[2,3]<-covg12
   covg[3,2]<-covg12

   covp[2,3]<-covg12+covg12y+(covg12ly+covg12l)/L[1]+ (covs12+covs12y)/T[1]+ (covs12ly+covs12l)/L[1]/T[1] + cove12/(L[1]*Rep[1]*T[1]) #JJM_index: added the term +cove12/(L[1]*Rep[1]*T[1])
   covp[3,2]<-covp[2,3]

   covg[4,5]<-covg12
   covg[5,4]<-covg12

   covp[4,5]<-covg12+covg12y+(covg12ly+covg12l)/L[2]+ (covs12+covs12y)/T[2]+ (covs12ly+covs12l)/L[2]/T[2] + cove12/(L[2]*Rep[2]*T[2])#JJM_index: added the term +cove12/(L[1]*Rep[1]*T[1])
   covp[5,4]<-covp[4,5]


   if (indexTrait==c("Optimum"))
   {
	 vecB<- solve(covp[2:5,2:5]) %*% covg[2:5,2:5] %*%  c(a1,a2,a1,a2) # JJM_index P-1Ga see wricke pag 339

     b1<-1
     b2<-vecB[2]/vecB[1]
     b3<-1 #JJM_index_2 According to Prof Utz
	   #b3<-vecB[3]/vecB[1] #Old code
     b4<-vecB[4]/vecB[3]# JJM_index_2 according to Prof Utz
	   #b4<-vecB[4]/vecB[1]# Old Code
   }
   else if(indexTrait==c("Base"))
   {
     b1<-1
     b2<-a2/a1
     b3<-1
     b4<-a2/a1
   }
      else if(indexTrait==c("Restricted"))
   {
     vecB<- (diag(4)- solve(covp[2:5,2:5]) %*% covg[2:5,c(3,5)] %*%
               solve(t(covg[2:5,c(3,5)]) %*% solve(covp[2:5,2:5]) %*% covg[2:5,c(3,5)]) %*% t(covg[2:5,c(3,5)])) %*%
       solve(covp[2:5,2:5]) %*%  covg[2:5,2:5] %*% c(a1,a2,a1,a2)  # JJM_index see wricke pag 343

     b1<-1
     b2<-vecB[2]/vecB[1]
     b3<-1 #JJM_index_2 According to Prof Utz
     #b3<-vecB[3]/vecB[1] #Old code
     b4<-vecB[4]/vecB[3]# JJM_index_2 according to Prof Utz
     #b4<-vecB[4]/vecB[1]# Old Code
   }


	   cov= diag(dim)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) #JJM_index_2: removed dim-1 and set to dim
	   #cov= diag(dim-1)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) # what about this code for the position [1,1], it correspond to the variance of the target trait or net merit a'*G*a
     #cov= diag(dim-1)* (b1^2*Vg + b2^2*Vg2 + 2*b1*b2*covg12) #JJM_index: That was the old code


     # need to pay attention for higher dimension it is not like that

     cov[2,2]= (b1^2*covp[2,2] + b2^2*covp[3,3] + 2*b1*b2*covp[2,3]) #JJM_index: this is fine anc orrespond to V(I)

     cov[2,1]= (a1*b1*covg[2,2] + a2*b2*covg[3,3] + a1*b2*covg[2,3] + a2*b1*covg[2,3]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338
     #cov[2,1]= (b1^2*Vg + b2^2*Vg2 + 2*b1*b2*covg12) #JJM_index: that was the old line of code

     cov[1,2]=cov[2,1]

     cov[3,3]= (b3^2*covp[4,4] + b4^2*covp[5,5] + 2*b3*b4*covp[4,5]) #JJM_index: this is fine anc orrespond to V(I)

     cov[3,1]= (a1*b3*covg[4,4] + a2*b4*covg[5,5] + a1*b4*covg[4,5] + a2*b3*covg[4,5]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338
     #cov[2,1]= (b1^2*Vg + b2^2*Vg2 + 2*b1*b2*covg12) #JJM_index: that was the old line of code

     cov[1,3]=cov[3,1]


	 covx1x1s <- Vg + (Vgl/max(L[1],L[2])) + (Vs)/max(T[1],T[2]) + (Vsl)/max(L[1],L[2])/max(T[1],T[2]) # See Longin 2007 paper II equation 2
	 #we divide by the max becuase the formula of Utz is: (Vgl*CommonLoc)/loc1*Loc2 and we assume that Loc1 are the common locs
	 covx1x2s <- covg12+(covg12l/max(L[1],L[2])) + (covs12/max(T[1],T[2])) + ((covs12l)/max(L[1],L[2])/max(T[1],T[2])) #new line according to suggestion of Prof. UTZ
	 covx2x1s <- covg12+(covg12l/max(L[1],L[2])) + (covs12/max(T[1],T[2])) + ((covs12l)/max(L[1],L[2])/max(T[1],T[2]))
	 covx2x2s <- Vg2 + (Vg2l/max(L[1],L[2])) + (Vs2)/max(T[1],T[2]) + (Vsl2)/max(L[1],L[2])/max(T[1],T[2])  #new line


	 cov[3,2]= (b1*b3*covx1x1s + b2*b4*covx2x2s + b1*b4*covx1x2s + b2*b3*covx2x1s)
   cov[2,3]=cov[3,2]

   covGeno=covg #JJM_index: we will need this covariance matrix to estimate the gain for each trait separately
   covPheno=covp

     #JJM_index_2 producing correlation matrices for estimating the gain for one trait at a time
     covT1=cov

	   covT1[1,1]<-covg[2,2]
     covT1[1,2]<-b1*covg[2,2]+b2*covg[2,3]
     covT1[2,1]<-covT1[1,2]
     covT1[1,3]<-b3*covg[4,4]+b4*covg[4,5]
     covT1[3,1]<-covT1[1,3]

     corrT1<-cov2cor(covT1)


     covT2=cov

	   covT2[1,1]<-covg[3,3]
     covT2[1,2]<-b1*covg[2,3]+b2*covg[3,3]
     covT2[2,1]<-covT2[1,2]
     covT2[1,3]<-b3*covg[4,5]+b4*covg[5,5]
     covT2[3,1]<-covT2[1,3]

     corrT2<-cov2cor(covT2)




}else if (covtype=="2traits_GS-PS-PS")   #SelIndexProject Scenario 5 assuming GA and PS in three stages for the two traits
{
  a1=ecoweight[1]
  a2=ecoweight[2]

  covg= diag(7)*Vg # diag is set to 7 to include the accuracies fo both traits and the two phenotypic stages
  covp= diag(7)*Vg # diag is set to 7 to include the accuracies fo both traits and the two phenotypic stages

  LT=L*T
  LTR=L*T*Rep

  if (length(VGCAandE2)==5& all(is.numeric(VGCAandE2)))
  {
    Vg2=V2[1]
    Vg2l=V2[2]
    Vg2y=V2[3]
    Vg2ly=V2[4]
    Ve2=V2[5]
  }else
  {
    stop("the value of VGCAandE2 has to be numeric with length 5")
  }

  if (length(VSCA2)==4& all(is.numeric(VSCA2)))
  {
    Vs2=VSCA2[1]
    Vsl2=VSCA2[2]
    Vsy2=VSCA2[3]
    Vsly2=VSCA2[4]

  }else
  {
    stop("the value of VSCA2 has to be numeric with length 4")
  }

  if (length(COVgca)==5& all(is.numeric(COVgca)))
  {
    covg12=COVgca[1]
    covg12l=COVgca[2]
    covg12y=COVgca[3]
    covg12ly=COVgca[4]
    cove12=COVgca[5]
  }else
  {
    stop("the value of COVgca has to be numeric with length 5") #JJM_index: VGCAandE2 replaced by COVgca
  }

  if (length(COVsca)==4& all(is.numeric(COVsca)))
  {
    covs12=COVsca[1]
    covs12l=COVsca[2]
    covs12y=COVsca[3]
    covs12ly=COVsca[4]

  }else
  {
    stop("the value of COVsca has to be numeric with length 4") #JJM_index: VGCAandE2 replaced by COVgca
  }

  covg[4,4]=Vg
  covp[4,4]=Vg + Vgy + (Vgl+Vgly)/L[2] + (Vs + Vsy)/T[2] + (Vsl+Vsly)/L[2]/T[2] + Ve/L[2]/Rep[2]/T[2]

  covg[5,5]=Vg2
  covp[5,5]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[2] + (Vs2 + Vsy2)/T[2] + (Vsl2+Vsly2)/L[2]/T[2] + Ve2/L[2]/Rep[2]/T[2]


  covg[6,6]=Vg
  covp[6,6]=Vg+Vgy+ (Vgl+Vgly)/L[3]+ (Vs + Vsy)/T[3] + (Vsl+Vsly)/L[3]/T[3] + Ve/L[3]/Rep[3]/T[3]

  covg[7,7]=Vg2
  covp[7,7]=Vg2 + Vg2y + (Vg2l+Vg2ly)/L[3] +  (Vs2 + Vsy2)/T[3] + (Vsl2+Vsly2)/L[3]/T[3]  + Ve2/L[3]/Rep[3]/T[3]


  covg[1,1]<- a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12   # JJM_index: this computation equals the V(H)=a'Ga see wricke pag 338
  covp[1,1]<- covg[1,1]+0
  # we will not need covp [1,1] temply. so i let the rests to be 0.

  covg[4,5]<-covg12
  covg[5,4]<-covg12

  covp[4,5]<-covg12+covg12y+(covg12ly+covg12l)/L[2]+ (covs12+covs12y)/T[2]+ (covs12ly+covs12l)/L[2]/T[2] + cove12/(L[2]*Rep[2]*T[2])#JJM_index: added the term +cove12/(L[1]*Rep[1]*T[1])
  covp[5,4]<-covp[4,5]

  covg[6,7]<-covg12
  covg[7,6]<-covg12

  covp[6,7]<-covg12+covg12y+(covg12ly+covg12l)/L[3]+ (covs12+covs12y)/T[3]+ (covs12ly+covs12l)/L[3]/T[3] + cove12/(L[3]*Rep[3]*T[3])
  covp[7,6]<-covp[6,7]




  #including the stage of genomic prediction
  covg[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
  covp[2,2]=Vg*maseff^2  #according to dekkers 2007 page 334 second paragraph
  #Dekkers, J. C. M. (2007). Prediction of response to marker-assisted and genomic selection using selection index theory. Journal of Animal Breeding and Genetics, 124(6), 331-341. doi:10.1111/j.1439-0388.2007.00701.x

  covg[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph
  covp[3,3]=Vg2*maseff2^2  #according to dekkers 2007 page 334 second paragraph

  rq1=maseff/sqrt(q12)
  rq2=maseff2/sqrt(q22)

  covg[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
  covg[3,2]=covg[2,3]

  covp[2,3]= rq1*rq2*maseff*maseff2*covg12  # derived from dekkers 2007 page 335, see Joses derivation in notebook
  covp[3,2]=covp[2,3]


  if (indexTrait==c("Optimum"))
  {
    vecB<- solve(covp[4:7,4:7]) %*%  covg[4:7,4:7] %*%   c(a1,a2,a1,a2) # JJM_index P-1Ga see wricke pag 339. Economic value for markers is set to 0 as they do not have any intrinsic economic value (lande and thopmosn)

    b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
    b2<- a2/a1   # According to Ceron-rojas 2015 G3 page 2157
    #Ceron-Rojas, J. J., Crossa, J., Arief, V. N., Basford, K., Rutkoski, J., Jarqu?n, D., . DeLacy, I. (2015). A Genomic Selection Index Applied to Simulated and Real Data. G3 (Bethesda, Md.), 5(10), 2155-64. doi:10.1534/g3.115.019869
    b3<- 1    #JJM_index_2 According to Prof Utz
    b4<- vecB[2]/vecB[1]# JJM_index_2 according to Prof Utz
    b5<- 1    #JJM_index_2 According to Prof Utz
    b6<- vecB[4]/vecB[3]# JJM_index_2 according to Prof Utz
  }
  else if(indexTrait==c("Base"))
  {
    b1<-1
    b2<-a2/a1
    b3<-1
    b4<-a2/a1
    b5<-1
    b6<-a2/a1
  }
  else if(indexTrait==c("Restricted"))
  {

    vecA<- (diag(2)- solve(covp[2:3,2:3]) %*% covg[2:3,c(3)] %*%
           solve(t(covg[2:3,c(3)]) %*% solve(covp[2:3,2:3]) %*% covg[2:3,c(3)]) %*% t(covg[2:3,c(3)])) %*%
           c(a1,a2) #new approach. idea of JM
    vecB<- (diag(4)- solve(covp[4:7,4:7]) %*% covg[4:7,c(5,7)] %*%
            solve(t(covg[4:7,c(5,7)]) %*% solve(covp[4:7,4:7]) %*% covg[4:7,c(5,7)]) %*% t(covg[4:7,c(5,7)])) %*%
            solve(covp[4:7,4:7]) %*%  covg[4:7,4:7] %*%   c(a1,a2,a1,a2)  # JJM_index see wricke pag 343

    b1<- 1   # According to Ceron-rojas 2015 G3 page 2157
    b2<- vecA[2]/vecA[1]   # According to Ceron-rojas 2015 G3 page 2157
    #Ceron-Rojas, J. J., Crossa, J., Arief, V. N., Basford, K., Rutkoski, J., Jarqu?n, D., . DeLacy, I. (2015). A Genomic Selection Index Applied to Simulated and Real Data. G3 (Bethesda, Md.), 5(10), 2155-64. doi:10.1534/g3.115.019869
    b3<- 1    #JJM_index_2 According to Prof Utz
    b4<- vecB[2]/vecB[1]# JJM_index_2 according to Prof Utz
    b5<- 1    #JJM_index_2 According to Prof Utz
    b6<- vecB[4]/vecB[3]# JJM_index_2 according to Prof Utz
  }


    #dim= length(L)+1, so in this case we will have dim = 3 because L = GS, PS so length = 2
    cov= diag(dim)* (a1^2*Vg + a2^2*Vg2 + 2*a1*a2*covg12) #JJM_index_2: removed dim-1 and set to dim. It corresponds to the variance of the target trait or net merit a'*G*a
    # need to pay attention for higher dimension it is not like that

    cov[3,3]= (b3^2*covp[4,4] + b4^2*covp[5,5] + 2*b3*b4*covp[4,5]) #JJM_index: this is fine anc orrespond to V(I)

    cov[3,1]= (a1*b3*covg[4,4] + a2*b4*covg[5,5] + a1*b4*covg[4,5] + a2*b3*covg[4,5]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338

    cov[1,3]=cov[3,1]

    cov[4,4]= (b5^2*covp[6,6] + b6^2*covp[7,7] + 2*b5*b6*covp[6,7]) #JJM_index: this is fine anc orrespond to V(I)

    cov[4,1]= (a1*b5*covg[6,6] + a2*b6*covg[7,7] + a1*b6*covg[6,7] + a2*b5*covg[6,7]) # here we need the covariance cov(H,I) which in matrix notation is b'*G*a see Wricke pag 338
    #cov[2,1]= (b1^2*Vg + b2^2*Vg2 + 2*b1*b2*covg12) #JJM_index: that was the old line of code

    cov[1,4]=cov[4,1]

    covx1x1s <- Vg + (Vgl/max(L[2],L[3])) +(Vs)/max(T[2],T[3]) + (Vsl)/max(L[2],L[3])/max(T[2],T[3]) #new line
    #we divide by the max because the formula of Utz is: (Vgl*CommonLoc)/loc1*Loc2 and we assume that Loc1 are the common locs
    covx1x2s <- covg12+(covg12l/max(L[2],L[3])) + (covs12/max(T[2],T[3])) + ((covs12l)/max(L[2],L[3])/max(T[2],T[3])) #new line according to suggestion of Prof. UTZ
    covx2x1s <- covg12+(covg12l/max(L[2],L[3])) + (covs12/max(T[2],T[3])) + ((covs12l)/max(L[2],L[3])/max(T[2],T[3]))
    covx2x2s <- Vg2 + (Vg2l/max(L[2],L[3])) + (Vs2)/max(T[2],T[3]) + (Vsl2)/max(L[2],L[3])/max(T[2],T[3])#new line


    cov[4,3]= (b3*b5*covx1x1s + b4*b6*covx2x2s + b3*b6*covx1x2s + b4*b5*covx2x1s)
    cov[3,4]=cov[4,3]


    #including the GS stage
    cov[2,2]= (b1^2*covp[2,2] + b2^2*covp[3,3] + 2*b1*b2*covp[2,3])  # as the usuall derivation of sel index
    cov[2,1]= (a1*b1*covg[2,2] + a2*b2*covg[3,3] + a1*b2*covg[2,3] + a2*b1*covg[2,3])   # as the usuall derivation of sel index
    cov[1,2]= cov[2,1]


    covx1x1s <- maseff^2*Vg  #according to dekkers 2007 table 1
    covx1x2s <- maseff^2*covg12  #according to dekkers 2007 table 1
    covx2x1s <- maseff2^2*covg12  #according to dekkers 2007 table 1
    covx2x2s <- maseff2^2*Vg2  #according to dekkers 2007 table 1


    cov[3,2]= (b1*b3*covx1x1s + b2*b4*covx2x2s + b1*b4*covx1x2s + b2*b3*covx2x1s) # According to UTZ 1969
    cov[2,3]=cov[3,2]

    cov[4,2]= (b1*b5*covx1x1s + b2*b6*covx2x2s + b1*b6*covx1x2s + b2*b5*covx2x1s) # According to UTZ 1969
    cov[2,4]=cov[4,2]

    covGeno=covg #JJM_index: we will need this covariance matrix to estimate the gain for each trait separately
    covPheno=covp


    #JJM_index_2 producing correlation matrices for estimating the gain for one trait at a time
    covT1=cov

    covT1[1,1]<-covg[4,4]
    covT1[1,2]<-b1*covg[2,2]+b2*covg[2,3]
    covT1[2,1]<-covT1[1,2]
    covT1[1,3]<-b3*covg[4,4]+b4*covg[4,5]
    covT1[3,1]<-covT1[1,3]
    covT1[1,4]<-b5*covg[6,6]+b6*covg[6,7]
    covT1[4,1]<-covT1[1,4]

    corrT1<-cov2cor(covT1)


    covT2=cov

    covT2[1,1]<-covg[5,5]
    covT2[1,2]<-b1*covg[2,3]+b2*covg[3,3]
    covT2[2,1]<-covT2[1,2]
    covT2[1,3]<-b3*covg[4,5]+b4*covg[5,5]
    covT2[3,1]<-covT2[1,3]
    covT2[1,4]<-b5*covg[6,7]+b6*covg[7,7]
    covT2[4,1]<-covT2[1,4]

    corrT2<-cov2cor(covT2)




}else
 {
  stop("covtype is not specified, see Longin's paperII")
 }

# calculate the optimal selection index = G^-1 /P

  if (covtype=="2traits_PS")
  {
    list(cov,c(b1,b2),corrT1, corrT2,covGeno, covPheno)  #JJM_index: please fell free to delete b3 and b4 as those were used just for comparison
    # covGeno is needed to compute the gain for each trait
  }
  else if (covtype=="2traits_GS")
  {
    list(cov,c(b1,b2),corrT1, corrT2,covGeno, covPheno)  #JJM_index: please fell free to delete b3 and b4 as those were used just for comparison
    # covGeno is needed to compute the gain for each trait
  }
  else if (covtype=="2traits_PS-PS")
  {
    list(cov,c(b1,b2,b3,b4),corrT1,corrT2,covGeno, covPheno)  #JJM_index_2: corrT1 and corrT2 are needed to compute the gain for each trait
  }
  else if (covtype=="2traits_GS-PS")
  {
    list(cov,c(b1,b2,b3,b4),corrT1, corrT2,covGeno, covPheno)  #JJM_index_2: corrT1 and corrT2 are needed to compute the gain for each trait
  }
  else if (covtype=="2traits_GS-PS-PS")
  {
    list(cov,c(b1,b2,b3,b4,b5,b6),corrT1, corrT2,covGeno, covPheno)  #JJM_index_2: corrT1 and corrT2 are needed to compute the gain for each trait
  }
  else
  {
    cov
  }


}



# end of covwithoutmas

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

          if (covtype=="LonginII-Parental-BS1")
          {
            cov=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep, covtype=covtype)

            output= list(cov2cor(cov),tempb,cov)
          }else if(covtype=="2traits_PS")
          {
            output=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep,
                                 covtype=covtype,indexTrait=indexTrait, V2=VGCAandE2,VSCA2=VSCA2, COVgca=COVgca,
                                 COVsca=COVsca, maseff=maseff, maseff2=maseff2, q12=q12, q22=q22)
            cov=output[[1]]
            tempb=output[[2]]
            corrT1=output[[3]]
            corrT2=output[[4]]
            covGeno=output[[5]]
            covPheno=output[[6]]
            output= list(cov2cor(cov),tempb,cov,corrT1,corrT2,covGeno,covPheno)
          }else if(covtype=="2traits_PS-PS")
          {
            output=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep,
                                 covtype=covtype,indexTrait=indexTrait, V2=VGCAandE2,VSCA2=VSCA2, COVgca=COVgca,
                                 COVsca=COVsca, maseff=maseff, maseff2=maseff2, q12=q12, q22=q22)
            cov=output[[1]]
            tempb=output[[2]]
            corrT1=output[[3]]
		       	corrT2=output[[4]]
		       	covGeno=output[[5]]
		       	covPheno=output[[6]]
		       	output= list(cov2cor(cov),tempb,cov,corrT1,corrT2,covGeno,covPheno)
          }
          else
          {
              tempb="empty"
	            cov=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep, covtype=covtype,index=index,indexTrait=indexTrait)
				           a1=ecoweight[1]
	                 a2=ecoweight[2]
	                 covlg=rhop*(Vg*VL)^0.5


              if (index==TRUE)
             {
                P=cov[-1,-1]
                G=matrix(rep(1,(dim-1)^2), nrow = dim-1, ncol=dim-1, byrow=TRUE) * Vg
                if (covtype=="LonginII-Parental" | covtype== "LonginII-Parental-gca" | covtype== "LonginII-Parental-perse")
									 {
									  G[1,1]=VL
										# genotypic covariance of stage 1

									  covlg=rhop*(Vg*VL)^0.5
									   for (i in 2:c(dim-1))
                     {
	                     G[1,i]=covlg
	                     G[i,1]=G[1,i]

	                     # covVLl,Vgl are assumed to be 0, if not 0, add the code here

	                    }

									 }





                for (i in 2:c(dim-1))
                {
                   tempb=matrix(c(a1,rep((1-a1)/(dim-2),i-1)),nrow=1,ncol=i,byrow=TRUE)

                   tempp=P[1:i,1:i]
                   tempg=G[1:i,1:i]
                   tempb= t((solve(tempp)%*% tempg )  %*% t( tempb))

                   tempb=tempb / sum(tempb)
                   	if (covtype=="LonginII-Parental-BS1")
									 {
									  	if (tempb[1]!=0)
									  	{
									  	tempb=tempb / tempb[1]
									  	}else
									  	{
									  	tempb=tempb / sum(tempb)
									  	}
									  }

                   P[i,i]=  tempb %*% tempp %*% t(tempb)
                   P[1:(i-1),i]= tempb %*% tempp[,1:(i-1)]
                   P[i,1:(i-1)]=t(P[1:(i-1),i])
                }

# calculate the covariance and give output
                   cov[2:dim,2:dim]=P

									if (covtype=="LonginII-Parental")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]

									# formulation have to be improved when dim >3

									 cov[1,3:dim]=a1*b1*VL+a2*b2*Vg+(a1*b2+a2*b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else if (covtype=="LonginII-Parental-gca")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]

									# formulation have to be improved when dim >3

									 cov[1,3:dim]=b2*Vg+(b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else if (covtype=="LonginII-Parental-perse")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                   # 	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]

									# formulation have to be improved when dim >3

									 cov[1,3:dim]=b1*VL+(b2)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }else 	if (covtype=="LonginII-Parental-BS1")
									 {
									  	a1=ecoweight[1]
                    	a2=ecoweight[2]
                    #	tempb=tempb / tempb[1]
                      b1=tempb[1]
											b2=tempb[2]

									# formulation have to be improved when dim >3

										cov[1,2]=a2*Vg+a1*covlg
									  cov[2,1]=cov[1,2]

									 cov[1,3:dim]=(a2*b1+a2*b2)*Vg+(a1*b2+a1*b1)*covlg
									  cov[3:dim,1]=cov[1,3:dim]
									 }


              }else if (index!=TRUE)
							{

							}

                   output= list(cov2cor(cov),tempb,cov)
        }
	   } else if(!is.na(maseff))
       {
             if (length(maseff)!=1)
             {
                 stop("if maseff is not NA the then the length of it has to be 1.")
             }

	     if(covtype=="2traits_GS")
	     {
	       output=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep,
	                            covtype=covtype,indexTrait=indexTrait, V2=VGCAandE2,VSCA2=VSCA2, COVgca=COVgca,
	                            COVsca=COVsca, maseff=maseff, maseff2=maseff2, q12=q12, q22=q22)

	       cov=output[[1]]
	       tempb=output[[2]]
	       corrT1=output[[3]]
	       corrT2=output[[4]]
	       covGeno=output[[5]]
	       covPheno=output[[6]]

	       output= list(cov2cor(cov),tempb,cov,corrT1,corrT2,covGeno,covPheno)

	       if (index==TRUE)
	       {
	         warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
	       }
	     }

	     else if(covtype=="2traits_GS-PS")
	       {
	       output=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep,
	                            covtype=covtype,indexTrait=indexTrait, V2=VGCAandE2,VSCA2=VSCA2, COVgca=COVgca,
	                            COVsca=COVsca, maseff=maseff, maseff2=maseff2, q12=q12, q22=q22)

	       cov=output[[1]]
	       tempb=output[[2]]
	       corrT1=output[[3]]
	       corrT2=output[[4]]
	       covGeno=output[[5]]
	       covPheno=output[[6]]

	       output= list(cov2cor(cov),tempb,cov,corrT1,corrT2,covGeno,covPheno)

	         if (index==TRUE)
	         {
	           warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
	         }
	       }

	     else if(covtype=="2traits_GS-PS-PS")
	        {
	       output=covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine, ecoweight=ecoweight,rhop=rhop,T=T,L=L,M=M,Rep=Rep,
	                            covtype=covtype,indexTrait=indexTrait, V2=VGCAandE2,VSCA2=VSCA2, COVgca=COVgca,
	                            COVsca=COVsca, maseff=maseff, maseff2=maseff2, q12=q12, q22=q22)

	       cov=output[[1]]
	       tempb=output[[2]]
	       corrT1=output[[3]]
	       corrT2=output[[4]]
	       covGeno=output[[5]]
	       covPheno=output[[6]]

	       output= list(cov2cor(cov),tempb,cov,corrT1,corrT2,covGeno,covPheno)

	       if (index==TRUE)
	        {
	         warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
	        }
	       }

	     else {

	             if (covtype=="LonginII-Parental")
             {
                cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
								covtype=covtype)
               	dim=dim-1
                corp=cov
								covlg=rhop*(Vg*VL)^0.5
	                cormas= diag(dim+1)
		            cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
		            cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
		            cormas[3:c(dim+1),1]=corp[2:c(dim),1]
	            	cormas[1,2]=maseff^2*corp[1,3]
		            cormas[2,1]=maseff^2*corp[3,1]
	            	cormas[2,2]=maseff^2*Vg
	              cormas[1,1]=corp[1,1]

                    cormas[2,3]=maseff^2*corp[2,3]
	 	           cormas[3,2]=cormas[2,3]

	          	for (i in 4:c(dim+1))
	      	   {
		           cormas[2,i]=maseff^2*Vg
	 	           cormas[i,2]=cormas[2,i]
	           }




             tempb="empty"
             if (index==TRUE)
               {
                warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }
             }else if(covtype=="LonginII-Parental-gca")
             {

                 cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                    covtype=covtype)
                 dim=dim-1
                 corp=cov
                 cormas= diag(dim+1)
                 cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
                 cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
                 cormas[3:c(dim+1),1]=corp[2:c(dim),1]
                 cormas[1,2]=maseff^2*Vg
                 cormas[2,1]=maseff^2*Vg
                 cormas[2,2]=maseff^2*Vg
                 cormas[1,1]=corp[1,1]

                 cormas[2,3]=maseff^2*corp[1,2]
                 cormas[3,2]=cormas[2,3]

                 for (i in 4:c(dim+1))
                 {
                   cormas[2,i]=maseff^2*Vg
                   cormas[i,2]=cormas[2,i]
                 }



               tempb="empty"
               if (index==TRUE)
               {
                 warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }

             }else if(covtype=="LonginII-Parental-perse")
             {

               cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                  covtype=covtype)
               dim=dim-1
               corp=cov
               cormas= diag(dim+1)
               cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
               cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
               cormas[3:c(dim+1),1]=corp[2:c(dim),1]
               cormas[1,2]=maseff^2*corp[1,3]
               cormas[2,1]=maseff^2*corp[1,3]
               cormas[2,2]=maseff^2*Vg
               cormas[1,1]=corp[1,1]

               cormas[2,3]=maseff^2*corp[1,3]
               cormas[3,2]=cormas[2,3]

               for (i in 4:c(dim+1))
               {
                 cormas[2,i]=maseff^2*Vg
                 cormas[i,2]=cormas[2,i]
               }




             tempb="empty"
              if (index==TRUE)
              {
                warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
              }


             }else if(covtype=="LonginII-Parental-BS1")
             {

               cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,VLine=VLine,ecoweight=ecoweight,rhop=rhop,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
                                  covtype=covtype)
               dim=dim-1
               corp=cov
               cormas= diag(dim+1)
               cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
               cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
               cormas[3:c(dim+1),1]=corp[2:c(dim),1]
               cormas[1,2]=maseff^2*corp[1,3]
               cormas[2,1]=maseff^2*corp[1,3]
               cormas[2,2]=maseff^2*Vg
               cormas[1,1]=corp[1,1]

               cormas[2,3]=maseff^2
               cormas[3,2]=cormas[2,3]

               for (i in 4:c(dim+1))
               {
                 cormas[2,i]=maseff^2*Vg
                 cormas[i,2]=cormas[2,i]
               }




               tempb="empty"
               if (index==TRUE)
               {
                 warning("Heffner's equation is kind of index, no optimal index calculate will be executed",                 call. = FALSE)
               }


             }else if (index!=TRUE & covtype!="LonginII-Parental"& covtype!="LonginII-Parental-gca"& covtype!="LonginII-Parental-perse"& covtype!="LonginII-Parental-BS1")
             {
                cov= covwithoutmas(V=VGCAandE,VSCA=VSCA,T=T[-1], L=L[-1],Rep=Rep[-1],M=M[-1],
								covtype=covtype)
               	dim=dim-1
                corp=cov

	              cormas= diag(dim+1)
		            cormas[3:c(dim+1),3:c(dim+1)]=corp[2:c(dim),2:c(dim)]
		            cormas[1,3:c(dim+1)]=corp[1,2:c(dim)]
		            cormas[3:c(dim+1),1]=corp[2:c(dim),1]
	            	cormas[1,2]=maseff^2*corp[2,1]  ####PleaseCheck  : if the Length of L is 2 and we are using mas, corp will be of dimentions 2x2. Then corp[3,1] does not exist. The construction of the
					                                                   # cov matrix (cov=corp) indicates that for covtype=LonginII corp[3,1]=corp[2,1]=corp[1,1]=Vg So we can replace corp[3,1] by
																	   # corp[2,1] or corp[1,1] or Vg . I replace it by corp[2,1] and the packege worked. Is this change correct?
	            	                     # modifyed by mi 2015.11.05, yes you can use corp[2,1]

		            cormas[2,1]=cormas[1,2]
	            	cormas[2,2]=maseff^2*Vg
	              cormas[1,1]=corp[1,1]

	          	for (i in 3:c(dim+1))
	      	   {
		           cormas[2,i]=maseff^2*Vg
	 	           cormas[i,2]=cormas[2,i]
	           }



             }

	           tempb="empty"
             output= list(cov2cor(cormas),tempb,cormas)
	   }
   }
if (detail==TRUE)
{
output
}else
{
output[[1]]
}

}












