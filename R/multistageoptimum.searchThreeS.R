
# package selectiongain

# modified at 28-06-2013, for 1MAS+2PS
# modified at 13.08.2015 by Jose Marulanda for Preselection on Nurseries for an uncorrelated trait.
# modified at 13.08.2015 by Jose Marulanda for a better use of the testcross seed production information
# All modified lines or added lines end with #JM

# modified at 2015-11-09,
# 1, for the warning of "budget too small",  changed from stop to warning, make more freedom
# 2, try to use matrix multiplication to relace vactor multimplication such that the RRO3.2.2 and MKL can take advantages.

# 3, added t2 free = FALSE,  if FALSE, the cost of using T3 and T2 testers will be accounted seperately

# modified at 26.07.2016 by Jose Marulanda to include one more stage of phenotypic selection, so, up to three (instead of two) stages of phenotypic selecion
# all chages marked as SP, because it was made for the Sweet Potato project
#

# modified at 23.08.2016, by JM: debugged: Allow nursery selection in schemes using MAS. Changes labeled as #JM_1

#


multistageoptimum.searchThreeS<-function (maseff=0.4, VGCAandE,
  VSCA=c(0,0,0,0), CostProd = c(0.5,1,1,1), CostTest = c(0.5,1,1,1), #SP
  Nf = 10, Budget = 10021, N2grid = c(11, 1211, 30, 20), #SP
  N3grid = c(11, 211, 5, 3),N4grid = c(11, 211, 5, 3), L2grid=c(1,3,1), L3grid=c(6,8,1), L4grid=c(6,8,1),#SP
  T2grid=c(1,1,1), T3grid=c(1,1,1),T4grid=c(1,1,1),R2=1,R3=1,R4=1,  alg = Miwa(),detail=FALSE,fig=FALSE, #SP
  alpha.nursery=1,cost.nursery=c(0,0) #JM
  ,t2free= FALSE, # Mi 2015-11-10
  parallel.search=FALSE
  )

{
  if (parallel.search)
  {
    no_cores <- detectCores() - 1
    #no_cores <- 1
    cl <- makeCluster(no_cores);
  #  clusterEvalQ(cl,library(selectiongain))
  # just using for testing grid2
    clusterEvalQ(cl,{library(selectiongain);source("multistageoptimum.grid.R")})
  }

#  source('H:/2015-Qingdao/2015-08-15-selectiongain/2_selectiongain_2.0.40_Modif_JM/selectiongainv49/R/multistageoptimum.grid.R')

# pre-define parameters
  Vgca=VGCAandE
  Vsca=VSCA

  L2limit=L2grid
  L3limit=L3grid
  L4limit=L4grid #SP
  T2limit=T2grid
  T3limit=T3grid
  T4limit=T4grid #SP

  alpha.nur=alpha.nursery # JM
  Cost.nur=cost.nursery # JM

  dim=(T4limit[2]-T4limit[1]+1)/T4limit[3]*(T3limit[2]-T3limit[1]+1)/T3limit[3]*(T2limit[2]-T2limit[1]+1)/T2limit[3]*(L4limit[2]-L4limit[1]+1)/L4limit[3]*(L3limit[2]-L3limit[1]+1)/L3limit[3]*(L2limit[2]-L2limit[1]+1)/L2limit[3] #SP

  gainmatrix=array(0,c(1,23)) # JM  # modifyed by X. Mi v49 #SP
  colnames(gainmatrix)<-c("Nf","Nini","alpha.nur","N1","N2","N3","N4","L2","L3","L4","T2","T3","T4","R2","R3","R4","Bini","B1","B2","B3","B4","Budget","Gain") # JM  # modifyed by X. Mi v49 #SP


# main function

    if (alpha.nur == 1 )                                           # JM
       {                                                           # JM
         CostProdMod1<-CostProd[1]+Cost.nur[1]                     # JM
	     warning("No nursery is used as alpha.nursery is set to 1. Then cost of production in Nursery added to CostProd[1]")  # JM
       }
    else
	   {
	   CostProdMod1<-CostProd[1]
	   }


  if (Budget< c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1]+N4grid[1]*L4grid[1]*T4grid[1])) #SP
  {
    warning("budget too small, try value => c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1]+N4grid[1]*L4grid[1]*T4grid[1])") #SP
  }

   # modified at 2015-11-09, from stop to warning, make more freedom

  if (Budget> c(N2grid[2]*L2grid[2]*T2grid[2]+N3grid[2]*L3grid[2]*T3grid[2]+N4grid[2]*L4grid[2]*T4grid[2]))#SP
  {
    warning("budget too great, try value =>  c(N2grid[2]*L2grid[2]*T2grid[2]+N3grid[2]*L3grid[2]*T3grid[2]+N4grid[2]*L4grid[2]*T4grid[2])")#SP
  }

  # modified at 2015-11-09, from stop to warning, make more freedom

# if (length(CostTest)!= 4) #SP
#  {
#    stop( "dimension of CostTest has to be 4") #SP
#  }
#  if (length(CostProd)!= 4) #SP
#  {
#    stop( "dimension of CostProd has to be 4")#SP
#  }

# modifyed at 2016-10-05, stopped otherwise the example of crop sci can not run

for (T4 in seq.int(T4limit[1],T4limit[2],T4limit[3])) #SP
{ #SP
for (T3 in seq.int(T3limit[1],T3limit[2],T3limit[3]))
{
  for (T2 in seq.int(T2limit[1],T2limit[2],T2limit[3]))
	{
    for (L4 in L4limit[1]:L4limit[2]) #SP
      # seq(T3limit[1],T3limit[2],T3limit[3])
    {
    for (L3 in L3limit[1]:L3limit[2])
	  # seq(T3limit[1],T3limit[2],T3limit[3])
    {
		   for (L2 in L2limit[1]:L2limit[2])
		   #seq(T3limit[1],T3limit[2],T3limit[3])

			{
			      allocation = c(Nf,0,alpha.nur,0,0,0,0,L2,L3,L4,T2,T3,T4,R2,R3,R4,0,0,0,0,0,0,0) # JM # modifyed by X. Mi v49 #SP
            gainmatrix=rbind(gainmatrix,allocation)
			}
		}
	}
}
}
}

# here begins the loop circle


theloop<-function(j,gainmatrix, N2grid= N2grid, N3grid= N3grid,N4grid= N4grid,maseff,t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=c(0,0) ) #SP

{

i=j+1

#alpha.nur=alpha.nursery
Cost.nur=cost.nursery

alpha.nur=gainmatrix[i,"alpha.nur"]

#BudgetDH=gainmatrix[i,"U-DH"]
#BudgetMAS=10000-BudgetDH
L4=gainmatrix[i,"L4"] #SP
L3=gainmatrix[i,"L3"]
L2=gainmatrix[i,"L2"]

T4=gainmatrix[i,"T4"] #SP
T3=gainmatrix[i,"T3"]
T2=gainmatrix[i,"T2"]

R4=gainmatrix[i,"R4"] #SP
R3=gainmatrix[i,"R3"]
R2=gainmatrix[i,"R2"]

N.fs=gainmatrix[i,"Nf"]

L1=1
T1=1
R1=1

if (!is.na(maseff))
{

corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3,L4),Rep=c(R1,R2,R3,R4),T=c(T1,T2,T3,T4),index=FALSE,maseff=maseff) #SP
corr.matrix=corr.longin.mas.index

# to temperaly fix the problem of dim for CostTestloop =4
# I deleted the conditionals as they are not longer needed because we created a new function multistageoptimum.searchThreeS

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3,CostTest[4]*L4*T4*R4) #SP


# CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM

if(t2free)
  {
    CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2),CostProd[4]*(T4-T3))  # JM #SP
  }else
  {
    CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*T3,CostProd[4]*T4) # Mi #SP
  }


result=multistageoptimum.grid(N.upper = c(100000,N2grid[2],N3grid[2],N4grid[2]), N.lower =  c(1,N2grid[1],N3grid[1],N4grid[1]), #SP
                                Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3],N4grid[3]), #SP
                                Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop, Nf = Nf, #SP
                                detail = FALSE, alg = Miwa(),fig=FALSE # JM
                                ,alpha.nursery = alpha.nur,cost.nursery = Cost.nur  #JM_1 the old line was: ,alpha.nursery = 1,cost.nursery = c(0,0) for this reason the nursery stage was not computed.
  )


gainmatrix[i,"Budget"]=Budget
gainmatrix[i,"Nini"]= result[1]                               # JM
gainmatrix[i,"Bini"]= result[1]*(Cost.nur[1]+Cost.nur[2])     # JM
gainmatrix[i,"N1"]= result[2]                                 # JM
gainmatrix[i,"B1"]= result[2]*(CostTest[1]+CostProdMod1)      # JM
gainmatrix[i,"N2"]= result[3]                                 # JM
gainmatrix[i,"N3"]= result[4]                                 # JM
gainmatrix[i,"N4"]= result[5]                                 #SP

gainmatrix[i,"B2"]= result[3]*( (L2*T2*R2*CostTest[2])+(CostProd[2]*T2)) # JM

if(t2free)
 {
  gainmatrix[i,"B3"]= result[4]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*(T3-T2))) # JM
  gainmatrix[i,"B4"]= result[5]*( (L4*T4*R4*CostTest[4])+(CostProd[4]*(T4-T3))) #SP
 }else
 {
  gainmatrix[i,"B3"]= result[4]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*T3))# Mi
  gainmatrix[i,"B4"]= result[5]*( (L4*T4*R4*CostTest[4])+(CostProd[4]*T4)) #SP
 }

gainmatrix[i,"Gain"]= result["Gain"] # 2016-10-06 modifyed by mi, to make consistant
}else

# aqui empieza el computo de sin marcadores moleculares #SP
{
corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L2,L3,L4),Rep=c(R2,R3,R4),T=c(T2,T3,T4),index=FALSE,maseff=maseff) #SP
corr.matrix=corr.longin.mas.index


CostTestloop=c(CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3,CostTest[4]*L4*T4*R4) #SP
#CostProdloop=c(CostProd[1]+(CostProd[2]*T2),CostProd[3]*(T3-T2))  # JM


  if(t2free)
  {
    CostProdloop=c(CostProd[1]+(CostProd[2]*T2),CostProd[3]*(T3-T2),CostProd[4]*(T4-T3))  # JM #SP #JM:Please Keep CostProd[1], otherwise will be a major bug in the nursery
  }else
  {
    CostProdloop=c(CostProd[1]+(CostProd[2]*T2),CostProd[3]*T3,CostProd[4]*T4) # Mi #SP #JM:Please Keep CostProd[1], otherwise will be a major bug in the nursery
  }


result=multistageoptimum.grid( N.upper = c(N2grid[2],N3grid[2],N4grid[2]), N.lower =  c(N2grid[1],N3grid[1],N4grid[1]), #SP
                               Vg=Vgca[1],corr = corr.matrix, width =  c(N2grid[3],N3grid[3],N4grid[3]), #SP
							   Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop,
							   Nf = Nf, detail = FALSE, alg = Miwa(),fig=FALSE  # JM
							   ,alpha.nursery = alpha.nur,cost.nursery = Cost.nur
							   )

							   
gainmatrix[i,"Budget"]=Budget
gainmatrix[i,"Nini"]= result[1]                              # JM
gainmatrix[i,"Bini"]= result[1]*(Cost.nur[1]+Cost.nur[2])    # JM
#gainmatrix[i,"N1"]= result[1]
#gainmatrix[i,"B1"]= result[1]*(CostTest[1]+CostProd[1])
gainmatrix[i,"N2"]= result[2]                                # JM
gainmatrix[i,"N3"]= result[3]                                # JM
gainmatrix[i,"N4"]= result[4]

gainmatrix[i,"B2"]= result[2]*( (L2*T2*R2*CostTest[2])+(CostProdMod1+(CostProd[2]*T2)))# JM


if(t2free)
{
  gainmatrix[i,"B3"]= result[3]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*(T3-T2))) # JM
  gainmatrix[i,"B4"]= result[4]*( (L4*T4*R4*CostTest[4])+(CostProd[4]*(T4-T3)))  #SP
}else
{
  gainmatrix[i,"B3"]= result[3]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*T3))# Mi
  gainmatrix[i,"B4"]= result[4]*( (L4*T4*R4*CostTest[4])+(CostProd[4]*T4)) #SP
}


gainmatrix[i,"Gain"]= result["Gain"]                                     #2016-10-06 modifyed by mi to keep consistant

}

gainmatrix[i,]

}

if (!parallel.search)
{

  for (j in 1:dim )
  {
    gainmatrix[j+1,]= theloop(j=j,gainmatrix,N2grid= N2grid, N3grid= N3grid,N4grid= N4grid,maseff,t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=cost.nursery) #SP

  }
}else if(parallel.search)
{
  #clusterExport(cl, "alpha.nursery")
   resulta<- parSapply(cl=cl, 1:dim, FUN=theloop,gainmatrix,N2grid= N2grid, N3grid= N3grid,N4grid= N4grid,maseff=maseff,t2free=t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=cost.nursery) #SP
   gainmatrix[1:dim+1,]<-t(resulta)
}




Output=round( gainmatrix,digits=1)
Output[,"Gain"]=round( gainmatrix[,"Gain"],digits=3)
output=Output[2:c(dim+1),c("Nf", "Nini","N1", "N2", "N3","N4","L2","L3","L4","T2","T3","T4","R2","R3","R4","Bini","B1","B2","B3","B4","Budget","Gain")]  # JM #SP

gainmax=max(output[,"Gain"])

gainlocation = which(output[,"Gain"]==gainmax,arr.ind =TRUE)


if (fig==TRUE)
{
  i=gainlocation[1]
#BudgetDH=gainmatrix[i,"U-DH"]
#BudgetMAS=10000-BudgetDH
L3=gainmatrix[i,"L3"]
L2=gainmatrix[i,"L2"]

T3=gainmatrix[i,"T3"]
T2=gainmatrix[i,"T2"]
R3=gainmatrix[i,"R3"]
R2=gainmatrix[i,"R2"]
N.fs=gainmatrix[i,"Nf"]

L1=1
T1=1
R1=1

corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3),Rep=c(R1,R2,R3),T=c(T1,T2,T3),index=FALSE,maseff=maseff)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)
# CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM

if(t2free)
{
  CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM
}else
{
  CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*T3) # Mi
}

result=multistageoptimum.grid( N.upper = c(100000,N2grid[2],N3grid[2]), N.lower =  c(1,N2grid[1],N3grid[1]),
                               Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3]),
							   Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop, Nf = Nf,
							   detail = detail, alg = Miwa(),fig=TRUE # JM
							   ,alpha.nursery=alpha.nursery,cost.nursery=cost.nursery)
  }

  if (detail==TRUE)
  {
     output
  }else  if (detail!=TRUE )
  {
     output[gainlocation[1],]
	}

}
