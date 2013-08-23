
# package selectiongain

# modified at 28-06-2013, for 1MAS+2PS


multistageoptimum.search<-function (maseff=0.4, VGCAandE, 
  VSCA, CostProd = c(0.5,1,1), CostTest = c(0.5,1,1), 
  Nf = 10, Budget = 10021, N2grid = c(11, 1211, 30), 
  N3grid = c(11, 211, 5), L2grid=c(1,3,1), L3grid=c(6,8,1),
  T2grid=c(1,2,1), T3grid=c(3,5,1),R2=1,R3=1,  alg = Miwa(),detail=FALSE,fig=FALSE)

{  
	if (Budget< c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1]))
	{
	  stop("budget too small, try value => c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1])")
	}
	
	Vgca=VGCAandE
	Vsca=VSCA
	
	L2limit=L2grid
	L3limit=L3grid
	T2limit=T2grid
	T3limit=T3grid
	
	dim=(T3limit[2]-T3limit[1]+1)/T3limit[3]*(T2limit[2]-T2limit[1]+1)/T2limit[3]*(L3limit[2]-L3limit[1]+1)/L3limit[3]*(L2limit[2]-L2limit[1]+1)/L2limit[3]

  gainmatrix=array(0,c(1,15))
  colnames(gainmatrix)<-c("Nf","N1","N2","N3","L2","L3","T2","T3","R2","R3","B1","B2","B3","Budget","Gain")
 
for (T3 in seq.int(T3limit[1],T3limit[2],T3limit[3]))
{
  for (T2 in seq.int(T2limit[1],T2limit[2],T2limit[3]))
	{
	  for (L3 in L3limit[1]:L3limit[2])
	  # seq(T3limit[1],T3limit[2],T3limit[3])
	  
		{
		   for (L2 in L2limit[1]:L2limit[2])
		   #seq(T3limit[1],T3limit[2],T3limit[3])
		   
			{
			      allocation = c(Nf,0,0,0,L2,L3,T2,T3,R2,R3,0,0,0,0,0)
            gainmatrix=rbind(gainmatrix,allocation)
			} 
		}
	}
}

for (j in 1:dim )
{  
i=j+1
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

corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3),Rep=c(R1,R2,R3),T=c(T1,T2,T3),index=FALSE,maseff=0.4)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)


result=multistageoptimum.grid( N.upper = c(100000,N2grid[2],N3grid[2]), N.lower =  c(1,N2grid[1],N3grid[1]), Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3]), Budget = Budget, CostProd =CostProd, CostTest = CostTestloop, Nf = 10, detail = FALSE, alg = Miwa(),fig=FALSE)
gainmatrix[i,"Budget"]=Budget
gainmatrix[i,"N1"]= result[1]
gainmatrix[i,"N2"]= result[2]
gainmatrix[i,"N3"]= result[3]
gainmatrix[i,"B1"]= result[1]*(CostTest[1]+CostProd[1])
gainmatrix[i,"B2"]= result[2]*( L2*T2*R2*CostTest[2]+CostProd[2])
gainmatrix[i,"B3"]= result[3]*( L3*T3*R3*CostTest[3]+CostProd[3])
gainmatrix[i,"Gain"]= result[4]
}


Output=round( gainmatrix,digits=1)
Output[,"Gain"]=round( gainmatrix[,"Gain"],digits=3)
output=Output[2:c(dim+1),c("Nf", "N1", "N2", "N3","L2","L3","T2","T3","R2","R3","B1","B2","B3","Budget","Gain")]
 
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

corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3),Rep=c(R1,R2,R3),T=c(T1,T2,T3),index=FALSE,maseff=0.4)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)

result=multistageoptimum.grid( N.upper = c(100000,N2grid[2],N3grid[2]), N.lower =  c(1,N2grid[1],N3grid[1]), Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3]), Budget = Budget, CostProd =CostProd, CostTest = CostTestloop, Nf = 10, detail = detail, alg = Miwa(),fig=TRUE)
  }

  if (detail==TRUE)
  {
     output
  }else  if (detail!=TRUE )
  {
     output[gainlocation,]
	}

}
