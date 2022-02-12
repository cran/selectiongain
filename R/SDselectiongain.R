# SDSG.r author: jose marulanda, first version 18-08-2015, for selectiongain package v2.0.47

# Function to compute the Standard deviation of selection gain for three stage selection
# Using asumptions mentioned by Longin et al 2015 TAG "Genomic selection in wheat: optimum..."
# CARE SHOULD BE TAKEN WITH REGARD TO VARIANCE COMPONENTS AND HERITABILITY ESTIMATION!
# Author: Jose Marulanda


SDselectiongain<-function (Ob, 
                           maseff=NA, 
                           VGCAandE=c(0,0,0,0,0), 
                           VSCA=c(0,0,0,0), 
                           VLine=c(0,0,0,0,0), 
                           years=1,
                           Genotypes="Hybrids")
  
{
  if (class(Ob)[1] =="numeric") 
  {Ob = data.frame(t(Ob))}
  
   # create two new variables corresponding to:
  objectnew<-data.frame(Ob,SDsg=rep(0,dim(Ob)[1]),AnualG=rep(0,dim(Ob)[1]),row.names = c(1:dim(Ob)[1]))
  #standard error of the selection gain
  #GainYear<-c() #Selection gain per year
  
  ####### Hybrid Schemes ###################################################################
  if (Genotypes=="Hybrids")
  {
    
    if (any(colnames(Ob) == "N3"))
    {
      ##########################################################################################
      # Standard deviation of selectiongain for two stage scheme without GS: HybridPSstandard ####
      if(sum(Ob[,"B3"])!=0 && is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N3"]/Ob[g,"N2"]
          alpha2=Ob[g,"Nf"]/Ob[g,"N3"]
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          gamma1=i1*(i1-k1)
          m2= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L2"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L2"])+(VSCA[1]/Ob[g,"T2"])+(VSCA[2]/(Ob[g,"T2"]*Ob[g,"L2"]))+    
            (VSCA[3]/Ob[g,"T2"])+(VSCA[4]/(Ob[g,"T2"]*Ob[g,"L2"]))+(VGCAandE[5]/(Ob[g,"L2"]*Ob[g,"R2"]*Ob[g,"T2"]))
          H2stage1=VGCAandE[1]/m2
          rhozx1=sqrt(H2stage1)
          vpg=VGCAandE[1]*(1-rhozx1^2*gamma1)
          m3= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L3"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L3"])+(VSCA[1]/Ob[g,"T3"])+(VSCA[2]/(Ob[g,"T3"]*Ob[g,"L3"]))+    
            (VSCA[3]/Ob[g,"T3"])+(VSCA[4]/(Ob[g,"T3"]*Ob[g,"L3"]))+(VGCAandE[5]/(Ob[g,"L3"]*Ob[g,"R3"]*Ob[g,"T3"]))
          H2stage2=VGCAandE[1]/m3
          rhozx2=sqrt(H2stage2)
          Pc<-min(Ob[g,"L2"],Ob[g,"L3"])
          Tc<-min(Ob[g,"T2"],Ob[g,"T3"])
          rhox1x2=(VGCAandE[1]+(Pc*VGCAandE[2])/(Ob[g,"L2"]*Ob[g,"L3"])+(Tc*VSCA[1])/(Ob[g,"T2"]*Ob[g,"T3"])+(Pc*Tc*VSCA[2])/(Ob[g,"T2"]*Ob[g,"T3"]*Ob[g,"L2"]*Ob[g,"L3"]))/(sqrt(m2)*sqrt(m3))
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          SDsg2<-sqrt((1-(rhopzx2^2)*(1-v2))*(vpg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      ##########################################################################################
      # Standard deviation of selectiongain for one stage scheme without GS: HybridPSrapid ####
      else   if(sum(Ob[,"B3"])==0 && sum(Ob[,"B2"])!=0 && is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"Nf"]/Ob[g,"N2"]
          
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          m2= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L2"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L2"])+(VSCA[1]/Ob[g,"T2"])+(VSCA[2]/(Ob[g,"T2"]*Ob[g,"L2"]))+    
            (VSCA[3]/Ob[g,"T2"])+(VSCA[4]/(Ob[g,"T2"]*Ob[g,"L2"]))+(VGCAandE[5]/(Ob[g,"L2"]*Ob[g,"R2"]*Ob[g,"T2"]))
          H2stage1=VGCAandE[1]/m2
          rhozx1=sqrt(H2stage1)
          SDsg2<-sqrt((1-(rhozx1^2)*(1-v1))*(VGCAandE[1]/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      
      ##########################################################################################
      # Standard deviation of selectiongain for three stage scheme with GS: HybridGSstandard ####
      
      else if(sum(Ob[,"B3"])!=0 && !is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N2"]/Ob[g,"N1"]
          alpha2=Ob[g,"N3"]/Ob[g,"N2"]
          alpha3=Ob[g,"Nf"]/Ob[g,"N3"]
          
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          k3=qnorm(1-alpha3,0,1)
          i3=(dnorm(k3,0,1))/alpha3
          
          
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          v3=(1-(i3*(i3-k3)))+((1-alpha3)*(i3-k3)^2)
          
          
          gamma1=i1*(i1-k1)
          gamma2=i2*(i2-k2)
          gamma3=i3*(i3-k3)
          
          m2= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L2"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L2"])+(VSCA[1]/Ob[g,"T2"])+(VSCA[2]/(Ob[g,"T2"]*Ob[g,"L2"]))+    
            (VSCA[3]/Ob[g,"T2"])+(VSCA[4]/(Ob[g,"T2"]*Ob[g,"L2"]))+(VGCAandE[5]/(Ob[g,"L2"]*Ob[g,"R2"]*Ob[g,"T2"]))
          H2stage2=VGCAandE[1]/m2
          
          m3= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L3"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L3"])+(VSCA[1]/Ob[g,"T3"])+(VSCA[2]/(Ob[g,"T3"]*Ob[g,"L3"]))+    
            (VSCA[3]/Ob[g,"T3"])+(VSCA[4]/(Ob[g,"T3"]*Ob[g,"L3"]))+(VGCAandE[5]/(Ob[g,"L3"]*Ob[g,"R3"]*Ob[g,"T3"]))
          H2stage3=VGCAandE[1]/m3    
          
          rhozx1=maseff
          rhozx2=sqrt(H2stage2)
          rhox1x2=rhozx1*rhozx2
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          vpg=VGCAandE[1]*(1-rhozx1^2*gamma1)
          rhozx3=sqrt(H2stage3)
          rhox1x3=rhozx1*rhozx3
          Pc<-min(Ob[g,"L2"],Ob[g,"L3"])
          Tc<-min(Ob[g,"T2"],Ob[g,"T3"])
          rhox2x3=(VGCAandE[1]+(Pc*VGCAandE[2])/(Ob[g,"L2"]*Ob[g,"L3"])+(Tc*VSCA[1])/(Ob[g,"T2"]*Ob[g,"T3"])+(Pc*Tc*VSCA[2])/(Ob[g,"T2"]*Ob[g,"T3"]*Ob[g,"L2"]*Ob[g,"L3"]))/(sqrt(m2)*sqrt(m3))
          rhopx2x3=(rhox2x3-(rhox1x2*rhox1x3*gamma1))/(sqrt(1-rhox1x2^2*gamma1)*sqrt(1-rhox1x3^2*gamma1))
          rhox3x1=(rhozx3-(rhozx1*rhox1x3*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x3^2*gamma1))   
          rhoppzx3=(rhox3x1-(rhopzx2*rhopx2x3*gamma2))/(sqrt(1-rhopzx2^2*gamma2)*sqrt(1-rhopx2x3^2*gamma2))
          vppg=vpg*(1-rhopzx2^2*gamma2)
          SDsg2<-sqrt((1-(rhoppzx3^2)*(1-v3))*(vppg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      ##########################################################################################
      # Standard deviation of selectiongain for two stage scheme with GS: HybridGSrapid    ####
      
      else if(sum(Ob[,"B3"])==0 && sum(Ob[,"B2"])!=0 && !is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N2"]/Ob[g,"N1"]
          alpha2=Ob[g,"Nf"]/Ob[g,"N2"]
          
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          gamma1=i1*(i1-k1)
          m2= VGCAandE[1]+(VGCAandE[2]/Ob[g,"L2"])+VGCAandE[3]+(VGCAandE[4]/Ob[g,"L2"])+(VSCA[1]/Ob[g,"T2"])+(VSCA[2]/(Ob[g,"T2"]*Ob[g,"L2"]))+    
            (VSCA[3]/Ob[g,"T2"])+(VSCA[4]/(Ob[g,"T2"]*Ob[g,"L2"]))+(VGCAandE[5]/(Ob[g,"L2"]*Ob[g,"R2"]*Ob[g,"T2"]))
          H2stage2=VGCAandE[1]/m2
          rhozx1=maseff
          rhozx2=sqrt(H2stage2)
          rhox1x2=rhozx1*rhozx2
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          vpg=VGCAandE[1]*(1-rhozx1^2*gamma1)
          SDsg2<-sqrt((1-(rhopzx2^2)*(1-v2))*(vpg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
    }  
    ##########################################################################################
    # Standard deviation of selectiongain for one stage scheme with GS: HybridGSonly     ####
    else 
    {
      for (g in c(1:dim(Ob)[1])) 
      {
        alpha1=Ob[g,"N2"]/Ob[g,"N1"]
        k1=qnorm(1-alpha1,0,1)
        i1=(dnorm(k1,0,1))/alpha1
        v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
        rhozx1=maseff
        SDsg2<-sqrt((1-(rhozx1^2)*(1-v1))*(VGCAandE[1]/Ob[g,"N2"]))
        objectnew[g,"SDsg"]<-SDsg2
        AnualGain<-Ob[g,"Gain"]/years
        objectnew[g,"AnualG"]<-AnualGain    
      }
      objectnew
    }
    ##########################################################################################
  }
  
  ####### Lines schemes  ##################################################################
  else if(Genotypes=="Lines")
  {
    if (any(colnames(Ob) == "N3"))
    {
      ##########################################################################################
      # Standard deviation of selectiongain for two stage scheme without GS: LinePSstandard ####
      if(sum(Ob[,"B3"])!=0 && is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N3"]/Ob[g,"N2"]
          alpha2=Ob[g,"Nf"]/Ob[g,"N3"]
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          gamma1=i1*(i1-k1)
          m2= (VLine[1] + (VLine[2]/Ob[g,"L2"]) + VLine[3] + (VLine[4]/Ob[g,"L2"]) + (VLine[5]/(Ob[g,"L2"]*Ob[g,"R2"])))
          H2stage1=VLine[1]/m2
          rhozx1=sqrt(H2stage1)
          vpg=VLine[1]*(1-rhozx1^2*gamma1)
          m3= (VLine[1] + (VLine[2]/Ob[g,"L3"]) + VLine[3] + (VLine[4]/Ob[g,"L3"]) + (VLine[5]/(Ob[g,"L3"]*Ob[g,"R3"])))
          H2stage2=VLine[1]/m3
          rhozx2=sqrt(H2stage2)
          Pc<-min(Ob[g,"L2"],Ob[g,"L3"])
          rhox1x2=(VLine[1]+(Pc*VLine[2])/(Ob[g,"L2"]*Ob[g,"L3"]))/(sqrt(m2)*sqrt(m3))
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          SDsg2<-sqrt((1-(rhopzx2^2)*(1-v2))*(vpg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      ##########################################################################################
      # Standard deviation of selectiongain for one stage scheme without GS: LinePSrapid ####
      else   if(sum(Ob[,"B3"])==0 && sum(Ob[,"B2"])!=0 && is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"Nf"]/Ob[g,"N2"]
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          v1=(1-(i1*(i1-k1)))+((1+alpha1)*(i1-k1)^2)
          m2= (VLine[1] + (VLine[2]/Ob[g,"L2"]) + VLine[3] + (VLine[4]/Ob[g,"L2"]) + (VLine[5]/(Ob[g,"L2"]*Ob[g,"R2"])))
          H2stage1=VLine[1]/m2
          rhozx1=sqrt(H2stage1)
          SDsg2<-sqrt((1-(rhozx1^2)*(1-v1))*(VLine[1]/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      
      ##########################################################################################
      # Standard deviation of selectiongain for three stage scheme with GS: LineGSstandard ####
      
      else if(sum(Ob[,"B3"])!=0 && !is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N2"]/Ob[g,"N1"]
          alpha2=Ob[g,"N3"]/Ob[g,"N2"]
          alpha3=Ob[g,"Nf"]/Ob[g,"N3"]
          
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          k3=qnorm(1-alpha3,0,1)
          i3=(dnorm(k3,0,1))/alpha3
          
          
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          v3=(1-(i3*(i3-k3)))+((1-alpha3)*(i3-k3)^2)
          
          
          gamma1=i1*(i1-k1)
          gamma2=i2*(i2-k2)
          gamma3=i3*(i3-k3)
          
          rhozx1=maseff
          vpg=VLine[1]*(1-rhozx1^2*gamma1)
          m2= (VLine[1] + (VLine[2]/Ob[g,"L2"]) + VLine[3] + (VLine[4]/Ob[g,"L2"]) + (VLine[5]/(Ob[g,"L2"]*Ob[g,"R2"])))
          H2stage2=VLine[1]/m2
          rhozx2=sqrt(H2stage2)
          rhox1x2=rhozx1*rhozx2
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          m3= (VLine[1] + (VLine[2]/Ob[g,"L3"]) + VLine[3] + (VLine[4]/Ob[g,"L3"]) + (VLine[5]/(Ob[g,"L3"]*Ob[g,"R3"])))
          H2stage3=VLine[1]/m3
          rhozx3=sqrt(H2stage3)
          rhox1x3=rhozx1*rhozx3
          Pc<-min(Ob[g,"L2"],Ob[g,"L3"])
          rhox2x3=(VLine[1]+(Pc*VLine[2])/(Ob[g,"L2"]*Ob[g,"L3"]))/(sqrt(m2)*sqrt(m3))
          rhopx2x3=(rhox2x3-(rhox1x2*rhox1x3*gamma1))/(sqrt(1-rhox1x2^2*gamma1)*sqrt(1-rhox1x3^2*gamma1))
          rhox3x1=(rhozx3-(rhozx1*rhox1x3*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x3^2*gamma1))
          vppg=vpg*(1-rhopzx2^2*gamma2)
          rhoppzx3=(rhox3x1-(rhopzx2*rhopx2x3*gamma2))/(sqrt(1-rhopzx2^2*gamma2)*sqrt(1-rhopx2x3^2*gamma2))
          SDsg2<-sqrt((1-(rhoppzx3^2)*(1-v3))*(vppg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
      ##########################################################################################
      # Standard deviation of selectiongain for two stage scheme with GS: LineGSrapid    #######
      
      else if(sum(Ob[,"B3"])==0 && sum(Ob[,"B2"])!=0 && !is.na(maseff))
      {
        for (g in c(1:dim(Ob)[1])) 
        {
          alpha1=Ob[g,"N2"]/Ob[g,"N1"]
          alpha2=Ob[g,"Nf"]/Ob[g,"N2"]
          k1=qnorm(1-alpha1,0,1)
          i1=(dnorm(k1,0,1))/alpha1
          k2=qnorm(1-alpha2,0,1)
          i2=(dnorm(k2,0,1))/alpha2
          v1=(1-(i1*(i1-k1)))+((1-alpha1)*(i1-k1)^2)
          v2=(1-(i2*(i2-k2)))+((1-alpha2)*(i2-k2)^2)
          gamma1=i1*(i1-k1)
          rhozx1=maseff
          vpg=VLine[1]*(1-rhozx1^2*gamma1)
          m2= (VLine[1] + (VLine[2]/Ob[g,"L2"]) + VLine[3] + (VLine[4]/Ob[g,"L2"]) + (VLine[5]/(Ob[g,"L2"]*Ob[g,"R2"])))
          H2stage2=VLine[1]/m2
          rhozx2=sqrt(H2stage2)
          rhox1x2=rhozx1*rhozx2
          rhopzx2=(rhozx2-(rhozx1*rhox1x2*gamma1))/(sqrt(1-rhozx1^2*gamma1)*sqrt(1-rhox1x2^2*gamma1))
          SDsg2<-sqrt((1-(rhopzx2^2)*(1-v2))*(vpg/Ob[g,"Nf"]))
          objectnew[g,"SDsg"]<-SDsg2
          AnualGain<-Ob[g,"Gain"]/years
          objectnew[g,"AnualG"]<-AnualGain    
        }
        objectnew
      }
      ##########################################################################################
      
    }  
    ##########################################################################################
    # Standard deviation of selectiongain for two stage scheme with GS: LineGSonly     ####
    else 
    {
      for (g in c(1:dim(Ob)[1])) 
      {
        alpha1=Ob[g,"N2"]/Ob[g,"N1"]
        k1=qnorm(1-alpha1,0,1)
        i1=(dnorm(k1,0,1))/alpha1
        v1=(1-(i1*(i1-k1)))+((1+alpha1)*(i1-k1)^2)
        rhozx1=maseff
        SDsg2<-sqrt((1-(rhozx1^2)*(1-v1))*(VLine[1]/Ob[g,"N2"]))
        objectnew[g,"SDsg"]<-SDsg2
        AnualGain<-Ob[g,"Gain"]/years
        objectnew[g,"AnualG"]<-AnualGain    
      }
      objectnew
    }
    ##########################################################################################
  }  
  
  
}
# end of function


