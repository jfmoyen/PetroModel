#### Nandedkar, Kd for amphibole
NanKd<-function(TK){
  
  elemNand<-c("Zr","Gd","Yb","Th","U")
  
  # Parameters for fit of log(Kd) = a*(10000/T(K))+b
  # Fitted by O. Laurent (2020)
  
  fit<-matrix(c(0.5018,0.4562,0.4683,0.8386,0.9578,
                -4.4388,-3.3665,-3.4920,-8.3997,-9.4929),
              byrow=T,
              ncol=5)
  
  colnames(fit)<-elemNand
  rownames(fit)<-c("a","b")
  
  kds<-10^(fit["a",]*10000/TK + fit["b",])
  names(kds)<-elemNand
  return(kds)
}

#### Kirkland, Kd for zircon
KirkKd<-function(TK){
  kdTh<- exp( -5.27 + 8081/TK)
  kdU <- exp( 0.636 + 4244/TK)
  
  kds<-c(kdTh,kdU)
  names(kds)<-c("Th","U")
  return(kds)
}

#### Claiborne, Kd for zircon
ClaibKd<-function(Tippm){
  elemsClaib<-c("Th","U","Nb","Y","Hf","Ce","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu")
  # Claiborne 2018. Fit Kd = a Ti ^b
  
  fit<-matrix(c(70.1,
                659,
                2.53,
                536,
                23108, # Hf, dubious....
                0.942,
                0.164,
                2.52,
                2.91,
                36.7,
                99.2,
                248,
                508,
                1027,
                1379,
                2223,
                3060,
                # b terms
                -0.979,
                -1.191,
                -0.959,
                -1.125,
                -1.168, # Hf
                -0.628,
                -0.675,
                -0.729,
                -0.473,
                -0.853,
                -0.88,
                -1.03,
                -1.008,
                -1.136,
                -1.077,
                -1.249,
                -1.291
                ),
              byrow=T,
              ncol=17)
  
  colnames(fit)<-elemsClaib
  rownames(fit)<-c("a","b")
  
  kds<-fit["a",]*Tippm^fit["b",]
  names(kds)<-elemsClaib
  return(kds)
  }


#### Claiborne, Kd for zircon (version temperature)
ClaibKdTemp<-function(TK){
  elemsClaib<-c("Th","U","Nb","Y","Hf","Ce","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu")
  # Claiborne 2018. Fit Kd = a exp (b/T)
  
  fit<-matrix(c(0.0126, # Th
                0.0465, # U
                0.0003, #Nb
                0.0036, #Y
                0.0965, #Hf
                0.00005, #Ce, dummy
                0.0001, #Nd
                0.0009, #Sm
                0.0032, #Eu
                0.0005, #Gd
                0.0021, #Tb
                0.0041, #Dy
                0.0038, #Ho
                0.0052, #Er
                0.0086, #Tm
                0.0044, #Yb
                0.0032, #Lu
                # b terms
                6696,
                7167,
                7241,
                9806,
                10206,
                4500, #Ce
                5867,
                6636,
                6026,
                9436,
                9160,
                9090,
                9948,
                10088,
                9990,
                10784,
                11358
  ),
  byrow=T,
  ncol=17)
  
  colnames(fit)<-elemsClaib
  rownames(fit)<-c("a","b")
  
  kds<-fit["a",]*exp(fit["b",]/TK)
  names(kds)<-elemsClaib
  return(kds)
}


#### Trail Ce/Ce*
# Trail et al  2012 GCA
TrailKd<-function(fO2,TK,DLa,DPr){
  lnCeD <- 0.1156 * log(fO2) + 13860/TK - 6.125
  
  KdCe<- exp(lnCeD) * sqrt(DLa*DPr)
  names(KdCe)<-"Ce"
  
  return(KdCe)
  
}



############# USE
kd.this<-kd
if(varKdAmph){
  # Alter the Kd value for amphibole
  modifKd <- NanKd(tprK)
  kd.this["Amph",names(modifKd)]<-modifKd
}

#  if(varKdZirc_Claiborne){
#     modifKd <- ClaibKd(tiz)
#      els<-intersect(names(modifKd),colnames(kd.this))
#      kd.this["Zrn",els]<-modifKd
#  }

if(varKdZirc_Claiborne){
  modifKd <- ClaibKdTemp(tprK)
  els<-intersect(names(modifKd),colnames(kd.this))
  kd.this["Zrn",els]<-modifKd
}


if(varKdZirc_Kirkland){
  # Alter the Kd value for zircon
  modifKd <- KirkKd(tprK)
  kd.this["Zrn",names(modifKd)]<-modifKd
}

if(varKdZirc_Trail){
  modifKd <- TrailKd(fug,tprK,kd.this["Zrn","La"],kd.this["Zrn","Pr"])
  kd.this["Zrn",names(modifKd)]<-modifKd
}

