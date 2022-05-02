####################################################
# A simple script to calculate activities of TiO2, SiO2 and fO2 for perpleX 
# Original MATLAB script by Jesse Walters <jesse.walters@maine.edu>
# Translated to R by J.-F. Moyen <jfmoyen@gmail.com>
####################################################

## Calculate a(TiO2) from mu(TiO2) and a(SiO2) from mu(SiO2)
# We need to solve
# aTiO2 = exp( -(G_ru - mu_TiO2)/ RT)
# G_ru = G0 + RTln(a_ru)

cat("This works in HSC convention, i.e.\n",
    "The tag HSC_conversion mmust be commented OFF in hpxxxver.dat:\n",
    "| HSC_conversion \n"
    )

############ FROM THE DATABASES ###############
# Thermodynamics of ru
# PPx
# ru        EoS = 8 | H=  -944160.0
# TiO2(1)
# GH = -959216.6  S0 = 50.5  V0 = 1.882  
# c1 = 90.4  c2 = .29E-2  c5 = -623.8  
# b1 = .224E-4  b5 = 457.0037  b6 = 2220000  b7 = -.19E-5  b8 = 4.24  
# end

# TC
# ru   1  2   1.0000 10   2.0000  0
# -944.36   0.05050   1.8820
# 0.0904    0.000002900       0.0   -0.6238
# 0.0000224  2220.00     4.24   -0.00190    0


# Thermodynamics of qz
# q         EoS = 8 | H=  -910710.0
# SiO2(1)
# GH = -923062.4  S0 = 41.43  V0 = 2.269  
# c1 = 92.9  c2 = -.642E-3  c3 = -714900  c5 = -716.1  
# b5 = 525.2346  b6 = 730000  b7 = -.82E-5  b8 = 6  
# transition = 1  type = 4  t1 = 847  t2 = 4.95  t3 = .1188  
# end

# TC
# q   1  1   1.0000 10   2.0000  0
# -910.71   0.04143   2.2690
# 0.0929   -0.000000642    -714.9   -0.7161
# 0.0000000   730.00     6.00   -0.00820    1    847   0.00495    0.1188
# ? ? n_at ? n_O ? 
#   Hf S0 v0
# c1 c2 c3 c5  = a b c d
# b1 b6 b8 b7 ? t1 t2 ? = alpha K0 dK0 d2K0 Tc0 Smax Vmax ?

# name	?	?	n_cat	?	n_O	?	
#   Hf	S0	v0					
# a	b	c	d				
# alpha	K0	dK0	d2K0	?	Tc0	Smax	Vmax

# Jesse W.
# -910720.00000000000000	41.43000000000000	0.00002269000000	
# 92.90000000000000	-0.00064200000000	-714900.00000000000000	-716.10000000000000	
# 0.00000000000000 73000000000.00000000000000	6.00000000000000	-0.00000000008200	
# 847.00000000000000	4.95000000000000	0.00000118800000	3.00000000000000

# Hf=Mt_thermo(:,1); %enthalpies of fortmation (J)
# S0=Mt_thermo(:,2); %entropies at 1 bar 298K (J/(K*mol))
# V0=Mt_thermo(:,3); %volume at 1 bar 298K (m^3/mol)

# %Heat capacity terms
# a=Mt_thermo(:,4);
# b=Mt_thermo(:,5);
# c=Mt_thermo(:,6);
# d=Mt_thermo(:,7);

# alpha=Mt_thermo(:,8); %thermal expansivity (1/K)
# K0=Mt_thermo(:,9); %bulk modulus (Pa)
# dK0=Mt_thermo(:,10); %first derv. bulk modulus
# d2K0=Mt_thermo(:,11); %second derv. bulk modulus (1/Pa)

# Tc0=Mt_thermo(:,12); %Temperature of the landau transition (K)
# Smax=Mt_thermo(:,13); %Max Entropy (J/(K*mol))
# Vmax=Mt_thermo(:,14); %Max volume (Pa)
# n=Mt_thermo(:,15); %number of atoms per formula unit 


# Define minerals lazily, from above (may also be used to read from file, one day)

# Order of values in TC database
tc_format<-c("V1","V2","nAt","V3","nOx","V4",
             "Hf","S0","V0",
             "aCp","bCp","cCp","dCp",
             "alpha","K0","dK0","d2K0","V5","Tc0","Smax","Vmax")

# Conversion factors, Thermocalc to SI
to_SI<-c(1,	1,	1,	1,	1,	1,	
         1.00E+03,	1.00E+03,	1.00E-05,					
         1.00E+03,	1.00E+03,	1.00E+03,	1.00E+03,				
         1.00E+00,	1.00E+08,	1.00E+00,	1.00E-08,	1.00E+00,	1.00E+00,	1.00E+03,	1.00E-05)

# From ds-62.txt
Qtz <- c( 1,  1 ,  1.0000, 10  , 2.0000 , 0,
  -910.71,   0.04143,   2.2690,
   0.0929,   -0.000000642,    -714.9,  -0.7161,
   0.0000000,   730.00,     6.00,   -0.00820 ,   1 ,   847 ,  0.00495,    0.1188)

Rt  <-c(   1,  2,   1.0000, 10,   2.0000,  0,
     -944.36 ,  0.05050,   1.8820,
     0.0904  ,  0.000002900   ,    0.0 ,  -0.6238,
    0.0000224,  2220.00  ,   4.24  , -0.00190,    0, 0, 0, 0 )

Qtz <- Qtz * to_SI
Rt  <- Rt  * to_SI

names(Qtz)<-tc_format
names(Rt)<-tc_format

R <- 0.0083143*1000 # Gas constant (J/(K*mol))

TEOS<-function(V0,alpha,K0,dK0,d2K0,Pbar,TK,S0,n){
  # THERMOCALC 2011 (Tait) EOS FOR SOLIDS
  # Adapted from matlab code by Jesse Walters
  
  aT   <- (1+dK0)/(1+dK0+K0*d2K0)
  bT   <- (dK0/K0)-(d2K0/(1+dK0))
  cnum <- 1+dK0+K0*d2K0
  cden <- ((dK0*dK0)+dK0)-(K0*d2K0)
  cT   <- cnum/cden
  
  # Einstein Temperature approximation
  theta <- 10636/(S0/n+6.44)
  x <- theta/TK
  x0 <- theta/298.15
  ex <- exp(x)
  ex0 <- exp(x0)
  
  # Einstein thermal energy 
  Eth  <- 3*n*8.3144598*theta*(0.5+(1/(ex-1)))
  Eth0 <- 3*n*8.3144598*theta*(0.5+(1/(ex0-1)))
  
  # Einstein heat capacity
  Cv0 <- 3*n*8.3144598*((x0*x0*ex0)/((ex0-1)*(ex0-1)))
  
  # thermal pressures
  Pt <- (alpha*K0*Eth)/Cv0
  Pt0 <- ((alpha*K0*Eth0)/Cv0)
  
  # thermal pressure relative to standard state
  Pth <- Pt-Pt0
  
  # intVdP
  P0 <- 100000;
  psubpth <- Pbar-P0-Pth

  intVdP <- (Pbar-P0)*V0*(1-aT+(aT*(((1-bT*Pth)^(1-cT))-((1+bT*(psubpth))^(1-cT)))/(bT*(cT-1)*(Pbar-P0))))
  # cat("TEOS:",intVdP,"\n",
  #     "n:",n,"\n",
  #     "Cv0:",Cv0,"\n",
  #     "Pth:",Pth,"\n"
  #     )
  return(intVdP)
}

Landau <- function(Tc0,Vmax,Smax,Pbar,TK){
  # Landau Model - Holland and Powell 2011
  
  # note: the HP2011 database uses a correcterd Landau model for some
  # minerals. This correction was not incorporated into the HP2011 paper;
  # however, this is code uses the same formulism as Thermocalc. For the
  # original HP1998 approach, please see the other code. In addition, thermal
  # expansivity and bulk modulus are not included. In practice, HP set the Vf
  # term of the volume integral to 1.
  
  # Adapted from matlab code by Jesse Walters
  
  P0 <- 1e5 # pressure of experimental data in pascals
  T0 <- 298.15 # temperature of experimental data
  
  # Critical Temperature
  
  Tc <- Tc0+(Vmax*(Pbar-P0))/Smax # critical T of landau transition at P of interest
  
  # Q - Order Parameter
  # When T>Tc then Q=0
  if(TK>Tc){ # for T's greater than the critical T
    Q <- 0
  } else {
    Q <- ((Tc-TK)/Tc0)^(0.25)
  }
  if(T0>Tc0){# for T's greater than the critical T 
    # CHANGED Tc -> Tc0 ???
    Q0 <- 0
  } else {
    Q0 <- ((Tc0-T0)/Tc0)^(1/4) # Q at reference T (298K)
  }

  # Excess Gibbs free energies
  Term1 <- Tc0*Smax*((Q0^2)-(1/3)*(Q0^6))
  Term2 <- Smax*(Tc*Q*Q-Tc0*(Q^6)*(1/3))
  Term3 <- TK*Smax*(Q0*Q0-Q*Q)
  Term4 <- (Pbar-P0)*Vmax*Q0*Q0
  
  # To deal with minerals that do not have a Landau correction, the Gibbs free
  # energy contribution of the Landau correction is set to zero for every row
  # with a critical temperature of zero. This requires the input data to have
  # zeros even if a critical temperature is not reported. 
  
  if(Tc0==0){
    Gl <- 0
  }else{
    Gl <- Term1-Term2-Term3+Term4
  }

  
  dGdT <- Smax*(Q*Q-Q0*Q0)  # molar entropy contribution
  dGdP <- -Vmax*(Q*Q-Q0*Q0) # molar volume contribution
  d2GdP2 <- -(Vmax^2)/(2*Smax*Tc0*Q*Q) # bulk modulus contribution
  d2GdT2 <- -Smax/(2*Tc0*Q*Q) # heat capacity contribution
  d2GdPdT <- Vmax/(2*Tc0*Q*Q) # thermal expansivity contribution
  
 
  return(list(Gl=Gl,
              dGdT=dGdT,
              dGdP=dGdP,
              d2GdP2=d2GdP2,
              d2GdT2=d2GdT2,
              d2GdPdT=d2GdPdT))
}

GibbsSolid_landau<-function(Hf,S0,V0,a,b,c,d,alpha,K0,dK0,d2K0,PGPa,TK,n,Vmax,Smax,Tc0,doLandau){
  
  # This function calculates the molar gibbs free energy for solid phases in
  # two parts: G=G0+RTlnK
  # where
  # G0=Hf-TS0-integral(Cp)dT-T*integral(Cp/T)dT+integral(Vsolid)dP
  #
  # Adapted from matlab code by Jesse Walters
  
  ## Integral(Cp)dT
  T0 <- 298.15
  Pbar <- PGPa*1e9 # converts GPa to pascals
  intCpdT <- (a*TK+0.5*b*TK*TK-c/TK+2*d*sqrt(TK))-(a*T0+0.5*b*T0*T0-c/T0+2.0*d*sqrt(T0))
  
  ## Integral(Cp/T)dT
  intCpoverTdT <- a*log(TK/298.15)+b*(TK-298.15)-(c/2)*(1/(TK*TK)-1/(298.15*298.15))-2*d*(1/sqrt(TK)-1/(sqrt(298.15)))
  
  # calls the TEOS to calculate the
  # volume contribution to the Gibbs free energy 
  intVdP <- TEOS(V0,alpha,K0,dK0,d2K0,Pbar,TK,S0,n) 
  
  # excess gibbs free energy from Landau
  ### ONLY if the mineral has a Landau transition !
  if(doLandau){
    Gl <- Landau(Tc0,Vmax,Smax,Pbar,TK)$Gl 
  }else{Gl <-0}
  
  Gs <- Hf+intCpdT-TK*(S0+intCpoverTdT)+intVdP+Gl;

  return(Gs)
}

Gibbs<-function(mineral,PGPa,TK){
  # Wrapper
  G <- GibbsSolid_landau(mineral["Hf"],
                    mineral["S0"],
                    mineral["V0"],
                    mineral["aCp"],
                    mineral["bCp"],
                    mineral["cCp"],
                    mineral["dCp"],
                    mineral["alpha"],
                    mineral["K0"],
                    mineral["dK0"],
                    mineral["d2K0"],
                    PGPa,
                    TK,
                    mineral["nAt"]+mineral["nOx"],
                    mineral["Vmax"],
                    mineral["Smax"],
                    mineral["Tc0"],
                    doLandau=mineral["Tc0"]>0)
  
  names(G)<-"G"
  return(G)
}

vGibbs <- Vectorize(Gibbs,c("PGPa","TK"))


###### fO2 #####
## O2 thermo

O2thermo<-c(0.00E+00,205.2,48.3,-0.000691,499200,-420.7,
            54.5963,-8.6392,0.918301,-3305.58,0.00230524,
            0.000693054,-8.38E-05,0,100000)
names(O2thermo)<-c("Hf", # enthalpies of fortmation (J)
                   "S0", # entropies at 1 bar 298K (J/(K*mol))
                   "a", # heat capacity terms
                   "b",
                   "c",
                   "d",
                   "Ca0", # cork parameters
                   "Ca1",
                   "Cb0",
                   "Cc0",
                   "Cc1",
                   "Cd0",
                   "Cd1",
                   "Tc", # critical T (K) 
                   "Pc"  # critical P (Pa)
)

## Compensated-Redlich-Kwong (CORK) Equation
## Following Holland and Powell 1991 
# Jesse Walters
# Rtlnf

CORK<-function(Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf,Pbar,TK){
  ## Equations 9
  
  Ca <- ((Ca0*(Tcf^(2.5)))/Pcf)+((Ca1*(TK*Tcf^(1.5)))/(Pcf))
  Cb <- (Cb0*Tcf)/Pcf
  Cc <- ((Cc0*Tcf)/(Pcf^(1.5)))+((TK*Cc1)/((Pcf^(1.5))))
  Cd <- ((Cd0*Tcf)/(Pcf*Pcf))+((TK*Cd1)/(Pcf*Pcf))
  
  ## Fugacity equation
  Pr <- Pbar-10000 # relative pressure (Pa)
  
  if(Tcf==0){RTlnf <- 0}else{
    RTlnf <- R*TK*log(1e-5*Pr)+(Cb*Pr)+Ca/(Cb*sqrt(TK))*(log((R*TK)+(Cb*Pr))-log((R*TK)+(2*Cb*Pr)))+((2/3)*Cc*Pr*sqrt(Pr))+((Cd/2)*Pr*Pr)
  }
  
  return(RTlnf)
}


## Molar Gibbs Free Energy Calculator for Pure Solid Phases
# Jesse Walters
#Gf
GibbsPure<-function(TK,PGPa,Hff,S0f,af,bf,cf,df,Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf)
{
  # This function calculates the molar gibbs free energy for solid phases in
  # two parts: G=G0+RTlnK
  # where
  # G0=Hf-TS0-integral(Cp)dT-T*integral(Cp/T)dT+integral(Vsolid)dP
  
  Pbar <- PGPa*1e9 # converts GPa to pascals
  
 ## Integral(Cp)dT
  T0 <- 298.15
  intCpdTf <- (af*TK+0.5*bf*TK*TK-cf/TK+2*df*sqrt(TK))-
              (af*T0+0.5*bf*T0*T0-cf/T0+2*df*sqrt(T0))
  
 ## Integral(Cp/T)dT
  intCpoverTdTf <- af*log(TK/298.15)+bf*(TK-298.15)-(cf/2)*(1/(TK*TK)-
              1/(298.15*298.15))-2*df*(1/sqrt(TK)-1/(sqrt(298.15)))
  
  ## calls the CORK EoS to calculate the volume contribution to the Gibbs free energy 
  RTlnf <- CORK(Ca0,Ca1,Cb0,Cc0,Cc1,Cd0,Cd1,Tcf,Pcf,Pbar,TK);
  
  Gf<- Hff + intCpdTf - TK*(S0f+intCpoverTdTf)+RTlnf
  
  return(Gf)
  
}

# Wrapper
GibbsO2<-function(PGPa,TK){
  # Wrapper
  G <- GibbsPure(TK,PGPa,
                 Hff=O2thermo["Hf"],
                 S0f=O2thermo["S0"],
                 af=O2thermo["a"],
                 bf=O2thermo["b"],
                 cf=O2thermo["c"],
                 df=O2thermo["d"],
                 Ca0=O2thermo["Ca0"],
                 Ca1=O2thermo["Ca1"],
                 Cb0=O2thermo["Cb0"],
                 Cc0=O2thermo["Cc0"],
                 Cc1=O2thermo["Cc1"],
                 Cd0=O2thermo["Cd0"],
                 Cd1=O2thermo["Cd1"],
                 Tcf=O2thermo["Tc"],
                 Pcf=O2thermo["Pc"]
  )
    
  
  names(G)<-"G"
  return(G)
}

fO2<-function(mu,PGpa,TK){
  G0<-GibbsO2(0.0001,TK)
  return( exp((mu - G0) / R / TK )  )
}

FMQ<-function(Pbar,TK){
  A <- 5.5976
  B <- 24505
  C <- 0.8099
  D <- 0.0937
  
  logfO2<-A - B / TK + C * log10(TK) + D * (Pbar - 1) / TK
  return(logfO2)
}


# Example of usage

 fO2(-500000,0.3,973)
 Gibbs(Qtz,0.5,873)
 Gibbs(Qtz,0.5,973)

