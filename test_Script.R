# library(R6)
# library(tidyverse)

source("C:\\Users\\moje4671\\Documents\\Recherche\\Zircon_modelling\\Src\\PetroClass_0.2.R")
setwd("C:/Users/moje4671/Documents/Geol Soft/PerpleX_6.9.0")

## KCG Diorite
KCGdio<-PerplexModel$new("KCG Diorite")
#runWerami36("./GrTypes/","KCGdio") %>% .[["phm"]] %>% KCGdio$readPhm(verbose = F)
KCGdio$readPhm("./GrTypes/KCGdio_1.phm",verbose = F)

KCGdio$countDuplicates()
KCGdio$classifyMineral()
KCGdio$completeTable()
KCGdio$refineEutectic()
KCGdio$plotPhaseMap(Tunit="C")

## KCG granite
KCGgrn<-PerplexModel$new("KCG Granite")
#runWerami36("./GrTypes/","KCGgrn") %>% .[["phm"]] %>% KCGgrn$readPhm(verbose = F)
KCGgrn$readPhm("./GrTypes/KCGgrn_21.phm")


KCGgrn$countDuplicates()
KCGgrn$mergeMinerals(mineral="Ilm(WPH)",outMinName="ilmenite")
KCGgrn$classifyMineral(cutoff=6.6)
KCGgrn$completeTable()
KCGgrn$refineEutectic()
KCGgrn$plotPhaseMap(Tunit="K")

KCGgrn$phm %>% filter(TK>850 & TK < 920) %>% select(all_of(c("TK","Name","wt",KCGgrn$majors))) %>% filter (wt>0 & Name =="ilmenite") %>% print(n=60)

## ACG diorite
ACGdio<-PerplexModel$new("ACG Diorite")
runWerami36("./GrTypes/","ACGdio") %>% .[["phm"]] %>% ACGdio$readPhm(verbose = F)

ACGdio$countDuplicates()
ACGdio$classifyMineral()
ACGdio$completeTable()
ACGdio$refineEutectic()
ACGdio$plotPhaseMap(Tunit="K")

## ACG granite
ACGgrn<-PerplexModel$new("ACG Granite")
#runWerami36("./GrTypes/","ACGgrn") %>% .[["phm"]] %>% ACGgrn$readPhm(verbose = F)
ACGgrn$readPhm("./GrTypes/ACGgrn_3.phm")

# There is really not much rutile... 
ACGgrn$phm %>% filter(Name == "ru" & wt > 0) %>% relocate(c("TK","wt"))
# Affect the marginal amount of ru to ilmenite
ACGgrn$renameMinerals("ru","Ilm(WPH)")
ACGgrn$countDuplicates()
ACGgrn$mergeMinerals(mineral="Ilm(WPH)",outMinName="ilmenite")
ACGgrn$classifyMineral()

ACGgrn$completeTable()
ACGgrn$refineEutectic()
ACGgrn$plotPhaseMap(Tunit="K")

# CPG
CPGgrn<-PerplexModel$new("CPG Granite")
runWerami36("./GrTypes/","CPG") %>% .[["phm"]] %>% CPGgrn$readPhm(verbose = F)

CPGgrn$countDuplicates()
CPGgrn$classifyMineral()
CPGgrn$phm %>% filter(Name == "ru" & wt > 0) %>% relocate(c("TK","wt"))
# Affect the marginal amount of ru to ilmenite
CPGgrn$renameMinerals("ru","Ilm(WPH)")
CPGgrn$countDuplicates()
CPGgrn$mergeMinerals(mineral="Ilm(WPH)",outMinName="ilmenite")
CPGgrn$classifyMineral()
CPGgrn$completeTable()
CPGgrn$refineEutectic()
CPGgrn$plotPhaseMap(Tunit="K")

# MPG
MPGgrn<-PerplexModel$new("MPG Granite")
runWerami36("./GrTypes/","MPG") %>% .[["phm"]] %>% MPGgrn$readPhm(verbose = F)


MPGgrn$phm %>% filter(Name == "ru" & wt > 0) %>% relocate(c("TK","wt"))
# Affect the marginal amount of ru to ilmenite
MPGgrn$renameMinerals(c("ru"),c("Ilm(WPH)"))
MPGgrn$classifyMineral(cutoff=6.6)
MPGgrn$countDuplicates()
MPGgrn$mergeMinerals(mineral="Ilm(WPH)",outMinName="ilmenite")

MPGgrn$completeTable()
MPGgrn$refineEutectic()
MPGgrn$plotPhaseMap(Tunit="K")

# tdj
tdj<-PerplexModel$new("Trondhjemite")
runWerami36("./GrTypes/","tdj") %>% .[["phm"]] %>% tdj$readPhm(verbose = F)

tdj$classifyMineral()
tdj$countDuplicates()
tdj$completeTable()
tdj$refineEutectic()

tdj$plotPhaseMap(Tunit="K")

# tdj
ton<-PerplexModel$new("Tonalite")
runWerami36("./GrTypes/","ton") %>% .[["phm"]] %>% ton$readPhm(verbose = F)

ton$classifyMineral()
ton$countDuplicates()
ton$completeTable()
ton$refineEutectic()

ton$plotPhaseMap(Tunit="K")


##
dabuPhFile <- "C:\\Users\\moje4671\\Documents\\Recherche\\Zircon_modelling\\Dabu\\perplex\\diorite03_1.phm"
dabu <-PerplexModel$new("Dabu diorite, 3% H2O")
dabu$readPhm(dabuPhFile,verbose = F)
dabu$countDuplicates()

dabu$completeTable()
dabu$plotPhaseMap(Tunit="K")