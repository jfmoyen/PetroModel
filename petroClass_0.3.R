library(R6)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(magrittr) # Needs to be imported explicitly, because https://github.com/tidyverse/magrittr/issues/194 (WtF ?)

## Code not yet packaged
source("U:\\Recherche\\Zircon_modelling\\Src\\activities_0.1.R")

################## Utilities  ###############################

## Utility function - translate the headers of a table
translateHeader = function(tab,Glossary,verbose=T){
  #' Utility function to convert the header of a table from old to new names
  #' @param tab: the table to convert
  #' @param glossary: conversion of perpleX names into "clean" names
  #' A two-column df, with old and matching new names
  #' @return the renamed version of the table
  colnames(Glossary)<-c("oldNames","newNames")
  hdr <- colnames(tab)
  for(i in 1:length(hdr)){
    if(hdr[i] %in% Glossary$oldNames){
      old<-hdr[i]
      new<-as.character(Glossary$newNames[Glossary$oldNames==old])
      hdr[i]<-new
      if(verbose){
        cat(old,"-->",new)
      }
      
    }else{
      if(verbose){
        cat(hdr[i],"not replaced")
      }
      
    }
    if(verbose){
      cat("\n")
    }
    
  }
  
  colnames(tab)<-hdr
  return(tab)
}


# Run any program of the perplex suite
runPerplex <- function(exe,prefix,datfile,com0){
  #' A function to call any of the perpleX executables
  #' 
  #' This is mostly used for side effects (vertex, werami...). Also, note that vertex takes a long time to run, so using this on vertex may not be a good idea.
  #' @param exe: name of the executable (program to run). Full-qualified path from current dir
  #' @param prefix: path from prog dir to file
  #' @param datfile: name of .dat file
  #' @param com0: commands passed to prog
  #' @return The console output of the program. It is also written to a file in the same directory as the .dat. Note that perpleX programs typically create other files, which is the main use of this function !
  
  cat(exe,": processing",datfile,"...")
  exeFull <- paste(exe,"exe",sep=".")
  ee<-system(exeFull,intern=T,
             input=c(paste(prefix,datfile,sep=""),com0))
  
  fn<-paste(prefix,datfile,"_",exe,"Output.txt",sep="")
            
  save(ee,file=fn,ascii=T)
  cat("...done\n")
  return(ee)
}

# Werami mode 36 - assume no complications !
runWerami36<-function(prefix,datfile){
  #' Run werami in mode 36 (all phase and system properties). This generates the phm file that is used by the rest of the class. Requires the working directory to be the one where the executables reside.
  #' 
  #' The function attempts to extract the name of the phm file from the textural output - it should work, but don't use indiscriminately.
  #' @param prefix: path from prog dir to file
  #' @param datfile: name of .dat file
  #' @return (invisibly) a list containing (i) all the text output of werami; (ii) the name of the phm file. Note that a side-effect is to generate the phm file in question ! 
  
  cat("Running Werami mode 36, check for possible complications\n")
  
  wer36com <- c("3","36","3","y","0")
  
  werOutput<-runPerplex("werami",prefix,datfile,wer36com)
  # Attempt to guess the output file name
  # Line containing the phm name
  werOutput %>% str_which(".phm") %>% min -> idxPhmName
  
  # Extract the phm name (naive)
  werOutput[idxPhmName] %>% str_split(":") %>% unlist %>% .[2] %>% str_trim("both") -> phmName
  
  cat("Guessing: the output file is named",phmName,"\n")
  
  invisible(list(output=werOutput,phm=phmName))
}


#  ###### Deprecated version of find duplicates
#  ## Get intervals (i.e. blocks of lines corresponding to one PT point)
#  # As an index, we use the line containing "system", i.e. the first of each block in current (11/2020) phm format
#  i<-grep("system",self$phm[[self$nameCol]])
# 
# # We build a table of intervals describing the phases (NOT system)
#  intervals<-cbind(i[1:length(i)]+1,c(i[2:length(i)],nrow(self$phm))-1)
#  
#  # Get duplicated minerals (ie, those that appear at least twice in at least one block)
#  ee<-sapply(1:nrow(intervals),function(i){
#    min.names<- self$phm[intervals[i,1]:intervals[i,2],self$nameCol]
#    
#    z<-table(min.names)
#    #z<-.countPhases(phm,intervals,i)
#    return(z)
#  })
#  
#  # Clean and build the duplicate table
#  x<-unlist(ee)
#  cat("Duplicates:\n")
#  ( tibble(n=names(x),v=x) %>% unique() %>% filter(v>1) -> self$duplicates ) %>% print
#   
#  invisible(self)



################## PetroPoint ###############################

PetroPoint <- R6Class(classname="PetroPoint", 
                      public=list(
                        
                        initialize = function() {
                          
                        }
                        
                        
                      ))


################## PetroModel ###############################
#' R6 Class Representing a generic model
#'
#' @description
#' The Mother of all models
#'
#' @details
#' Not much to see here now, this is mostly an empty shell

PetroModel <- R6Class(classname="PetroModel", 
                      public=list(
                        
                        #' @field A user-friendly name for the model
                        modName = "Anonymous model",
                        
                        #' @field Names of major elements
                        majors = c("H2O" ,  "MgO"  , "Al2O3", "SiO2" , "K2O" ,  "CaO"  , "FeO" , "O2", "TiO2" , "Na2O" ),
                        
                        #' @field Names of trace elements
                        traces = c("Rb","Sr"),
                        
                        #' @field Number of points in the model
                        nPoints= 0,
                        
                        #' @field Table of (static) partition coefficients
                        #' NOTE that kd is a MATRIX and not a tibble !
                        kd = NULL,
                        
                        #' @field Composition of the source
                        #' NOTE that kd is a VECTOR and not a tibble !
                        c0 = NULL,
                        
                        ## The justification is that they are used in a proper matrix operation during melting...
                        
                        #' @description
                        #' Create a new PetroModel object.
                        #' @param name A human-readable name to describe the model
                        #' @return A new `PetroModel` object.
                        initialize = function(name = "Anonymous model"){
                          self$modName <- name
                          
                          invisible(self)
                          
                        }
                        
                      ))

################## perplexModel ###############################
#' R6 Class Representing a perpleX-generated model
#'
#' @description
#' This is a parser for phm files
#'
#' @details
#' RTFC
#' You should probably do the actions in the following sequence :
#' load -> correct duplicates -> complete the table -> fix eutectic

PerplexModel <- R6Class("PerplexModel", 
                             inherit = PetroModel,
                             public = list(
                               
                               #' @field phm table
                               phm = NULL,
                               
                               ## The following define the characteristics of the phm file
                               
                               #' @field Header rows in the phm file, above the colnames
                               hdrRows=8,
                               
                               #' @field Translation of phm to clean names: original (perpleX) names
                               ppxNames = c("H2O,wt%","MgO,wt%","Al2O3,wt%","SiO2,wt%","K2O,wt%","CaO,wt%","FeO,wt%","O2,wt%","TiO2,wt%","Na2O,wt%",
                                           "mu[TiO2],J/mol","mu[SiO2],J/mol","mu[O2],J/mol"),
                               
                               #' @field Translation of phm to clean names: new (clean) names
                               Names = c("H2O","MgO","Al2O3","SiO2","K2O","CaO","FeO","O2","TiO2","Na2O",
                                           "muTiO2","muSiO2","muO2"),
                               
                               # For now, DO NOT RENAME "H,J/mol","rho,kg/m3","cp,J/K/mol" or bad things will happen (in refineEutectic)
                         
                               # The following are conventional cols, required with their proper name
                               
                               #' @field Colum storing pressure variable - will be "Pbar"
                               PCol = "P(bar)",
                               
                               #' @field Colum storing temperature variable - will be "TK"
                               TCol = "T(K)",
                               
                               #' @field In this column, is T stored as C or K?
                               TisK = T,
                               
                               #' @field Column storing phase name - will be "Name"
                               nameCol = "Name",
                               
                               #' @field Column storing weight proportions phases - will be "wt"
                               propCol = "wt,%",
                               
            
                               #' @field All existing phases
                               phNames = NULL,
                               
                               #' @field Name of the melt phase
                               meltName = "melt(HGP)",
                               
                               #' @field Name of the (ternary) feldspar name
                               feldsparName = "feldspar",
                               
                               #' @field Duplicated phases (and count)
                               duplicates = NULL,
                               
                               #' @field Phase assemblage (and count) for all model points
                               phaseAssemblage = NULL,
                               
                               #' @field These variables are phase proportions
                               propCols = c("vol,%","wt","mol,%"),
                               
                               
                               #' @field All model points and their properties
                               modPoints = NULL,
                               
                               #' @field Column identifying unique model points
                               pointCols = c("TK","Pbar"), ## Counter ?
                               
                               #' @field Global variables, defined on a per-point (not per phase) basis
                               globVar = c("Counter","mu[H2O],J/mol","mu[MgO],J/mol","mu[Al2O3],J/mo","muSiO2",
                                           "mu[K2O],J/mol","mu[CaO],J/mol","mu[FeO],J/mol","muO2","muTiO2",
                                           "mu[Na2O],J/mol","nom_ox"),
                               
                               ## What have we done so far ?
                               #' @field Flag to record whether we did check for duplicates
                               dedup = F,
                               #' @field Flag to record whether we did add points at the eutectic
                               eutectic = F,
                               #' @field Flag to record whether we did rename perpleX names to clean ones
                               rename = F,
                               #' @field Flag to record whether we did complete (tidy up) the table
                               complete = F,
                               
                               # initialize = function() {
                               #   Nothing special for now, we just use the super's
                               #   
                               # }
  
                               # 
                               #' @description Fill a model from a perplex phm file (werami mode 36)
                               #' @param phmFileName (char) Fully-qualified path to the phm file (relative to current directory)
                               #' @param fixnames (boolean) Attempt to fix perpleX names into clean ones, using the dictionary supplied (see fields ppxNames, Names)
                               #' @param verbose (boolean) Blah, blah, blah...
                               #' @param TK, Pbar (numeric) Manual values for temperature and pressure
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self. It contains phm, a tibble that keeps all the info from the phm (werami output) file
                               readPhm = function(phmFileName,fixnames=T,verbose=T,TK=NULL,Pbar=NULL) {
                                 cat("Reading file",phmFileName,"...")
                                 # Read the file
                                 phm<-read.table(phmFileName,skip=self$hdrRows,header=T,stringsAsFactors=F,check.names=F) %>%
                                   as_tibble
                                   
                                 # The tidyverse version is subtly different from the base R version, and does not work as well on phm format
                                 #phm2<-read_table2(phmFileName,skip=self$hdrRows)
                                 
                                 # Convert ppx names to clean names
                                 if(fixnames){
                                   glossary<-data.frame(self$ppxNames,self$Names)
                                   phm <- translateHeader(phm,glossary,verbose=verbose)
                                   self$rename <- T
                                 }
                                 
                                 # Add P or T if missing
                                 # 
                                 if(!is.null(TK)){
                                   phm %<>% bind_cols(manualT=TK)
                                   self$TCol <- "manualT"
                                   cat("Warning, T defined manually, did you correct self$pointCols ?\n")
                                 }
                                 
                                 if(!is.null(Pbar)){
                                   phm %<>% bind_cols(manualP=Pbar)
                                   self$PCol <- "manualP"
                                   cat("Warning, P defined manually, did you correct self$pointCols ?\n")
                                 }
                                 
                                 # Enforce conventional names
                                 convNames<-c("TK","Pbar","Name","wt") # All the functions below expect cols with this name!
                                 orNames <- c(self$TCol,self$PCol,self$nameCol,self$propCol)
                                 glossary<-data.frame(orNames,convNames)
                                 phm <- translateHeader(phm,glossary,verbose=verbose)
                                 
                                 # If needed, convert T to Kelvin
                                 if(!self$TisK){
                                   phm %<>% mutate(TK = TK + 273.15)
                                 }
                                 
                                 # We generate a UniqueId - an identifier for each model point
                                 phm %<>% unite("uniqueId",all_of(self$pointCols),remove=F)
                                 
                                 # And we flag the lines as coming from Werami
                                 phm %<>% add_column(from="werami") 
                             
                                 # Assign to the real thing
                                 self$phm <- phm
                                 
                                 # Build a table of model points
                                 self$modPoints <- phm %>% select(all_of(c("uniqueId",self$pointCols,"Pbar","TK"))) %>% distinct
                             
                                 # Extract melt proportion and add to the global table
                                 phm %>%
                                   filter(Name == self$meltName) %>%
                                   select(uniqueId,wt) %>%
                                   {.} -> melt_prop
                                 
                                 self$modPoints %>%
                                   left_join(melt_prop,by="uniqueId") %>%
                                   replace_na(list(wt=0)) %>%
                                   rename(FF=wt) %>%
                                   {.} -> self$modPoints
                                 
                                 # Count the model points, while we're at it...
                                 self$nPoints <- nrow(self$modPoints)
                                 
                                 # Add glob variables to it
                                 self$modPoints %<>% left_join(
                                   phm %>%  filter(Name=="system") %>% select(all_of(c("uniqueId",self$globVar))) ,
                                   by="uniqueId"
                                   )
                                 
                                 #.. and remove the globvars from phm (we leave P and T, for convenience)
                                 self$phm %<>% select(all_of(setdiff(names(self$phm),self$globVar)))
                                 
                                 # All phase names - a useful variable to have somewhere
                                 self$phNames <- setdiff(unique(phm$Name),"system")
                                 
                                 # Calculate the composition of the solid, cs
                                 solids <- setdiff(self$phNames,self$meltName)
                                 self$mergeMinerals(minerals = solids,outMinName = "solid", drop=F, updatePhList = F)
                                 
                                 cat("..done!\n")
                                 invisible(self)
                               },
                               
                               
                               #' @description Calculate Ti activity
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self with added aTiO2 column
                               aTiO2 = function() {
                                 
                                 self$modPoints %<>% mutate(aTiO2 = exp( - ( vGibbs(Rt,Pbar*1e-4,TK ) - muTiO2) / R / TK ) )
                                 
                                 #The above code is rather brittle, better check...
                                 if( any(self$phm$aTiO2 > 1)  ){
                                   cat("Something is wrong here, activities shouldn't be greater than 1 ! \nMost likely this has to do with HSC_conversion tag (should be OFF)\n")
                                 }
                               
                                 invisible(self)
                               },
                               
                               #' @description Calculate Si activity
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self with added aSiO2 column
                               aSiO2 = function() {
                                 
                                 self$modPoints %<>% mutate(aSiO2 = exp( - ( vGibbs(Qtz,Pbar*1e-4,TK ) - muSiO2) / R / TK ) )
                                 
                                 #The above code is rather brittle, better check...
                                 if( any(self$phm$aSiO2 > 1)  ){
                                   cat("Something is wrong here, activities shouldn't be greater than 1 ! \nMost likely this has to do with HSC_conversion tag (should be OFF)\n")
                                 }
                                 
                                 invisible(self)
                               },
                               
                               #' @description Calculate O fugacity
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self with added fO2 and Delta_FMQ columns
                               fO2 = function() {
                                 
                                 self$modPoints %<>% mutate(fO2 = exp((muO2 - GibbsO2(0.0001,TK) ) / R / TK ) ,
                                                      Delta_FMQ = log10(fO2) - FMQ(Pbar,TK) )
                                 
                                 #The above code is rather brittle, better check...
                                 if( any(self$phm$fO2 > 1)  ){
                                   cat("Something is wrong here, activities shouldn't be greater than 1 ! \nMost likely this has to do with HSC_conversion tag (should be OFF)\n")
                                 }
                                 
                                 invisible(self)
                                 },
                               
                               #' @description Calculate Ti in zircon
                               #' @details Requires activities of Si and Ti, of course
                               #' We calculate the Ti[ppm] in zircon (even where there is no Zrn present, filter as required !)
                               #' And the "apparent" temeprature making the common assumptions on aTiO2 and aSi2O
                               #' @param assumedATiO2, assumedATiO2 : assumed activities (commonly 0.7 and 1 resp.)
                               #' @return modified version of self with added TiZrn and appTTi column
                               TiinZr = function(assumedATiO2=0.7,assumedASiO2=1) {
                                 
                                 self$modPoints %<>% 
                                   mutate(TiZrn = 10^( 5.711 - 4800 / TK - log10(aSiO2) + log10(aTiO2) ),
                                          appTTi = - 4800 / (log10(TiZrn) + log10(assumedASiO2) - log10(assumedATiO2) - 5.711) )
                            
                                 invisible(self)
                               },
                               
                               #' @description Read a table containing WR compositions and select one to be the source for the model
                               #' @param  srcCompo (named numeric vector) A fall-back solution allowing the user to supply the properties of the source directly
                               #' (note that if all fails, one may still do self$c0 <- ...)
                               #' @param srcFile (char) the name of a tab-delimited file containing rock compositions
                               #' @param srcName (char) in this file, the name of the rock to use as the source
                               #' @details srcFile and srcName are supplied as a (not very robust) convenience to extract surce composition from a (sane) GCDkit-like file. 
                               #' It may not work. Don't fret and make your own vector, load it with self$getKd(srcCompo = ...)
                               #' @return modified version of self.
                               setSource = function(srcCompo=NULL,srcFile,srcName) {
                                 
                                 if(!is.null(srcCompo)){ # A fall-back solution, the responsibility is with the user
                                   self$c0 <- srcCompo
                                 }else{
                                   # We try to be clever and read from a file...
                                   c0 <- read.table(srcFile,sep="\t")
                                   c0 <- c0[srcName,]
                                   c0 <- c0[sapply(c0,is.numeric)]
                                   self$c0 <- as_vector(c0) # note that we use the tidyverse version as_vector,not the base as.vector
                                   }
                                 
                                # browser()
                                 
                                 invisible(self)
                               },
                               
                               #' @description Clean the source and the Kd file to keep only certain elements
                               #' @param intersect (boolean) keep only elements present BOTH in kd and C0, else rely on self$traces
                               #' @param clean.kd (boolean) trim the kd matrix
                               #' @param clean.c0 (boolean) trim the c0 vector
                               #' @details To keep a custom list of elements, write it in self$traces and use intersect = F
                               #' @return modified version of self.
                               cleanTrc = function(intersect=T,clean.kd=T,clean.c0=T) {
                                 
                                 if(intersect){
                                   trc <- intersect(colnames(self$kd),names(self$c0) )
                                 }else{
                                   trc <- self$traces
                                 }
                                 
                                 self$traces <- trc
                                 if(clean.kd){
                                   self$kd <- self$kd[,trc,drop=F]
                                 }
                                 if(clean.c0){
                                   self$c0 <- self$c0[trc]
                                 }
                                 
                                 invisible(self)
                               },
                               
                               
                               #' @description Read a table containing Kd values and, if required, check against existing phases in the model
                               #' @param kdFileName (char) Fully-qualified path to the kd file (relative to current directory)
                               #' @param fixnames (boolean) Fix mineral names as found in the kd files to new names, probably the perpleX compatible names of the model, 
                               #' or new names defined with self$renameMinerals() names into clean ones, 
                               #' using the dictionary supplied (see fields inMinNames, outMinNames)
                               #' @param inMinNames, outMinNames (char) List of names to translate; inMinNames are found in the Kd file, outMinNames are the new ones
                               #' @param check (boolean) If true (default), check that all the model names (i.e. self$phNames) 
                               #' have a line in the Kd file. If not, replace by the value of missingKd, probably 0 (default) or NA.
                               #' @param missingKd (num). What should we replace missing Kds with ?
                               #' @param extraMin. Extra minerals to add (or keep) to the Kd file, 
                               #' even if they are not present in the model. Typically used to preserve Kds for Zrn and Mnz.
                               #' @param keepOnly (boolean). If true, clean the kd table, keeping only the required minerals (i.e. c(self$phNames,extraMin) )
                               #' @return modified version of self adding the kd matrix.
                               readKd = function(kdFileName,fixnames=F,inMinNames=NULL,outMinNames=NULL,
                                                 check=T,missingKd=0,extraMin=c("Zrn","Mnz"),
                                                 keepOnly=T) {
                                 # cat("Reading file",kdFileName,"...")
                                 
                                 kd <- as.matrix(read.table(kdFileName) )
                                 
                                 if(fixnames){
                                   glossary <- data.frame(inMinNames,outMinNames)
                                   # A bit hackish because translateHeader works on colnames...
                                   kd <- t(translateHeader(t(kd),glossary,verbose=T))
                                 }
                                 
                                 properMins <- setdiff(c(self$phNames,extraMin),self$meltName )
                                 
                                 if(check){
                                   for(minName in properMins){
                                    if( !(minName %in% rownames(kd)) ){
                                      cn <-c(rownames(kd),minName)
                                      cat("WARNING:",minName,"is not found in Kd file. Replacing by",missingKd,"\n")
                                      kd %<>% rbind(missingKd)
                                      rownames(kd) <- cn
                                    }  
                                   }
                                  }
                                 
                                 if(keepOnly){
                                   kd <- kd[properMins,,drop=F]
                                 }
                                 
                                 self$kd <- kd
                                 
                                 invisible(self)
                               },
                               
                               #' @description Identify duplicate phases - i.e. those occuring more than once at a given point
                               #' @details  Typically, feldspar (ksp and plag). Other phases wth legitimate solvus may also be duplicate, but more often than not this is likely to be a perpleX issue - see solvus_tolerance, reach_increment, etc. in perplex_options.dat
                               #' @return modified version of self, updating self$phaseAssemblage, self$duplicates and self$dedup
                               countDuplicates = function(){
                               
                                 self$phaseAssemblage<- self$phm %>% 
                                   group_by(uniqueId) %>%  
                                   count(Name) %>% # Count the occurence of each phase
                                   spread(Name,n,fill=0) %>% # tabulate
                                   left_join(self$modPoints,by="uniqueId") %>% # Append P-T data
                                   ungroup() # Cleanup...
                              
                                 # Now we interpret this table to yield a list of duplicate phases
                                 self$duplicates <- self$phaseAssemblage %>% 
                                                  select(all_of(self$phNames)) %>%
                                                  select_if(~any(. > 1) ) %>% map(max) %>% as_tibble()
                                 
                                 cat("List of duplicate phases:\n")
                                 print(self$duplicates)
              
                                 self$dedup <- T
                                 invisible(self)
                               },
                               
                               #' @description Merge duplicate phases - i.e. those occuring more than once at a given point
                               #' @details The phases are simply merged, by a weighted average of their compositions. 
                               #' Beware - the reasons for having duplicates phases are generally something that you should into, either 
                               #' (i) legitimate solvi or (ii) poorly resolved phase boundaries in perpleX. 
                               #' Both should require your attention. Use as a last resort !  
                               #' 
                               #' If minerals contains more than one element, all of the minerals so named are aggregated.
                               #' @param minerals (char vector) The name(s) of mineral(s) to merge
                               #' @param outMinName New name to give to the "merged" mineral (unchanged, by default)
                               #' @param drop (boolean) drop or keep the original values ?
                               #' @param updatePhList (bbolean) Update the list of phases ?
                               #' @details both drop and updatePhList should be TRUE when merging proper minerals, 
                               #' and FALSE when calculating aggregate (e.g. solid) composition.
                               #' @return modified version of self, updating self$phm and self$phases
                               mergeMinerals = function(minerals,outMinName=minerals[1],drop=T,updatePhList=T){
                                 
                                 #if(!self$dedup){cat("WARNING: Run duplicate identification <model>$countDuplicates() first!\n")}
                                 cat("Merging and averaging all of",minerals,"as",outMinName,"...")
                              
                                 ## Make lists of columns to proceed
                                 
                                 # These colums are copied without alteration
                                 keep_col <- c("Name",self$pointCols,"from")
                                
                                 # For these cols (proportions), we add all the values found
                                 sum_col <- self$propCols
                                 
                                 # For the rest, we do a weigted mean
                                 mean_col <- setdiff(names(self$phm)[sapply(self$phm, is.numeric)],c(keep_col,sum_col))
                                 
                                 # Calculate the new minerals based on these specs
                                 self$phm %>% 
                                   filter(Name %in% minerals) %>%
                                   group_by(uniqueId) %>% 
                                   summarise(across( all_of(keep_col), max ),
                                             across( all_of(mean_col), weighted.mean, wt),
                                             across( all_of(sum_col), sum),
                                             .groups = "drop"
                                   ) %>%
                                   mutate(from="merged") %>%
                                   mutate(Name = outMinName) %>%
                                   {.} -> newmins
                                 
                                 # Update the main table
                                 # remove the old mins
                                 if(drop){
                                   self$phm %<>% 
                                     filter(!(Name %in% minerals))
                                 }
                                 # Add the new ones
                                 self$phm %>% 
                                   bind_rows(newmins) %>%
                                   {.} -> self$phm
                                 
                                 # Update a list of phases 
                                 if(updatePhList){
                                   self$phNames <- setdiff(self$phNames,minerals)
                                   self$phNames <- c(self$phNames,outMinName)
                                 }
                                 
                                 
                                 cat("..done\n")
                                 invisible(self)
                               },
                               
                               #' @description Merge duplicate phases - i.e. those occuring more than once at a given point
                               #' @details The phases are simply merged, by a weighted average of their compositions. Beware - the resons for having duplicates phases are generally something that you should into, either (i) legitimate solvi or (ii) poorly resolved phase boundaries in perpleX. Both should require your attention. Use as a last resort !  
                               #' @param mineral The mineral name to merge
                               #' @param outMinName New name to give to the "merged" mineral (unchanged, by default)
                               #' @return modified version of self, updating self$phm and self$phases
                               mergeMineralsOkz = function(mineral,outMinName=mineral){
                                 
                                 if(!self$dedup){cat("WARNING: Run duplicate identification <model>$countDuplicates() first!\n")}
                                 cat("Merging and averaging all instances of",mineral,"...")
                                 
                                 ## Make lists of columns to proceed
                                 
                                 # These colums are copied without alteration
                                 keep_col <- c("Name",self$pointCols,"from")
                                 
                                 # For these cols (proportions), we add all the values found
                                 sum_col <- self$propCols
                                 
                                 # For the rest, we do a weigted mean
                                 mean_col <- setdiff(names(self$phm)[sapply(self$phm, is.numeric)],c(keep_col,sum_col))
                                 
                                 # Calculate the new minerals based on these specs
                                 self$phm %>% 
                                   filter(Name == mineral) %>%
                                   group_by(uniqueId) %>% 
                                   summarise(across( all_of(keep_col), max ),
                                             across( all_of(mean_col), weighted.mean, wt),
                                             across( all_of(sum_col), sum),
                                             .groups = "drop"
                                   ) %>%
                                   mutate(from="merged") %>%
                                   mutate(Name = outMinName) %>%
                                   {.} -> newmins
                                 
                                 # Update the main table
                                 self$phm %>% 
                                   filter(Name != mineral) %>% # drop the old mineral
                                   bind_rows(newmins) %>%
                                   {.} -> self$phm
                                 
                                 # Update a list of phases 
                                 self$phNames <- setdiff(self$phNames,mineral)
                                 self$phNames <- c(self$phNames,outMinName)
                                 
                                 cat("..done\n")
                                 invisible(self)
                               },
                               
                               #' @description A generalized function to classify minerals on simple chemical criteria
                               #' @details A naive implementation, plain cutoff value on a single column. However, one may first calculate an extra col before classifying, if required.  
                               #' @details The default values are appropriate for feldspar (8.7 wt% K2O corresponds to Or50 on the Ab-Or join)
                               #' @param mineral Name of the mineral to classify
                               #' @param element on which the classification is based
                               #' @param cutoff The cut-off value to separate species
                               #' @param name_lower, name_upper the names to give to minerals resp. above and below the cutoff
                               #' @return modified version of self, updating self$phm by changing the mineral name to either name_lower or name_upper
                               classifyMineral = function(mineral="feldspar",element="K2O",cutoff=8.7,name_lower="plag",name_upper="ksp"){
                                
                                 cat("Classifying",mineral,"\n")
                                 upper_list <- self$phm[[element]] >  cutoff & 
                                               self$phm$Name == mineral &
                                               self$phm$from != "completed"
                                 lower_list <- self$phm[[element]] <= cutoff &
                                               self$phm$Name == mineral &
                                               self$phm$from != "completed"
                                 
                            
                                 self$phm[lower_list,"Name"] <- name_lower
                                 self$phm[upper_list,"Name"] <- name_upper
                                 
                                 self$phNames <- setdiff(self$phNames,mineral)
                                 if(any(upper_list)){
                                   self$phNames <- c(self$phNames,name_upper)
                                 }
                                 if(any(lower_list)){
                                   self$phNames <- c(self$phNames,name_lower)
                                 }
                                 
                                 invisible(self)
                               },
                               
                               
                               #' @description Rename one or several minerals throughout
                               #' @details As a by product, this also sorts minerals such that the first will appear at the botton of phase maps
                               #' @param minerals (initial) Name of the mineral to rename
                               #' @param outMinNames (new) Name of the mineral renamed
                               #' @param verbose More blah-blah
                               #' @return modified version of self, updating self$phm by changing the mineral names and updating self$phNames
                               renameMinerals = function(minerals,outMinNames=minerals,verbose=T){
                                 
                                 # As a by product, we reverse-sort the minerals (in phase plot) so by rev-ing the input we ensure that mins will plot with the first below
                                 minerals <- rev(minerals)
                                 #### FIXME ####
                                 outMinNames <- rev(outMinNames) # OUPS, will fail of the order is not the same in both lists.
                                 
                                 .renameOne<-function(mineral,outMinName){
                                   
                                   
                                   if(!(mineral %in% self$phNames) ){
                                     cat("No",mineral,"found... skipping\n")
                                   }else{
                                     if(verbose){
                                       cat("Renaming",mineral,"to",outMinName,"\n")}
                                       
                                       m_list <- self$phm$Name == mineral
                                       
                                       self$phm[m_list,self$nameCol] <- outMinName
                                       
                                       self$phNames <- setdiff(self$phNames,mineral)
                                       self$phNames <- c(self$phNames,outMinName)
                                   }
                                   
                                 }
                                 
                                 walk2(minerals,outMinNames,.renameOne)
                                 
                                 self$dedup<-F
                                 invisible(self)
                               },
                               
                               #' @description Completes (tidies) the phm table
                               #' @details This method ensures that the phm table contains one row per unique combination of phase and model point. Phases absent from a point get a 0 mode.
                               #' @return modified version of self, updating self$phm by adding required rows.
                               completeTable = function(){
                               
                                 cat("Be careful, don't try to rename/merge/identify minerals after that !\n")
                                 cat("Completing the table with all phases for all points:\n")
                                 # fillin <- c (self$nameCol,self$pointCols)
                                 
                                 # in principle, it should be possible to pass the names of the cols with
                                 # complete(all_of(self$pointCols))
                                 # however as of 11/2020 this does not work
                                 # https://github.com/tidyverse/tidyr/issues/1033
                                 # Workaround :
                                
                                 self$phm %>% complete(Name,
                                     nesting(uniqueId,!!!rlang::syms(self$pointCols)) ) %>%
                                     mutate( across( all_of(self$propCols) ,replace_na, 0 )) %>%
                                     mutate(from = replace_na(from, "completed")) %>%
                                     {.} -> self$phm
                                 
                                 self$complete <- T
                                 
                                 invisible(self)
                               },
                               
                               #' @description Generates more model point near the eutectic
                               #' @details Eutectics, by definition, occur at constant T but many things change nonetheless (melt amount in particular). For some applications it may be advantageous to have more point at the eutectic. We densify the data table accordingly, to have approx. one new point per wt% of melt abundance change. This is expected to work only for a 1D-model !
                               #' @return modified version of self, updating self$phm by adding more points at eutectic.
                               refineEutectic = function(FStep=1){
                                 
                                 cat("Attempting to fix eutectic:\n")
                                 if(length(self$pointCols)>1){cat("WARNING: are you trying to do eutectic on a 2D grid ? This would be a BAD idea ! \n")}
                                 if(!self$complete){cat("WARNING: you are trying to work on a table without phase proportions for all points.\nThis is likely to fail! Run <model>$complete() first \n")}
                               
                                 ## Find eutectic (specifically the colder point with non-0 melt amount)
                                 # Bracket it between "above" and below, collect identifiers
                                 
                                 ## The system above eutectic
                                 self$phm %>% 
                                   filter(Name == self$meltName) %>%
                                   filter(wt >0) %>%
                                   filter(.data$TK == min(.data$TK)) %>%
                                   pull(uniqueId) %>%
                                   {.} -> IdAbove
                                 
                                 phmAbove<- self$phm %>% filter(uniqueId == IdAbove)
                                 
                                 # Phase list
                                 phmAbove %>% 
                                   filter(wt > 0) %>%
                                   pull(Name) %>%
                                   {.} -> phAbove
                                 
                                 # Melt amount
                                 phmAbove %>% 
                                   filter(Name == self$meltName) %>%
                                   pull(wt) %>%
                                   {.} -> Fabove
                                 
                                 # Temperature
                                 phmAbove %>% pull(TK) %>%
                                   unique %>%
                                   {.} -> Tabove
                                 
                                 ## The system below eutectic
                                 self$phm %>% 
                                   filter(Name == self$meltName) %>%
                                   filter(wt <= 0) %>%
                                   filter(.data$TK == max(.data$TK)) %>%
                                   pull(uniqueId) %>%
                                   {.} -> IdBelow
                                 
                                 phmBelow<- self$phm %>% filter(uniqueId == IdBelow)
                                 
                                 # Phase list
                                 phmBelow %>% 
                                   filter(wt > 0) %>%
                                   pull(Name) %>%
                                   {.} -> phBelow
                                 
                                 # Temperature
                                 phmBelow %>% pull(TK) %>%
                                   unique %>%
                                   {.} -> Tbelow
                                 
                                 ## Changes at eutectic
                                 phOut <- setdiff(phAbove,phBelow)
                                 phIn <- setdiff(phBelow,phAbove)
                                 cat("Eutectic found between",Tabove,"and",Tbelow,";\n",
                                     "Phases out at eutectic:",phOut,"; phases in:",phIn,";\n",
                                     "melt crystallized:",Fabove,"%\n")
                                
                                 ## Melt amount(s) for which we wish to calculate
                                 XRange<-seq(0,1,length.out = max(3,ceiling(Fabove/FStep)) )
                                 XRange<-XRange[2:(length(XRange)-1)]
                                 
                                 cat(length(XRange),"steps will be added.\n")
                                 orL<-nrow(self$phm) # For reporting
                                 
                                ## Make lists of columns to proceed
                                phmEutectic <- bind_rows(phmAbove,phmBelow)
                                
                                totalPhases <- length(self$phNames)
                                
                                # The system inherits mostly its properties from the high-T side
                                # Except H rho and Cp, for which we average (a slight simplification, we are not sure that "below" is the actual foothill)
                                # And the variance, that we recalculate
                                ##### FIXME ##### 
                                # (hard coded var names: This fails if these variables have been renamed)
                                mean_col <- c(self$propCols,"H,J/mol","rho,kg/m3","cp,J/K/mol")
                                # variance_col<-"Counter"
                                
                                # These two will be set manually
                                special_col <- c("from","TK")
                                
                                # The rest is inherited
                                inherit_col <- setdiff(names(phmEutectic),c(mean_col,special_col)) #,variance_col
                                
                                for(X in XRange){
                                  
                                  # We make a new block (this) by merging the properties of each mineral (group_by) below and above solidus
                                  # Merge rules differ depending on column
                                  
                                  phmEutectic %>%
                                    add_column(X=c(rep(X,totalPhases+2),rep(1-X,totalPhases+2) )) %>%
                                    group_by(Name) %>%
                                    summarise( across( all_of(mean_col), weighted.mean, (1-X) ),.groups = "drop" ) %>% 
                                    left_join( phmAbove %>% select(all_of(inherit_col) ), by="Name") %>%
                                    mutate(from="eutectic",
                                           TK = Tabove - X/100) %>%
                                    unite("uniqueId",all_of(self$pointCols),remove=F) %>%
                                    {.} -> this
                                  
                                  # Add to data
                                  self$phm %<>% bind_rows(this)
                               
                                  # update the list of points (and their properties)
                                  newId = this$uniqueId[1]
                                  
                                  self$modPoints %>% filter(uniqueId==IdAbove) %>%
                                    mutate(TK = Tabove - X/100, uniqueId=newId, FF = Fabove * (1-X) ) %>%
                                    {.} -> newpoint
                                    
                                  # newpoint<-this[1,]%>%select(all_of(c(self$pointCols,"uniqueId")))
                                  self$modPoints %<>% add_row(newpoint)
                                  self$nPoints <- self$nPoints + 1
                                  
                                }
                                
                                 cat("phm table has grown by",nrow(self$phm)-orL,"lines\n")
                                 
                                 self$eutectic <- T
                                 invisible(self)
                               },
                               
                               #' @description Plots a phase proportion map as a function of temperature
                               #' @details  Assumes a 1D model with T as the main variable. Not generalized. In fact, not flexible at all - the colours are fixed, etc.
                               #' @param Tunit How should T be expressed ? "K" (default) or "C"
                               #' @return (invisibly) A ggplot graph
                               plotPhaseMap = function(Tunit="K",draw=T){
                                 # Assumes a 1D-file with T as the main variable
                                 
                                 if(Tunit=="C"){
                                   xlab <- "T (C)"
                                   Toffset <- -273
                                 }else{
                                   xlab <- "T (K)"
                                   Toffset <- 0
                                 }
                                
                                 TRange <- c(min(self$phm$TK),max(self$phm$TK)) + Toffset
                                 
                                 # Reorder phases to ensure we always have H2O and melt at the top
                                 special <- self$meltName
                                 mycols <- "palegoldenrod"
                                 if(any(self$phNames == "H2O")){
                                   special <- c("H2O",special)
                                   mycols <- c("black",mycols)
                                 }
                                 
                                 otherphases <- setdiff(self$phNames,special)
                                 
                                 phaseOrder <- c(special, otherphases)
                                 mycols<-c(mycols, brewer.pal(length(otherphases),"Paired") )
                                 #mycols<-c(mycols,hue_pal(h=c(0,360),l=80,c=80)(length(otherphases)))
                                
                                 p<- self$phm %>%
                                   filter(Name %in% self$phNames) %>%
                                   arrange(TK) %>%
                                   mutate(Phase = factor(Name, levels= phaseOrder)) %>%
                                   ggplot()+
                                   geom_area(mapping = aes(x=TK + Toffset,y=wt,fill=Phase ))+ 
                                   labs(y= "Phase proportion", x = xlab, title = self$modName )+
                                   scale_x_continuous(limits = TRange, expand = c(0, 0)) +
                                   scale_y_continuous(limits = c(0,100.01), expand = c(0, 0)) +
                                  # scale_fill_brewer(palette = "Paired")+
                                   scale_fill_manual(values = mycols )+
                                   theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=14,face="bold"))
                                 
                                 if(draw){
                                   print(p)  
                                 }
                                 
                                 invisible(p)
                               },
                               
                               #' @description Calculate (mostly) trace element partitioning during melting
                               #' @details The returned object contains only "partial" phm table, so one may need to run completeTable() after
                               #' @return modified version of self, updating self$phm and self$modPoint
                               #' @details the actual value depends on the melting function supplied
                               #' @param meltModel (fn) a melting function with formals function(self=self,melt_prop,phase_prop,thisPhm,thisGlob,...)
                               #' @param makeMinerals (char vector) A list of new minerals to create (probably because the melting model will form them). They will be empty shells in phm and present in phNames
                               #' @param makePropCols (char vector) New proportion columns (added to self$propCols, used for eutectic)
                               #' @param ... extra parameters passed to meltModel()
                                doMelt = function(meltModel,makeMinerals="",makePropCols=NULL,...){
                               
                                 .meltHere <- function(thisPoint){
                                   
                                   # Get the blocks of phases and system variables for this point
                                   self$phm %>% filter(uniqueId == thisPoint) %>%
                                   {.} -> thisPhm
                                   
                                   self$modPoints %>% filter(uniqueId == thisPoint) %>%
                                   {.} -> thisGlob
                                   
                                   # Convenience
                                   melt_prop <- thisGlob %>% pull(FF)
                                   # cat("F=",melt_prop,"..\n")
                                   
                                   thisPhm %>% filter( !(Name %in% c(self$meltName,"system","solid") ) & wt >0 ) %>%
                                     pull(wt,name=Name) %>%
                                     {.} ->  phase_prop
                                   
                                   # Do it
                                   model<-meltModel(self,melt_prop,phase_prop,thisPhm,thisGlob,...)
                                 
                                   thisPhm %>% inner_join(model$tr_point_data,by="Name") %>%
                                     {.} -> newPhm
                                   
                                   thisGlob %>% bind_cols(model$tr_glob_data) %>%
                                     {.} -> newGlob
                                   
                                   return(list(phData=newPhm,globData=newGlob))
                                 }
                                 ######### End inner function ########
                                 
                                 # Make empty new minerals as requested
                                 for(mx in makeMinerals){
                                   self$modPoints %>%
                                   select(all_of(c("uniqueId",self$pointCols))) %>%
                                     mutate(Name = mx,from="melting_routine") %>%
                                     {.} -> newMins 
                                   
                                   self$phm %<>% bind_rows(newMins)
                                   
                                   self$phNames <- c(self$phNames,mx)
                                 }
                                 
                                 self$propCols <- c(self$propCols,makePropCols)
                                   
                                   #browser()
                                   
                                   map(self$modPoints$uniqueId,.meltHere) %>%
                                    transpose %>%
                                    {.} -> ee
                                   
                                   ee$globData %>% tibble %>% unnest %>%
                                   {.} -> self$modPoints
                                     
                                   ee$phData %>% tibble %>% unnest %>%
                                   {.} -> self$phm
                                     
                                     
                                 # newphm <- NULL
                                 # newmod <- NULL
                                 #   
                                 # for(idx in self$modPoints$uniqueId){
                                 #   cat(idx,"..")
                                 # 
                                 #   np <- .meltHere(idx)
                                 #   newPoint <- np$phData
                                 #   newGlob <- np$globData
                                 # 
                                 #  # Update the model tibble
                                 #  newphm %>%
                                 #    bind_rows(newPoint) %>%
                                 #    {.} -> newphm
                                 #   
                                 #  # update the glob var tibble
                                 #  newmod %>%
                                 #    bind_rows(newGlob) %>%
                                 #    {.} -> newmod
                                 # 
                                 #  } # End of for loop
                                 # 
                                 # self$phm <- newphm
                                 # ## self$completeTable()
                                 # self$modPoints <- newmod
                                 
                                 invisible(self)
                               }
                               
                             )) # End of clas defn





