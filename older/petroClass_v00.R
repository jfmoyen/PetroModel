library(R6)
library(tidyverse)
library(scales)
library(RColorBrewer)


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
                        
                        #' @field Number of points in the model
                        nPoints= 0,
                        
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
                               
                               #' @field Header rows in the phm filme, above the colnames
                               hdrRows=8,
                               
                               #' @field Translation of phm to clean names: original (perpleX) names
                               ppxNames = c("H2O,wt%","MgO,wt%","Al2O3,wt%","SiO2,wt%","K2O,wt%","CaO,wt%","FeO,wt%","O2,wt%","TiO2,wt%","Na2O,wt%",
                                           "mu[TiO2],J/mol","mu[SiO2],J/mol","mu[O2],J/mol","P(bar)","T(K)","wt,%"),
                               
                               #' @field Translation of phm to clean names: new (clean) names
                               Names = c("H2O","MgO","Al2O3","SiO2","K2O","CaO","FeO","O2","TiO2","Na2O",
                                           "muTiO2","muSiO2","muO2","Pbar","TK","wt"),
                               
                               # For now, DO NOT RENAME "H,J/mol","rho,kg/m3","cp,J/K/mol" or bad things will happen (in fixEutectic)
                         
                               #' @field Colum storing temperature variable
                               TCol = "TK",
                               
                               #' @field In this column, is T stored as C or K?
                               TisK = T,
                               
                               #' @field Column storing phase name
                               nameCol = "Name",
            
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
                               propCol = c("vol,%","wt","mol,%"),
                               
                               #' @field Preferred unit for phase proportions
                               propUnit = "wt",
                               
                               #' @field Column identifying unique model points
                               pointCol = c("TK","Pbar"),
                               
                               #' @field All model points
                               modPoints = NULL,
                            
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
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self. It contains phm, a tibble that keeps all the info from the phm (werami output) file
                               readPhm = function(phmFileName,fixnames=T,verbose=T) {
                                 cat("Reading file",phmFileName,"...")
                                 # Read the file
                                 phm<-as.tibble(read.table(phmFileName,skip=self$hdrRows,header=T,stringsAsFactors=F,check.names=F))
                                 
                                 # The tidyverse version is subtly different from the base R version, and does not work as well on phm format
                                 #phm2<-read_table2(phmFileName,skip=self$hdrRows)
                                 
                                 # Convert ppx names to clean names
                                 if(fixnames){
                                   glossary<-data.frame(self$ppxNames,self$Names)
                                   phm <- translateHeader(phm,glossary,verbose=verbose)
                                   self$rename <- T
                                 }
                                 
                                 # Make a list of model points
                                 # Add a convenience column called uniqueId - a unique identifier for each point
                                 # Workaround, should be distinct(phm,all_of(self$pointCol))
                                 # But all_of() does not work as advertised
                                 pts<- distinct(phm,!!!rlang::syms(self$pointCol) ) 
                                 
                                 ids<- pts %>% unite("uniqueId",all_of(self$pointCol))
                                  
                                 pts %>% add_column(ids) -> self$modPoints
                                  
                                 self$nPoints <- nrow(self$modPoints)
                                 
                                 # Assign to the main variable
                                 uniqueId <- phm %>% unite("uniqueId",all_of(self$pointCol)) %>% select(uniqueId)
                                
                                 self$phm <- phm %>% add_column(uniqueId) %>% add_column(from="werami") 
                                 
                                 # All phase names - a useful variable to have somewhere
                                 self$phNames <- setdiff(unique(phm$Name),"system")
                                 
                                 cat("..done!\n")
                                 invisible(self)
                               },
                               
                               retrievePointXXX = function(pointLabel){
                                 ee <- self$phm %>% filter(uniqueId == pointLabel)
                                 return(ee)
                               },
                               
                               #' @description Retrieve the information corresponding to one model point
                               #' @param values a vector or one-line tibble or data.frame containing the coordinates of the point
                               #' @return A tibble with all the values for this point in the phm table
                               retrievePointOKZ = function(values){
                                 # Give each value its proper name
                                 # and force it to be a vector
                                 names(values) <- self$pointCol
                                 values <- unlist(values)
                                 
                                 ee<-self$phm
                                 
                                 # We subset sequentially the table for all the point-defining variables
                                 for(pdv in self$pointCol){
                                   
                                   ee %>% 
                                     filter(.data[[pdv]] == values[pdv]) -> ee
                                 }
                                 return(ee)
                               },
                               
                               #' @description For a model point, compute the phase assemblage table
                               #' @param values a vector or one-line tibble or data.frame containing the coordinates of the point
                               #' @return A tibble with all the values for this point in the phm table
                               listPhases = function(pointLabel){
                                 
                                 ee<- self$phm %>% filter(uniqueId == pointLabel) %>% 
                                   pull(self$nameCol) %>% 
                                   table() %>%
                                   bind_rows
                                 
                                 return(ee)
                               },
                               
                              
                               #' @description Identify duplicate phases - i.e. those occuring more than once at a given point
                               #' @details  Typically, feldspar (ksp and plag). Other phases wth legitimate solvus may also be duplicate, but more often than not this is likely to be a perpleX issue - see solvus_tolerance, reach_increment, etc. in perplex_options.dat
                               #' @return modified version of self, updating self$duplicates and self$dedup
                               countDuplicates = function(){
                               
                                 cat("Counting duplicates:... \n")
                                 # Map phase count to all model points
                                 self$phaseAssemblage<-  self$modPoints %>% pull(uniqueId) %>% 
                                   map(self$listPhases) %>% 
                                   tibble %>% unnest(cols=c(.)) %>%
                                   mutate_all(as.integer) %>% replace(is.na(.), 0) %>%
                                   bind_cols(self$modPoints)
                                 
                                 # Now we interpret this table to yield a list of duplicate phases
                                self$duplicates <- self$phaseAssemblage %>% 
                                                  select(all_of(self$phNames)) %>%
                                                  select_if(~any(. > 1) ) %>% map(max) %>% as_tibble()
                                 
                                 print(self$duplicates)
                               
                                 cat("..done !\n")
                                
                                 self$dedup <- T
                                 invisible(self)
                               },
                               
                               
                               mergeMinerals = function(mineral,outMinName=mineral){
                                 
                                 if(!self$dedup){cat("WARNING: Run duplicate identification <model>$countDuplicates() first!\n")}
                                 cat("Merging and averaging all instances of",mineral,"...")
                               
                                 # Get the average composition for a phase
                                 .average <- function(pointLabel){
                                   
                                   # Get the compositions of the mineral versions present here
                                   self$phm %>% filter(uniqueId == pointLabel) %>% 
                                                  filter(.data[[self$nameCol]] == mineral) -> to_process
                                   
                                   # Build an empty template
                                   newmin <- to_process[1,]
                                     
                                   # Everything is just the weighted average of the different versions...
                                   wta <- to_process %>%  summarise(across( where(is.numeric), weighted.mean, .data[[self$propUnit]]))
                                   
                                   # Except proportions, that are the sums
                                   tot <- to_process %>%  summarise(across( all_of(self$propCol), sum))
                                   
                                   # and point-defining variables that are the same !
                                   # No, averaging them is not good enough, rounding error can cause them to become very slightly different
                                   
                                   # Put together
                                   newmin[,colnames(wta)] <- wta[,colnames(wta)]
                                   newmin[,self$pointCol] <- to_process[1,self$pointCol] # P and T !
                                   newmin[,self$propCol]<-tot
                                   
                                   # Adjust names etc.
                                   newmin[,self$nameCol]<-outMinName
                                   newmin[,"from"] <- "merged"

                                   return(newmin)
                                 }
                                 
                                 
                                 # Identify the model points where two phases coexist
                                 treatMe <- self$phaseAssemblage[[mineral]] > 1
                                 treatPoints <- self$modPoints[treatMe,] %>% pull(uniqueId)

                                 # For each of those, calculate the new composition
                                 newmins <- treatPoints %>% map(.average) %>% bind_rows

                                 # Remove the old compos and add the new ones instead
                                 self$phm %>% 
                                   filter( !(uniqueId %in% treatPoints & .data[[self$nameCol]] == mineral ) ) %>%
                                   add_row(newmins) -> self$phm
                                 
                                 # Update also the other names
                                 self$phm[ self$phm[[self$nameCol]] == mineral, self$nameCol] <- outMinName
                                 
                                 self$phNames <- setdiff(self$phNames,mineral)
                                 self$phNames <- c(self$phNames,outMinName)
                                 
                                 cat("..done\n")
                                 invisible(self)
                               },
                               
                               #' @description A generalized function to classify minerals on simple chemical criteria
                               #' @details A naive implementation, plain cutoff value on a single column. However, one may first calculate an extra col before classifying, if required.  
                               #' @param mineral Name of the mineral to classify
                               #' @param element on which the classification is based
                               #' @param cutoff The cut-off value to separate species
                               #' @param name_lower, name_upper the names to give to minerals resp. above and below the cutoff
                               #' @return modified version of self, updating self$phm by changing the mineral name to either name_lower or name_upper
                               classifyMineral = function(mineral="feldspar",element="K2O",cutoff=8.7,name_lower="plag",name_upper="ksp"){
                                
                                 cat("Classifying",mineral,"\n")
                                 upper_list <- self$phm[[element]] >  cutoff & 
                                               self$phm[[self$nameCol]] == mineral &
                                               self$phm$from != "completed"
                                 lower_list <- self$phm[[element]] <= cutoff &
                                               self$phm[[self$nameCol]] == mineral &
                                               self$phm$from != "completed"
                                 
                            
                                 self$phm[lower_list,self$nameCol] <- name_lower
                                 self$phm[upper_list,self$nameCol] <- name_upper
                                 
                                 self$phNames <- setdiff(self$phNames,mineral)
                                 self$phNames <- c(self$phNames,name_lower,name_upper)
                                 
                                 invisible(self)
                               },
                               
                               #' @description Attempts to attribute the feldspar or feldspars to either plag or ksp
                               #' @details A naive implementation, everything above Or50 (= 8.7 % K2O, customizable) is a ksp  
                               #' @param cutoff The cut-off value to separate Ksp from plag. Default 50% on An-Or join 
                               #' @return modified version of self, updating self$phm by changing the feldspar name to eithjer "plag" or "ksp"
                               fixFeldspar = function(cutoff=8.7){
                         ### DEPRECATED
                                     kspList <-self$phm$K2O >  cutoff & self$phm[[self$nameCol]] == "feldspar"
                                     plaglist<-self$phm$K2O <= cutoff & self$phm[[self$nameCol]] == "feldspar"
                                     
                                     self$phm[plaglist,self$nameCol] <- "plag"
                                     self$phm[kspList,self$nameCol] <- "ksp"
                                     
                                     self$phNames <- setdiff(self$phNames,"feldspar")
                                     self$phNames <- c(self$phNames,"plag","ksp")
                                    
                                 invisible(self)
                               },
                               
                               #' @description Attempts to attribute the feldspar or feldspars to either plag or ksp
                               #' @details  This method requires that duplicates where idedentified first (see self$countDuplicates). The classification is purely based on K2O, with a 2 % threshold. 
                               #' @return modified version of self, updating self$phm by changing the feldspar name to eithjer "plag" or "ksp"
                               fixFeldsparOKZ = function(){
                          ### VERY DEPRECATED
                                 if (!self$dedup){stop("Please identify duplicates first")}
                                 
                                 if(  any(names(self$duplicates) == self$feldsparName) ){
                                   cat("Two feldspars found, classifying on the basis of K2O contents\n")
                                   
                                   ### new
                                   # .felsClass<-function(mp){
                                   #   self$retrievePoint(mp) %>% 
                                   #     pull(self$nameCol) %>% 
                                   #     table() %>% 
                                   #     return()
                                   # }
                                   # 
                                   # self$modPoints %>% 
                                   #   transpose %>% map(.felsClass) %>% 
                                   ### end new
                                   
                                   ### K2O hard-coded !
                                   kspList <-self$phm$K2O >  2 & self$phm[[self$nameCol]] == "feldspar"
                                   plaglist<-self$phm$K2O <= 2 & self$phm[[self$nameCol]] == "feldspar"
                                   
                                   self$phm[plaglist,self$nameCol] <- "plag"
                                   self$phm[kspList,self$nameCol] <- "ksp"
                                   
                                   self$phNames <- setdiff(self$phNames,"feldspar")
                                   self$phNames <- c(self$phNames,"plag","ksp")
                                   
                                 }else{
                                   cat("Only one feldspar - classyfying:\n")
                                   hasKsp<-F
                                   hasPlg<-F
                                   
                                   if(any(self$phm$K2O >  2 & self$phm[[self$nameCol]] == "feldspar")){hasKsp<-T}
                                   if(any(self$phm$K2O <  2 & self$phm[[self$nameCol]] == "feldspar")){hasPlg<-T}
                                   
                                   if(hasKsp & hasPlg){
                                     cat("...something is weird, composition varies a lot..  I'd rather not do it automatically. Sorry !\n")
                                   }else{
                                     self$phNames <- setdiff(self$phNames,"feldspar")
                                     if(hasKsp){
                                       self$phm[self$phm[[self$nameCol]] == self$feldsparName, self$nameCol] <- "ksp"
                                       cat("Feldspar is classified as ksp\n")
                                       self$phNames <- c(self$phNames,"ksp")
                                     }
                                     if(hasPlg){
                                       self$phm[self$phm[[self$nameCol]] == self$feldsparName, self$nameCol] <- "plag"
                                       cat("Feldspar is classified as plag\n")
                                       self$phNames <- c(self$phNames,"plag")
                                     }
                                     
                                   }
                                 }
                                 invisible(self)
                               },
                               
                               #' @description Completes (tidies) the phm table
                               #' @details This method ensures that the phm table contains one row per unique combination of phase and model point. Phases absent from a point get a 0 mode.
                               #' @return modified version of self, updating self$phm by adding required rows.
                               completeTable = function(){
                               
                                 cat("Completing the table with all phases for all points:\n")
                                 fillin <- c (self$nameCol,self$pointCol)
                                 # in principle, it should be possible to pass the names of the cols with
                                 # complete(phm, all_of(fillin))
                                 # however as of 11/2020 this does not work
                                 # https://github.com/tidyverse/tidyr/issues/1033
                                 # Workaround :
                                
                                 self$complete <- T
                                
                                 self$phm <- self$phm %>% 
                                   complete( !!!rlang::syms(fillin)) %>% # should read: complete(all_of(fillin))
                                   mutate( across( all_of(self$propCol) ,replace_na, 0 )) %>%
                                   mutate(from = replace_na(from, "completed")) %>%
                                   select(!uniqueId)
                                   
                                 # Regenerate unique Id (Ok, wasteful, we could just add it to the new lines...)
                                 uniqueId <- self$phm %>% unite("uniqueId",all_of(self$pointCol)) %>% select(uniqueId) 
                                 
                                 self$phm <- self$phm %>% add_column(uniqueId)
                                 
                                 cat("..done !\n")
                                 
                                 invisible(self)
                               },
                               
                               #' @description Generates more model point near the eutectic
                               #' @details Eutectics, by definition, occur at constant T but many things change nonetheless (melt amount in particular). For some applications it may be advantageous to have more point at the eutectic. We densify the data table accordingly, to have approx. one new point per wt% of melt abundance change. This is expected to work only for a 1D-model !
                               #' @return modified version of self, updating self$phm by adding more points at eutectic.
                               fixEutectic = function(FStep=1){
                                 
                                 cat("Attempting to fix eutectic:\n")
                                 if(length(self$pointCol)>1){cat("WARNING: are you trying to do eutectic on a 2D grid ? This would be a BAD idea ! \n")}
                                 if(!self$complete){cat("WARNING: you are trying to work on a table without phase proportions for all points.\nThis is likely to fail! Run <model>$complete() first \n")}
                               
                                 # The table must be ordered by increasing T !
                                 self$phm <- arrange(self$phm,!!!rlang::syms(self$TCol))
                                 
                                 # Find eutectic (specifically the melt-line that contains the less non-0 melt amount)
                                 
                                 self$phm %>% 
                                   filter(.data[[self$nameCol]] == self$meltName) %>%
                                   filter(.data[[self$propUnit]] >0) %>%
                                   filter(.data[[self$propUnit]] == min(.data[[self$propUnit]])) ->
                                   lastDrop
                                 
                                 # Where are we, in terms of F and T
                              
                                 Fabove <- lastDrop %>% pull(self$propUnit)
                                 Tabove <- lastDrop %>% pull(self$TCol)
                                 ptabove<- filter(self$modPoints, .data[[self$TCol]] == Tabove)
                                 
                                 Tbelow <- self$modPoints %>% 
                                            filter(.data[[self$TCol]] < Tabove) %>%
                                            pull(TK) %>%
                                            max
                                 ptbelow<- filter(self$modPoints, .data[[self$TCol]] == Tbelow)
                                 
                                 cat("Eutectic between",Tabove,"and",Tbelow,"\n")
                                 # The system sub and super solidus
                                 # This assumes that there is only one point per temperature (ie, no 2D grids)
                                 phmAbove<- self$phm %>% filter(uniqueId == ptabove$uniqueId)
                                 phmBelow<- self$phm %>% filter(uniqueId == ptbelow$uniqueId)
                                 
                                 # phmAbove<-self$retrievePoint(ptabove)
                                 # phmBelow<-self$retrievePoint(ptbelow)
                               
                                 # Describe the eutectic
                                 # This fails if the data has been completed already !
                                 phAbove <- phmAbove %>% 
                                   filter(.data[[self$propUnit]] >0) %>%
                                   pull(self$nameCol)

                                 phBelow <- phmBelow %>% 
                                   filter(.data[[self$propUnit]] >0) %>%
                                   pull(self$nameCol)
                                 
                                 phOut <- setdiff(phAbove,phBelow)
                                 phIn <- setdiff(phBelow,phAbove)
                                 
                                 cat("Phases out at eutectic:",phOut,"; phases in:",phIn,"; melt crystallized:",Fabove,"%\n")
                                
                                  # Melt amount for which we wish to calculate
                                 
                                 XRange<-seq(0,1,length.out = Fabove/FStep)
                                 XRange<-XRange[2:(length(XRange)-1)]
                                 
                                 cat(length(XRange),"steps will be added.\n")
                                 
                                 # Clone the super-solidus block, to have a tibble of the right shape
                                 # Also, most of the phase properties don't change here, so it's convenient to start from that
                                 template <- phmAbove
                                
                                 # if a phase appears on the eutectic, it was not there above ! Adjust the template
                                 for(p in phIn){
                                   template <- template %>% add_row(phmBelow[phmBelow[[self$nameCol]]==p,])
                                    }
                                 
                                 # The following data is generated from eutectic interpolation
                                 template$from<-"eutectic"
                           
                                 # Now to work... 
                                 orL<-nrow(self$phm)
                                 
                                 for(X in XRange){
                                   # Generate a new data block
                                   block<-template
                                   
                                   # The system inherits mostly its properties from the high-T side
                                   # Except H rho and Cp, for which we average (a slight simplification, we are not sure that "below" is the actual foothill)
                                   # This fails if these variables have been renamed
                                   ##### FIXME ##### 
                                   # (hard coded var names)
                                  
                                   for(v in c("H,J/mol","rho,kg/m3","cp,J/K/mol")){ 
                                     block[block[[self$nameCol]]=="system",v]<- (1-X)* phmAbove[phmAbove[[self$nameCol]]=="system",v] + X * phmBelow[phmBelow[[self$nameCol]]=="system",v]
                                   }
                                   
                                   # ... NOT the T dependent values, such as mu. these should stay at the above T.
                                   
                                   # For most of the stuff (e.g., all phases, exc. system) nothing changes
                                   # just the props need to be adjusted !
                                   
                                   for(ph in self$phNames){
                                     for(v in self$propCol){
                                     
                                      block[block[[self$nameCol]]==ph,v]<-
                                        (1-X)* phmAbove[phmAbove[[self$nameCol]]==ph,v] +
                                         X * phmBelow[phmBelow[[self$nameCol]]==ph,v]
                                                           }
                                   }
                                   
                                   # Adjust the variance
                                   ##### FIXME #####
                                   ## Hard-coded
                                   block$Counter[!is.na(block$Counter)]<-length(block$Counter[!is.na(block$Counter)]) -1
                                   
                                   # Temperature has been overwritten by the previous calculations - fix it
                                   # We assume the "true" eutectic temp is the "above" temp, to ensure consistent mu etc.
                                   # Also to avoid further trouble we jitter very slightly the eutectic temperatures so that they are all different
                                   
                                   block[[self$TCol]] <- Tabove - X / 100
                                   
                                   # Create a dummy identifier
                                
                                   id<- block %>% unite("uniqueId",all_of(self$pointCol)) %>% pull(uniqueId)
                                   block$uniqueId <- id
                                   
                                   # And add to the data table...
                                   self$phm <- add_row(self$phm,block)
                                   
                                   # Add the same identifier to the list of model points
                                   newpoint<-block[1,]%>%select(all_of(c(self$pointCol,"uniqueId")))
                                   self$modPoints <- self$modPoints %>% add_row(newpoint)
                                  }
                                 
                                 cat("phm table has grown by",nrow(self$phm)-orL,"lines\n")
                                 
                                 self$eutectic <- T
                                 invisible(self)
                               },
                               
                               #' @description Plots a phase proportion map as a function of temperature
                               #' @details  Assumes a 1D model with T as the main variable. Not generalized. In fact, not flexible at all - the colours are fixed, etc.
                               #' @param Tunit How should T be expressed ? "K", "C" or "as.is" -- the default, i.e. relying on self$TisK.
                               #' @return (invisibly) A ggplot graph
                               plotPhaseMap = function(Tunit="as.is"){
                                 # Assumes a 1D-file with T as the main variable
                                 
                                 if(Tunit == "as.is"){
                                  graphTisK <- self$TisK
                                 }else{
                                   if(Tunit == "K"){
                                     graphTisK <- T
                                   }else{
                                     graphTisK <- F
                                   }
                                 }
                                 
                                 if(graphTisK){
                                   xlab <- "T (K)"
                                   if(self$TisK){
                                     Toffset <- 0
                                   }else{
                                     Toffset <- 273
                                   }
                                   
                                 }else{
                                   xlab <- "T (C)"
                                   if(self$TisK){
                                     Toffset <- -273
                                   }else{
                                     Toffset <- 0
                                   }
                                 }
                                 
                                 TRange <- c(min(self$phm[[self$TCol]]),max(self$phm[[self$TCol]])) + Toffset
                                 
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
                                   filter(.data[[self$nameCol]] %in% self$phNames) %>%
                                   arrange(.data[[self$TCol]]) %>%
                                   mutate(Phase = factor(.data[[self$nameCol]], levels= phaseOrder)) %>%
                                   ggplot()+
                                   geom_area(mapping = aes(x=.data[[self$TCol]] + Toffset,y=.data[[self$propUnit]],fill=Phase ))+ 
                                   labs(y= "Phase proportion", x = xlab, title = self$modName )+
                                   scale_x_continuous(limits = TRange, expand = c(0, 0)) +
                                   scale_y_continuous(limits = c(0,100.01), expand = c(0, 0)) +
                                  # scale_fill_brewer(palette = "Paired")+
                                   scale_fill_manual(values = mycols )+
                                   theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=14,face="bold"))
                                 print(p)
                                 
                                 invisible(p)
                               }
                               
                             ))


## KCG granite 4 test
KCGgrn<-PerplexModel$new("KCG Granite")
KCGgrn$readPhm("./GrTypes/KCGgrn_1.phm",verbose = F)
 KCGgrn$countDuplicates()
 KCGgrn$mergeMinerals(mineral="Ilm(WPH)",outMinName="ilmenite")
 #KCGgrn$classifyMineral(mineral="Ilm(WPH)",element="TiO2",cutoff=30,name_lower="ilm1",name_upper="ilm2")
 KCGgrn$classifyMineral(cutoff=6.6)
 KCGgrn$plotPhaseMap(Tunit="K")
 

# 
# 
KCGgrn$completeTable()
KCGgrn$fixEutectic()








