library(R6)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(magrittr) # Needs to be imported explicitly, because https://github.com/tidyverse/magrittr/issues/194 (WtF ?)

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
                               
                               #' @field Header rows in the phm file, above the colnames
                               hdrRows=8,
                               
                               #' @field Translation of phm to clean names: original (perpleX) names
                               ppxNames = c("H2O,wt%","MgO,wt%","Al2O3,wt%","SiO2,wt%","K2O,wt%","CaO,wt%","FeO,wt%","O2,wt%","TiO2,wt%","Na2O,wt%",
                                           "mu[TiO2],J/mol","mu[SiO2],J/mol","mu[O2],J/mol","P(bar)","T(K)","wt,%"),
                               
                               #' @field Translation of phm to clean names: new (clean) names
                               Names = c("H2O","MgO","Al2O3","SiO2","K2O","CaO","FeO","O2","TiO2","Na2O",
                                           "muTiO2","muSiO2","muO2","Pbar","TK","wt"),
                               
                               # For now, DO NOT RENAME "H,J/mol","rho,kg/m3","cp,J/K/mol" or bad things will happen (in refineEutectic)
                         
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
                               pointCol = c("TK","Pbar"), ## Counter ?
                               
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
                                 
                                 # We generate a UniqueId - an identifier for each model point
                                 phm %<>% unite("uniqueId",all_of(self$pointCol),remove=F)
                                 
                                 # And we flag the lines as coming from Werami
                                 phm %<>% add_column(from="werami") 
                                 
                                 # Assign to the real thing
                                 self$phm <- phm
                                 
                                 # Build a table of model points
                                 self$modPoints <- phm %>% select(all_of(c("uniqueId",self$pointCol))) %>% distinct
                                   
                                 # Count the model points, while we're at it...
                                 self$nPoints <- nrow(self$modPoints)
                                 
                                 # All phase names - a useful variable to have somewhere
                                 self$phNames <- setdiff(unique(phm$Name),"system")
                                 
                                 cat("..done!\n")
                                 invisible(self)
                               },
                               
                               #' @description Read a table containing Kd values and, if required, check against existing phases in the model
                               #' @param kdFileName (char) Fully-qualified path to the kd file (relative to current directory)
                               #' @param fixnames (boolean) Attempt to fix perpleX names into clean ones, using the dictionary supplied (see fields ppxNames, Names)
                               #' @param verbose (boolean) Blah, blah, blah... 
                               #' @details Should be able to work on any kind of perpleX phm output, either 1D or 2D.
                               #' @return modified version of self. It contains phm, a tibble that keeps all the info from the phm (werami output) file
                               readKd = function(kdFileName) {
                                 cat("Reading file",kdFileName,"...")
                                 
                                 kd <- read_table(kdFileName,na=0) 
                                 
                                 return(self)
                               },
                               
                               
                               
                               #' @description Identify duplicate phases - i.e. those occuring more than once at a given point
                               #' @details  Typically, feldspar (ksp and plag). Other phases wth legitimate solvus may also be duplicate, but more often than not this is likely to be a perpleX issue - see solvus_tolerance, reach_increment, etc. in perplex_options.dat
                               #' @return modified version of self, updating self$phaseAssemblage, self$duplicates and self$dedup
                               countDuplicates = function(){
                               
                                 self$phaseAssemblage<- self$phm %>% 
                                   group_by(uniqueId) %>%  
                                   count(get(self$nameCol)) %>% # Count the occurence of each phase
                                   spread(`get(self$nameCol)`,n,fill=0) %>% # tabulate
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
                               #' @details The phases are simply merged, by a weighted average of their compositions. Beware - the resons for having duplicates phases are generally something that you should into, either (i) legitimate solvi or (ii) poorly resolved phase boundaries in perpleX. Both should require your attention. Use as a last resort !  
                               #' @param mineral The mineral name to merge
                               #' @param outMinName New name to give to the "merged" mineral (unchanged, by default)
                               #' @return modified version of self, updating self$phm and self$phases
                               mergeMinerals = function(mineral,outMinName=mineral){
                                 
                                 if(!self$dedup){cat("WARNING: Run duplicate identification <model>$countDuplicates() first!\n")}
                                 cat("Merging and averaging all instances of",mineral,"...")
                               
                                 ######### Experimental
                                 # browser()
                                 
                                 ## Make lists of columns to proceed
                                 
                                 # These colums are copied without alteration
                                 keep_col <- c(self$nameCol,self$pointCol,"from")
                                
                                 # For these cols (proportions), we add all the values found
                                 sum_col <- self$propCol
                                 
                                 # For the rest, we do a weigted mean
                                 mean_col <- setdiff(names(phm)[sapply(phm, is.numeric)],c(keep_col,sum_col))
                                 
                                 # Calculate the new minerals based on these specs
                                 self$phm %>% 
                                   filter(get(self$nameCol) == mineral) %>%
                                   group_by(uniqueId) %>% 
                                   summarise(across( all_of(keep_col), max ),
                                             across( all_of(mean_col), weighted.mean, !!!rlang::syms(self$propUnit) ),
                                             across( all_of(sum_col), sum),
                                             .groups = "drop"
                                   ) %>%
                                   mutate(from="merged") %>%
                                   mutate("{self$nameCol}" := outMinName) %>%
                                   {.} -> newmins
                                 
                                 # Update the main table
                                 self$phm %>% 
                                   filter(get(self$nameCol) != mineral) %>% # drop the old mineral
                                   bind_rows(newmins) %>%
                                   {.} -> self$phm
                                 
                                 
                                 
                                 
                                 # self$phm %>% group_by(uniqueId) %>% 
                                 #   summarise(across( all_of(self$propCol), sum),
                                 #             across( where(is.numeric), weighted.mean, .data[[self$propUnit]])
                                 #             )
                                 # 
                                 # 
                                 # self$phm %>% group_by(uniqueId) %>% summarise(across( where(is.numeric), weighted.mean, .data[[self$propUnit]]))
                                 # 
                                 # self$phm %>% group_by(uniqueId) %>% 
                                 #   summarise(across( all_of(c(2,4)), sum),
                                 #             across( all_of(c(4,7)), mean)
                                 #   )
                                 
                                 
                                 ############
                                 # 
                                 # # Get the average composition for a phase
                                 # .average <- function(pointLabel){
                                 #   
                                 #   # Get the compositions of the mineral versions present here
                                 #   self$phm %>% filter(uniqueId == pointLabel) %>% 
                                 #                  filter(.data[[self$nameCol]] == mineral) -> to_process
                                 #   
                                 #   # Build an empty template
                                 #   newmin <- to_process[1,]
                                 #     
                                 #   # Everything is just the weighted average of the different versions...
                                 #   wta <- to_process %>%  summarise(across( where(is.numeric), weighted.mean, .data[[self$propUnit]]))
                                 #   
                                 #   # Except proportions, that are the sums
                                 #   tot <- to_process %>%  summarise(across( all_of(self$propCol), sum))
                                 #   
                                 #   # and point-defining variables that are the same !
                                 #   # No, averaging them is not good enough, rounding error can cause them to become very slightly different
                                 #   
                                 #   # Put together
                                 #   newmin[,colnames(wta)] <- wta[,colnames(wta)]
                                 #   newmin[,self$pointCol] <- to_process[1,self$pointCol] # P and T !
                                 #   newmin[,self$propCol]<-tot
                                 #   
                                 #   # Adjust names etc.
                                 #   newmin[,self$nameCol]<-outMinName
                                 #   newmin[,"from"] <- "merged"
                                 # 
                                 #   return(newmin)
                                 # }
                                 # 
                                 # 
                                 # # Identify the model points where two phases coexist
                                 # treatMe <- self$phaseAssemblage[[mineral]] > 1
                                 # treatPoints <- self$modPoints[treatMe,] %>% pull(uniqueId)
                                 # 
                                 # # For each of those, calculate the new composition
                                 # newmins <- treatPoints %>% map(.average) %>% bind_rows
                                 # 
                                 # # Remove the old compos and add the new ones instead
                                 # self$phm %>% 
                                 #   filter( !(uniqueId %in% treatPoints & .data[[self$nameCol]] == mineral ) ) %>%
                                 #   add_row(newmins) -> self$phm
                                 # 
                                 # # Update also the other names
                                 # self$phm[ self$phm[[self$nameCol]] == mineral, self$nameCol] <- outMinName
                                 # 
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
                                               self$phm[[self$nameCol]] == mineral &
                                               self$phm$from != "completed"
                                 lower_list <- self$phm[[element]] <= cutoff &
                                               self$phm[[self$nameCol]] == mineral &
                                               self$phm$from != "completed"
                                 
                            
                                 self$phm[lower_list,self$nameCol] <- name_lower
                                 self$phm[upper_list,self$nameCol] <- name_upper
                                 
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
                                       
                                       m_list <- self$phm[[self$nameCol]] == mineral
                                       
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
                                 # fillin <- c (self$nameCol,self$pointCol)
                                 
                                 # in principle, it should be possible to pass the names of the cols with
                                 # complete(all_of(self$pointCol))
                                 # however as of 11/2020 this does not work
                                 # https://github.com/tidyverse/tidyr/issues/1033
                                 # Workaround :
                                
                                 self$phm %>% complete(!!!rlang::syms(self$nameCol),
                                     nesting(uniqueId,!!!rlang::syms(self$pointCol)) ) %>%
                                     mutate( across( all_of(self$propCol) ,replace_na, 0 )) %>%
                                     mutate(from = replace_na(from, "completed")) %>%
                                     {.} -> self$phm
                                 
                                 self$complete <- T
                                 
                                # browser()
                                #  self$phm <- self$phm %>% 
                                #    complete( !!!rlang::syms(fillin)) %>% # should read: complete(all_of(fillin))
                                #    mutate( across( all_of(self$propCol) ,replace_na, 0 )) %>%
                                #    mutate(from = replace_na(from, "completed")) %>%
                                #    select(!uniqueId)
                                #    
                                #  # Regenerate unique Id (Ok, wasteful, we could just add it to the new lines...)
                                #  uniqueId <- self$phm %>% unite("uniqueId",all_of(self$pointCol)) %>% select(uniqueId) 
                                #  
                                #  self$phm <- self$phm %>% add_column(uniqueId)
                                #  
                                # cat("..done !\n")
                                 
                                 invisible(self)
                               },
                               
                               #' @description Generates more model point near the eutectic
                               #' @details Eutectics, by definition, occur at constant T but many things change nonetheless (melt amount in particular). For some applications it may be advantageous to have more point at the eutectic. We densify the data table accordingly, to have approx. one new point per wt% of melt abundance change. This is expected to work only for a 1D-model !
                               #' @return modified version of self, updating self$phm by adding more points at eutectic.
                               refineEutectic = function(FStep=1){
                                 
                                 cat("Attempting to fix eutectic:\n")
                                 if(length(self$pointCol)>1){cat("WARNING: are you trying to do eutectic on a 2D grid ? This would be a BAD idea ! \n")}
                                 if(!self$complete){cat("WARNING: you are trying to work on a table without phase proportions for all points.\nThis is likely to fail! Run <model>$complete() first \n")}
                               
                                 # The table must be ordered by increasing T ! NOT ?
                                 # self$phm <- arrange(self$phm,!!!rlang::syms(self$TCol))
                                 
                                 ## Find eutectic (specifically the colder point with non-0 melt amount)
                                 # Bracket it between "above" and below, collect identifiers
                                 
                                 ## The system above eutectic
                                 self$phm %>% 
                                   filter(.data[[self$nameCol]] == self$meltName) %>%
                                   filter(.data[[self$propUnit]] >0) %>%
                                   filter(.data[[self$TCol]] == min(.data[[self$TCol]])) %>%
                                   pull(uniqueId) %>%
                                   {.} -> IdAbove
                                 
                                 phmAbove<- self$phm %>% filter(uniqueId == IdAbove)
                                 
                                 # Phase list
                                 phmAbove %>% 
                                   filter(.data[[self$propUnit]] > 0) %>%
                                   pull(self$nameCol) %>%
                                   {.} -> phAbove
                                 
                                 # Melt amount
                                 phmAbove %>% 
                                   filter(.data[[self$nameCol]] == self$meltName) %>%
                                   pull(self$propUnit) %>%
                                   {.} -> Fabove
                                 
                                 # Temperature
                                 phmAbove %>% pull(self$TCol) %>%
                                   unique %>%
                                   {.} -> Tabove
                                 
                                 ## The system below eutectic
                                 self$phm %>% 
                                   filter(.data[[self$nameCol]] == self$meltName) %>%
                                   filter(.data[[self$propUnit]] <= 0) %>%
                                   filter(.data[[self$TCol]] == max(.data[[self$TCol]])) %>%
                                   pull(uniqueId) %>%
                                   {.} -> IdBelow
                                 
                                 phmBelow<- self$phm %>% filter(uniqueId == IdBelow)
                                 
                                 # Phase list
                                 phmBelow %>% 
                                   filter(.data[[self$propUnit]] > 0) %>%
                                   pull(self$nameCol) %>%
                                   {.} -> phBelow
                                 
                                 # Temperature
                                 phmBelow %>% pull(self$TCol) %>%
                                   unique %>%
                                   {.} -> Tbelow
                                 
                                 ## Changes at eutectic
                                 phOut <- setdiff(phAbove,phBelow)
                                 phIn <- setdiff(phBelow,phAbove)
                                 cat("Phases out at eutectic:",phOut,"; phases in:",phIn,"; melt crystallized:",Fabove,"%\n")
                                
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
                                mean_col <- c(self$propCol,"H,J/mol","rho,kg/m3","cp,J/K/mol")
                                variance_col<-"Counter"
                                
                                # These two will be set manually
                                special_col <- c("from",self$TCol)
                                
                                # The rest is inherited
                                inherit_col <- setdiff(names(phmEutectic),c(mean_col,special_col,variance_col))
                                
                                for(X in XRange){
                                  
                                  # We make a new block (this) by merging the properties of each mineral (group_by) below and above solidus
                                  # Merge rules differ depending on column
                                  
                                  phmEutectic %>%
                                    add_column(X=c(rep(X,totalPhases+1),rep(1-X,totalPhases+1) )) %>%
                                    group_by(!!!rlang:::syms(self$nameCol) ) %>%
                                    summarise( across( all_of(mean_col), weighted.mean, X ),.groups = "drop" ) %>% 
                                    left_join( phmAbove %>% select(all_of(inherit_col) ), by=self$nameCol  ) %>%
                                    mutate(from="eutectic",
                                           "{self$TCol}" := Tabove - X/100,
                                           "{variance_col}" :=sum(wt>0)-1) %>%
                                    unite("uniqueId",all_of(self$pointCol),remove=F) %>%
                                    {.} -> this
                                  
                                  # Add to data
                                  self$phm %<>% bind_rows(this)
                               
                                  # update the list of points
                                  newpoint<-this[1,]%>%select(all_of(c(self$pointCol,"uniqueId")))
                                  self$modPoints %<>% add_row(newpoint)
                                  
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


### TESTS
# stop("The rest is just for testing!\n")
# 
# ACGdio<-PerplexModel$new("ACG Diorite")
# runWerami36("./GrTypes/","ACGdio") %>% .[["phm"]] %>% ACGdio$readPhm(verbose = F)
# ACGdio$renameMinerals("q","Quartz")
# ACGdio$plotPhaseMap(Tunit="K")
# 
# ACGdio$renameMinerals(c("q","Cpx(HGP)","Bi(HGP)"),c("Quartz","Cpx","Bi"))
#                       
# # ACGdio$countDuplicates()
# # ACGdio$classifyMineral()
# # ACGdio$completeTable()
# # ACGdio$refineEutectic()
# # 
#  ACGdio$plotPhaseMap(Tunit="K")
# # 
# # stop()
# # ACGdio$renameMineral("q","Quartz")
# # ACGdio$renameMineral("sph","Quartz")
# # ACGdio$countDuplicates()
# # ACGdio$mergeMinerals("Quartz")

