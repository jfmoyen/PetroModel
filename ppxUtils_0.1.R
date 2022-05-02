### Clean a werami tab file with phase proportions
# WARNING # in case of solvus, need manual cleaning - do NOT use

translateHeader<-function(tab,Glossary){
  #' @param tab: the table to convert
  #' @param glossary: conversion of perpleX names into "clean" names
  
  hdr <- colnames(tab)
  for(i in 1:length(hdr)){
    if(hdr[i] %in% Glossary$oldNames){
      old<-hdr[i]
      new<-as.character(Glossary$newNames[Glossary$oldNames==old])
      hdr[i]<-new
      cat(old,"-->",new)
    }else{
      cat(hdr[i],"not replaced")
    }
    cat("\n")
  }
  
  colnames(tab)<-hdr
  return(tab)
}

cleanWerFile<-function(basename,number=1,glossary="",firstLine=9,rnPrefix=""){
  #' @param basename: fully qualified, base name for files
  #' @param number: order number of tab-file
  #' @param glossary: conversion of perpleX names into "clean" names
  #' @param firstLine: line where the real data starts (in perpleX 6.9.0, line 9)
  #' @param rnPrefix prefix for the row names
  
  # Build the full file name
  tabfile<-paste(basename,number,sep="_")
  tabfile<-paste(tabfile,".tab",sep="")
  
  # Actually read the file
  ee<-scan(tabfile,what=character(),blank.lines.skip=F,sep="\n")
  
  ####  Process the headers
  hdr<-unlist(strsplit(ee[firstLine]," +"))
  
  # Warn if duplicate names are found
  if(any(table(hdr)>1)){
    dups<-which(table(hdr)>1)
    cat("WARNING: duplicate columns found;",names(dups),"\nExpect trouble !")
  }
  
  # if a glossary was supplied, attempt to convert the col names
  if(glossary!=""){
    for(i in 1:length(hdr)){
      if(hdr[i] %in% glossary$oldNames){
        old<-hdr[i]
        new<-as.character(glossary$newNames[glossary$oldNames==old])
        hdr[i]<-new
        cat(old,"-->",new)
      }else{
        cat(hdr[i],"not replaced")
      }
      cat("\n")
    }
  }
  
  #### Now process the real data
  qq<-strsplit(ee[(firstLine+1):length(ee)]," +")
  foo<-array(unlist(qq),dim=c(length(hdr)+1,length(qq)))
  foo<-foo[-1,]
  foo[foo=="NaN"]<-0
  foo<-apply(foo,1,as.numeric)
  
  colnames(foo)<-hdr
  rownames(foo)<-paste(rnPrefix,(1:nrow(foo)),sep="")
  
  return(foo)
}


runVertex<-function(prefix,datfile,com0){
  # prefix: path from vertex dir to file
  # datfile: name of .dat file
  # com0: commands passed to vertex
  
  cat("Vertex: processing",datfile,"...")
  ee<-system("vertex.exe",intern=T,
             input=c(paste(prefix,datfile,sep=""),com0))
  o
  save(ee,file=paste(prefix,datfile,"_vertexOutput.txt",sep=""),ascii=T)
  cat("...done\n")
}