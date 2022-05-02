######## Define melting functions here ################


cat("Melting routines\n")

BatchMelting <- function(self=self,melt_prop,phase_prop,thisPhm,thisGlob){
  
  if(melt_prop==100){
    
    tibble(Name=self$meltName,as_tibble_row(self$c0))  %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    DD <- vector()
    tr_glob_data <- as_tibble_row(DD)
  }else{
    bpm<-BatchPM(kd=self$kd,
                 c0=self$c0,
                 pm=melt_prop,
                 min.props=phase_prop,
                 cmins=matrix(),melt.arg=list(),dont=character(0)) 
    
    cL <- bpm$cL
    cmins <- bpm$cmins
    cS <- bpm$cS
    DD <- bpm$DD
    names(DD)<-paste("D",names(DD),sep="_")
    
    tibble(Name=rownames(cmins),as_tibble(cmins)) %>%
      bind_rows( tibble(Name="solid",as_tibble_row(cS)) ) %>%
      bind_rows( tibble(Name=self$meltName,as_tibble_row(cL)) ) %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    tr_glob_data <- as_tibble_row(DD)
  }
  return(list(tr_point_data=tr_point_data,tr_glob_data=tr_glob_data))
}

BatchMeltingWithZrn <- function(self=self,melt_prop,phase_prop,thisPhm,thisGlob,satModel="Boenhke"){
  
  
  tprK <- thisGlob %>% pull(TK)
  
  if(melt_prop>0){
    # We take the actual compo of the liquid
    liq.maj <- thisPhm %>% filter(Name == self$meltName) %>% select(self$majors)  %>% as_vector
  }else{
    # We use the compo of the coldest liquid
    liq.maj <- self$phm %>% filter(Name==self$meltName,wt>0) %>% filter(TK==min(TK)) %>% select(self$majors)

  }
  
  if(melt_prop==100){
    
    # This will fail if Zrn is a liquidus phase... hopefully it should not happen too often.
    tibble(Name=self$meltName,phWtSys=100,as_tibble_row(self$c0))  %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    DD <- vector()
    tr_glob_data <- as_tibble_row(DD)
  }else{
    bpm<-BatchPM(kd=self$kd,
                 c0=self$c0,
                 pm=melt_prop,
                 min.props=phase_prop,
                 cmins=matrix(),melt.arg=list(),dont=character(0)) 

        zz<-correctZrnSat(kd=self$kd,
                      c0=self$c0,
                      pm=melt_prop,
                      min.props=phase_prop,
                      melt.arg=list(TT=tprK,
                                    mjrs=liq.maj,
                                    trc=bpm$cL
                      ),
                      SatModel=satModel,
                      cmins=bpm$cmins,dont=character(0)
    )
# investigate Zr budget
       
        # if(melt_prop==0){
        # 
        #   print(zz$cmins[,"Zr"])
        #   print(zz$min.props)
        #   print(zz$cmins[,"Zr"] * zz$min.props )
        #   
        # }
        
    #cat("Zrc",zz$min.props["Zrn"])
    if(melt_prop>0){
       cL <- c(zz$cL,phWtSys=melt_prop)
    }else{
      cL<-vector()
     }
    
    cmins <- cbind(zz$cmins,phPropSol=zz$min.props,phWtSys=zz$min.props*(100-melt_prop) )
    cS <- zz$cS
    DD <- zz$DD
    names(DD)<-paste("D",names(DD),sep="_")
    
    tibble(Name=rownames(cmins),as_tibble(cmins)) %>%
      bind_rows( tibble(Name="solid",as_tibble_row(cS)) ) %>%
      bind_rows( tibble(Name=self$meltName,as_tibble_row(cL)) ) %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    tr_glob_data <- as_tibble_row(c(DD,zz$sat))
  }
  return(list(tr_point_data=tr_point_data,tr_glob_data=tr_glob_data))
}


BatchMeltingWithZrnKdmodOkz <- function(self=self,melt_prop,phase_prop,thisPhm,thisGlob,satModel="Boenhke",KdtoChange=NULL,Kdfn=NULL){
  
  tprK <- thisGlob %>% pull(TK)
  
  if(melt_prop>0){
    # We take the actual compo of the liquid
    liq.maj <- thisPhm %>% filter(Name == self$meltName) %>% select(self$majors)  %>% as_vector
  }else{
    # We use the compo of the coldest liquid
    liq.maj <- self$phm %>% filter(Name==self$meltName,wt>0) %>% filter(TK==min(TK)) %>% select(self$majors)
  }
  
  # Tweak Kd
  # Works for T-fitted Kd
  
  modifiedKd <- Kdfn(tprK)
  self$kd[KdtoChange,names(modifiedKd)]<-modifiedKd
  
  KdChanges <- as_tibble_row( c(KdtoChange=KdtoChange,modifiedKd) )
  
  if(melt_prop==100){
    
    # This will fail if Zrn is a liquidus phase... hopefully it should not happen too often.
    tibble(Name=self$meltName,phWtSys=100,as_tibble_row(self$c0))  %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    DD <- vector()
    tr_glob_data <- as_tibble_row(DD)
  }else{
    bpm<-BatchPM(kd=self$kd,
                 c0=self$c0,
                 pm=melt_prop,
                 min.props=phase_prop,
                 cmins=matrix(),melt.arg=list(),dont=character(0)) 
    
    zz<-correctZrnSat(kd=self$kd,
                      c0=self$c0,
                      pm=melt_prop,
                      min.props=phase_prop,
                      melt.arg=list(TT=tprK,
                                    mjrs=liq.maj,
                                    trc=bpm$cL
                      ),
                      SatModel=satModel,
                      cmins=bpm$cmins,dont=character(0)
    )
    
    if(melt_prop>0){
      cL <- c(zz$cL,phWtSys=melt_prop)
    }else{
      cL<-vector()
    }
    
    cmins <- cbind(zz$cmins,phPropSol=zz$min.props,phWtSys=zz$min.props*(100-melt_prop) )
    cS <- zz$cS
    DD <- zz$DD
    names(DD)<-paste("D",names(DD),sep="_")
    
    tibble(Name=rownames(cmins),as_tibble(cmins)) %>%
      bind_rows( tibble(Name="solid",as_tibble_row(cS)) ) %>%
      bind_rows( tibble(Name=self$meltName,as_tibble_row(cL)) ) %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    tr_glob_data <- as_tibble_row(c(DD,zz$sat,KdChanges=nest(KdChanges)))
  }
  return(list(tr_point_data=tr_point_data,tr_glob_data=tr_glob_data))
}


BatchMeltingWithZrnKdmod <- function(self=self,melt_prop,phase_prop,thisPhm,thisGlob,satModel="Boenhke",KdtoChange=NULL,Kdfn=NULL){
 
  tprK <- thisGlob %>% pull(TK)
  
  #localKd <- self$kd
  
  if(melt_prop>0){
    # We take the actual compo of the liquid
    liq.maj <- thisPhm %>% filter(Name == self$meltName) %>% select(self$majors)  %>% as_vector
  }else{
    # We use the compo of the coldest liquid
    liq.maj <- self$phm %>% filter(Name==self$meltName,wt>0) %>% filter(TK==min(TK)) %>% select(self$majors)
    }
  
  # Tweak Kd
  # Works for T-fitted Kd
  
  map(Kdfn,exec,TK=tprK) %>%  # We execute each fn in Kdfn
      map(as_tibble_row) %>% 
      reduce(bind_rows) %>%  # collapse the list to tibble
      bind_cols(KdtoChange=KdtoChange) %>% # Add the name of the mineral that was changed
      {.} -> KdChanges
  
  ## Testing
  # Kdfn <- list(KdKirklandZrnT,KdClaiborneZrnT)
  # KdtoChange=c("Zrn","Zrn1")
  # Implement the changes now !
  
  .tweak<- function(KdtoChange ,Element ,Value ){
    self$kd[KdtoChange,Element]<-Value
  }
  KdChanges %>% gather(key="Element",value="Value",-KdtoChange) %>% drop_na() %>% 
       pwalk(.tweak)
  
 
    if(melt_prop==100){
    
    # This will fail if Zrn is a liquidus phase... hopefully it should not happen too often.
    tibble(Name=self$meltName,phWtSys=100,as_tibble_row(self$c0))  %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data
    
    DD <- vector()
    tr_glob_data <- as_tibble_row(DD)
  }else{
    bpm<-BatchPM(kd=self$kd,
                 c0=self$c0,
                 pm=melt_prop,
                 min.props=phase_prop,
                 cmins=matrix(),melt.arg=list(),dont=character(0)) 
    
    zz<-correctZrnSat(kd=self$kd,
                      c0=self$c0,
                      pm=melt_prop,
                      min.props=phase_prop,
                      melt.arg=list(TT=tprK,
                                    mjrs=liq.maj,
                                    trc=bpm$cL
                      ),
                      SatModel=satModel,
                      cmins=bpm$cmins,dont=character(0)
    )

    if(melt_prop>0){
      cL <- c(zz$cL,phWtSys=melt_prop)
    }else{
      cL<-vector()
    }
    
    cmins <- cbind(zz$cmins,phPropSol=zz$min.props,phWtSys=zz$min.props*(100-melt_prop) )
    cS <- zz$cS
    DD <- zz$DD
    names(DD)<-paste("D",names(DD),sep="_")
    
    tibble(Name=rownames(cmins),as_tibble(cmins)) %>%
      bind_rows( tibble(Name="solid",as_tibble_row(cS)) ) %>%
      bind_rows( tibble(Name=self$meltName,as_tibble_row(cL)) ) %>%
      bind_rows( tibble(Name="system",as_tibble_row(self$c0)) ) %>%
      {.} -> tr_point_data

    tr_glob_data <- as_tibble_row(c(DD,zz$sat,KdChanges=nest(KdChanges)))
  }
  cat(self$kd["Zrn","Th"],"..")
  return(list(tr_point_data=tr_point_data,tr_glob_data=tr_glob_data))
}

