##' @export

metaData=function(Microbdata,CovData,linkIDname,testCov=NULL,ctrlCov=NULL,testMany=F,ctrlMany=F){
  results=list()
  
  if(length(linkIDname)==0){
    stop("linkIDname is missing.")
  }
  if(testMany&ctrlMany){
    stop("The options testMany and ctrlMany cannot be TRUE at the same time.")
  }
  if(sum(testCov%in%ctrlCov)>0){
    cat("Warnings: Variables in both testCov and ctrlCov will be treated as testCov.")
  }
  
  if(is.data.frame(Microbdata))MdataWithId=data.matrix(Microbdata)
  if(is.character(Microbdata)){
    nCharac=nchar(Microbdata)
    if(substr(Microbdata,(nCharac-2),nCharac)=="csv"){
      MdataWithId=data.matrix(read.csv(file=Microbdata,header=T,na.strings=c("","NA")))
    }
    if(substr(Microbdata,(nCharac-2),nCharac)=="tsv"){
      MdataWithId=data.matrix(read.table(file=Microbdata, sep='\t',header=T,na.strings=c("","NA")))
    }
  }
  
  missPropMData=sum(is.na(MdataWithId[,linkIDname]))/nrow(MdataWithId)
  if(missPropMData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the Microbiome data file. 
        Double check the data format.","\n")
  }
  Mdata=MdataWithId[,!colnames(MdataWithId)%in%linkIDname,drop=F]
  microbName=colnames(Mdata)
  newMicrobNames=paste0("microb",seq(length(microbName)))
  results$Mprefix="microb"
  
  colnames(Mdata)=newMicrobNames
  MdataWithId_new=cbind(MdataWithId[,linkIDname,drop=F],Mdata)
  results$microbName=microbName
  results$newMicrobNames=newMicrobNames
  rm(microbName,newMicrobNames)
  
  if(is.data.frame(CovData))CovarWithId=data.matrix(CovData)
  if(is.character(CovData)){
    nCharac=nchar(CovData)
    if(substr(CovData,(nCharac-2),nCharac)=="csv"){
      CovarWithId=data.matrix(read.csv(file=CovData,header=T,na.strings=c("","NA")))
    }
    if(substr(CovData,(nCharac-2),nCharac)=="tsv"){
      CovarWithId=data.matrix(read.table(file=CovData, sep='\t',header=T,na.strings=c("","NA")))
    }
  }
  
  missPropCovData=sum(is.na(CovarWithId[,linkIDname]))/nrow(CovarWithId)
  if(missPropCovData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the covariates data file. 
        Double check the data format.","\n")
  }
  
  Covariates=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=F]
  
  if(length(testCov)==0){
    if(!testMany){
      stop("No covariates are specified for estimating associations of interest.")
    }else{
      cat("Associations are being estimated for all covariates since no covariates are specified for testCov.","\n")
      testCov=colnames(Covariates)
    }
  }
  results$testCov=testCov
  
  if(!is.numeric(Covariates[,testCov,drop=F])){
    stop("There are non-numeric variables in the covariates for association test.")
  }
  
  ctrlCov=ctrlCov[!ctrlCov%in%testCov]
  
  xNames=colnames(Covariates)
  if(length(ctrlCov)==0 & ctrlMany){
    cat("No control covariates are specified, 
        all variables except testCov are considered as control covariates.","\n")
    ctrlCov=xNames[!xNames%in%testCov]
  }
  results$ctrlCov=ctrlCov[!ctrlCov%in%testCov]
  
  Covariates=data.matrix(Covariates[,c(testCov,ctrlCov),drop=F])
  
  xNames=colnames(Covariates)
  nCov=length(xNames)
  
  if(sum(is.na(Covariates))>0){
    cat("Samples with missing covariate values are removed from the analysis.","\n")
  }
  
  if(!is.numeric(Covariates[,ctrlCov,drop=F])){
    cat("Warnings: there are non-numeric variables in the control covariates","\n")
    nCtrlCov=length(ctrlCov)
    numCheck=unlist(lapply(seq(nCtrlCov),function(i)is.numeric(Covariates[,ctrlCov[i]])))+0
    for(i in which(numCheck==0)){
      Covariates[,ctrlCov[i]]=as.numeric(factor(Covariates[,ctrlCov[i]]))
    }
  }
  colnames(Covariates)=xNames
  binCheck=unlist(lapply(seq(nCov),function(i)dim(table(Covariates[,xNames[i]]))))
  
  if(length(which(binCheck==2))>0){
    Covariates=Covariates[,c(xNames[binCheck!=2],xNames[binCheck==2]),drop=F]
    binaryInd=length(which(binCheck!=2))+1
    cat(length(which(binCheck==2)),"binary covariates are detected.","\n")
  }else{binaryInd=NULL}
  # find the position of the first binary predictor and the rest are all binary
  results$binaryInd=binaryInd
  
  results$xNames=colnames(Covariates)  
  xNewNames=paste0("x",seq(length(xNames)))
  colnames(Covariates)=xNewNames
  results$covsPrefix="x"
  results$xNewNames=xNewNames
  results$testCovInd=which((results$xNames)%in%testCov)
  results$testCovInOrder=results$xNames[results$testCovInd]
  rm(xNames,xNewNames,testCov,ctrlCov)
  
  CovarWithId_new=cbind(CovarWithId[,linkIDname,drop=F],Covariates)
  
  data=merge(MdataWithId_new, CovarWithId_new,by=linkIDname,all.x=T)
  results$data=na.omit(data)
  rm(data)
  cat("Data dimensions (after removing missing data if any):","\n")
  cat(dim(results$data)[1],"samples","\n")
  cat(ncol(Mdata),"OTU's or microbial taxa","\n")
  cat(length(results$testCovInOrder),"testCov variables","\n")
  cat(length(results$ctrlCov),"ctrlCov variables","\n")
  cat(length(which(binCheck==2)),"binary covariates in both testCov and ctrlCov variables","\n")
  rm(Mdata,Covariates,binCheck)
  return(results)
  }
