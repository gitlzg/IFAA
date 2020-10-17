
metaData=function(MicrobData,CovData,linkIDname,testCov=NULL,ctrlCov=NULL,
                  testMany=T,ctrlMany=F,MZILN=F){
  results=list()

  if(length(linkIDname)==0){
    stop("linkIDname is missing.")
  }

  if(length(testCov)>0 | length(ctrlCov)>0){
    if(sum(c(testCov,ctrlCov)%in%colnames(CovData))!=length(c(testCov,ctrlCov))){
      stop("Error: some covariates are not available in the data.")
    }
  }

  if(sum(testCov%in%ctrlCov)>0){
    cat("Warnings: Variables appeared in both testCov list and ctrlCov list will be treated as testCov.","\n")
  }

  # read microbiome data
  if(is.matrix(MicrobData))MdataWithId=data.matrix(MicrobData)
  if(is.data.frame(MicrobData))MdataWithId=data.matrix(MicrobData)
  if(is.character(MicrobData)){
    nCharac=nchar(MicrobData)
    if(substr(MicrobData,(nCharac-2),nCharac)=="csv"){
      MdataWithId=data.matrix(read.csv(file=MicrobData,header=T,na.strings=c("","NA")))
    }
    if(substr(MicrobData,(nCharac-2),nCharac)=="tsv"){
      MdataWithId=data.matrix(read.table(file=MicrobData, sep='\t',header=T,na.strings=c("","NA")))
    }
  }

  if(length(colnames(MdataWithId))!=ncol(MdataWithId))
    stop("Microbiome data lack variable names.")

  missPropMData=sum(is.na(MdataWithId[,linkIDname]))/nrow(MdataWithId)
  if(missPropMData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the Microbiome data file.
             Double check the data format.","\n")
  }

  # read covariate data
  if(is.matrix(CovData))CovarWithId=data.matrix(CovData)
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

  if(length(colnames(CovarWithId))!=ncol(CovarWithId))
    stop("Covariate data lack variable names.")

  missPropCovData=sum(is.na(CovarWithId[,linkIDname]))/nrow(CovarWithId)
  if(missPropCovData>0.8){
    cat("Warning: There are over 80% missing values for the linkId variable in the covariates data file.
             Double check the data format.","\n")
  }

  Covariates1=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=F]

  # determine testCov and ctrlCov
  if(length(testCov)==0){
    if(!testMany){
      stop("No covariates are specified for estimating associations of interest.")
    }else{
      cat("Associations are being estimated for all covariates since no covariates are specified for testCov.","\n")
      testCov=colnames(Covariates1)
    }
  }
  results$testCov=testCov

  ctrlCov=ctrlCov[!ctrlCov%in%testCov]

  xNames=colnames(Covariates1)
  rm(Covariates1)

  if(length(ctrlCov)==0 & ctrlMany){
    cat("No control covariates are specified,
          all variables except testCov are considered as control covariates.","\n")
    ctrlCov=xNames[!xNames%in%testCov]
  }
  results$ctrlCov=ctrlCov[!ctrlCov%in%testCov]

  # merge data to remove missing
  CovarWithId1=CovarWithId[,c(linkIDname,testCov,ctrlCov)]

  allRawData=data.matrix(na.omit(merge(CovarWithId1,MdataWithId,by=linkIDname,all.x=F,all.y=F)))

  CovarWithId=allRawData[,(colnames(allRawData)%in%colnames(CovarWithId1)),drop=F]
  Covariates=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=F]
  rm(CovarWithId1)

  if(!is.numeric(Covariates[,testCov,drop=F])){
    stop("There are non-numeric variables in the covariates for association test.")
  }

  MdataWithId=allRawData[,(colnames(allRawData)%in%colnames(MdataWithId))]
  Mdata_raw=MdataWithId[,!colnames(MdataWithId)%in%linkIDname,drop=F]
  rm(allRawData)

  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads=sum(colSums(Mdata_raw)==0)
  if(numTaxaNoReads>0){
    Mdata_raw=Mdata_raw[,!(colSums(Mdata_raw)==0)]
    cat("There are",numTaxaNoReads,"taxa without any sequencing reads and
        excluded from the analysis","\n")
  }
  rm(numTaxaNoReads)

  numSubNoReads=sum(rowSums(Mdata_raw)==0)
  if(numSubNoReads>0){
    cat("There are",numSubNoReads,"subjects without any sequencing reads and
        excluded from the analysis","\n")
    Mdata_raw=Mdata_raw[!(rowSums(Mdata_raw)==0),]
    MdataWithId=MdataWithId[!(rowSums(Mdata_raw)==0),]
  }
  rm(numSubNoReads)

  Mdata=Mdata_raw
  rm(Mdata_raw)

  microbName=colnames(Mdata)
  newMicrobNames=paste0("microb",seq(length(microbName)))
  results$Mprefix="microb"

  colnames(Mdata)=newMicrobNames
  MdataWithId_new=cbind(MdataWithId[,linkIDname,drop=F],Mdata)
  results$microbName=microbName
  results$newMicrobNames=newMicrobNames
  rm(microbName,newMicrobNames)

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
    results$varNamForBin=xNames[binCheck==2]
    results$BinVars=length(results$varNamForBin)
    for(i in results$varNamForBin){
      mini=min(Covariates[,i],na.rm=T)
      maxi=max(Covariates[,i],na.rm=T)
      if(!(mini==0 & maxi==1)){
        Covariates[Covariates[,i]==mini,i]=0
        Covariates[Covariates[,i]==maxi,i]=1
        cat("Binary covariate",i,"is not coded as 0/1 which may generate analysis bias. It has been changed to 0/1. The changed covariates data can be extracted from the result file.","\n")
      }
    }
    #cat(length(which(binCheck==2)),"binary covariates are detected.","\n")
  }else{
    results$BinVars=0
    binaryInd=NULL
    results$varNamForBin=NULL
  }

  # find the position of the first binary predictor and the rest are all binary
  results$binaryInd=binaryInd
  results$xNames=colnames(Covariates)
  xNewNames=paste0("x",seq(length(xNames)))
  colnames(Covariates)=xNewNames
  results$covsPrefix="x"
  results$xNewNames=xNewNames

  results$testCovInd=which((results$xNames)%in%testCov)
  results$testCovInOrder=results$xNames[results$testCovInd]
  results$testCovInNewNam=results$xNewNames[results$testCovInd]
  rm(xNames,xNewNames)

  CovarWithId_new=cbind(CovarWithId[,linkIDname,drop=F],Covariates)

  data=merge(MdataWithId_new, CovarWithId_new,by=linkIDname,all.x=F,all.y=F)
  results$covariatesData=CovarWithId_new
  colnames(results$covariatesData)=c(linkIDname,results$xNames)
  rm(MdataWithId_new,CovarWithId_new)
  results$data=na.omit(data)
  rm(data)
  cat("Data dimensions (after removing missing data if any):","\n")
  cat(dim(results$data)[1],"samples","\n")
  cat(ncol(Mdata),"OTU's or microbial taxa","\n")

  if(!MZILN)cat(length(results$testCovInOrder),"testCov variables in the analysis","\n")
  if(MZILN)cat(length(results$testCovInOrder),"covariates in the analysis","\n")

  if(length(results$testCovInOrder)>0){
    if(!MZILN)print("These are the testCov variables:")
    if(MZILN)print("These are the covariates:")
    print(testCov)
  }
  rm(testCov)
  if(!MZILN){
    cat(length(results$ctrlCov),"ctrlCov variables in the analysis ","\n")
    if(length(results$ctrlCov)>0){
      print("These are the ctrlCov variables:")
      print(ctrlCov)
    }
    rm(ctrlCov)
  }
  cat(results$BinVars,"binary covariates in the analysis","\n")
  if(results$BinVars>0){
    print("These are the binary covariates:")
    print(results$varNamForBin)
  }
  rm(Mdata,Covariates,binCheck)
  return(results)
}
