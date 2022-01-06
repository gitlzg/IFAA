AIcalcu=function(
  data,
  ref,
  Mprefix,
  covsPrefix
){
  results=list()

  # get the original sample size
  nSub=nrow(data)
  MVarNamLength=nchar(Mprefix)

  # get taxa data
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  nTaxa = length(microPositions)
  nNorm=nTaxa-1
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  # rearrange taxa names
  otherTaxaNames=taxaNames[(taxaNames!=ref)]
  taxaNames=c(otherTaxaNames,ref)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  predNames=colnames(data)[predPositions]
  nPredics=length(predNames)
  rm(predics,predPositions)

  # taxa data
  w=data[,taxaNames]

  # extract x data
  xData=data[,predNames]
  rm(data,predNames)

  # transform data using log-ratio, creat Ai and Li
  l=rep(NA,nSub)
  lLast=rep(NA,nSub)
  taxa.non0=list()
  taxa.0=list()
  logRatiow=list()
  A=list()
  for (i in 1:nSub){
    taxa.nonzero=which(w[i,]!=0)
    lLast[i]=max(taxa.nonzero)
    taxa.zero=which(w[i,]==0)
    taxa.non0[[i]]=w[i,taxa.nonzero]
    taxa.0[[i]]=w[i,taxa.zero]
    if(length(taxa.nonzero)>0){
      last.nonzero=max(taxa.nonzero)
      logwi=as.numeric(log(w[i,taxa.nonzero]))
      l[i]=length(logwi)
      if(l[i]>1){
        logRatiow[[i]]=logwi[1:(l[i]-1)]-logwi[l[i]]
        zero.m=matrix(0,nrow=l[i]-1,ncol=nNorm)
        if(last.nonzero==nTaxa){
          zero.m[cbind(1:(l[i]-1),taxa.nonzero[1:(l[i]-1)])]=1
        } else {
          zero.m[cbind(1:(l[i]-1),taxa.nonzero[1:(l[i]-1)])]=1
          zero.m[,taxa.nonzero[l[i]]]=-1
        }
        A[[i]]=as(zero.m,"sparseMatrix")
        rm(zero.m)
      } else {
        logRatiow[[i]]=NA
        A[[i]]=NA
      }
    } else {
      l[i]=0
      logRatiow[[i]]=NA
      A[[i]]=NA
    }
  }

  # obtain the list of samples whose have at least 2 non-zero taxa
  twoList=which(l>1)
  lengthTwoList=length(twoList)

  rm(w)

  results$xData=xData
  rm(xData)

  results$logRatiow=logRatiow
  rm(logRatiow)
  results$A=A
  rm(A)
  results$twoList=twoList
  rm(twoList)
  results$taxaNames=taxaNames
  rm(taxaNames)
  results$lengthTwoList=lengthTwoList
  results$lLast=lLast
  results$l=l
  results$nTaxa=nTaxa
  results$nNorm=nNorm
  results$nSub=nSub
  results$nPredics=nPredics
  return(results)
}
