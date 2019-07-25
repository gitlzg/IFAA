
##' @export


dataAndInit=function(
  data,
  ref,
  Mprefix,
  covsPrefix,
  binPredInd,
  independence,
  identity,
  treeCorre=F,
  taxaSparsThresh=0,
  lowSDThresh=0,
  lowSDquantilThresh=0
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

  if(length(binPredInd)>0){
    firstBinPredNam=paste0(covsPrefix,binPredInd)
    binPredStart=which(predNames%in%firstBinPredNam)
    allBinPred=predNames[binPredStart:nPredics]
    nBinPred=length(allBinPred)

    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    taxaAndBinIndex=vector()
    for(i in 1:nNorm){
      for(j in 1:nBinPred){
        twoColumns.ij=data[,c(taxaNames[i],allBinPred[j])]
        nNonZero=length(which(twoColumns.ij[,1]>0))
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        if(sumOfBin%in%c(0,nNonZero)){
          index.ij=(i-1)*(nPredics+1)+binPredInd+j
          taxaAndBinIndex=c(taxaAndBinIndex,index.ij)
        }
      }
    }
    results$taxaAndBinIndex=taxaAndBinIndex
    rm(taxaAndBinIndex,allBinPred)
  }

  # check the sd and sd quantile threshold
  sdAllTaxa=rep(0,nNorm)
  for (i in 1:nNorm){
    taxaAll.i=w[,taxaNames[i]]
    posTaxaAll.i=taxaAll.i[(taxaAll.i>0)]
    if(length(posTaxaAll.i)>1){sdAllTaxa[i]=sd(posTaxaAll.i)}
  }
  taxaLowSd=taxaNames[(sdAllTaxa<min(lowSDThresh,quantile(sdAllTaxa,probs=lowSDquantilThresh)))]
  rm(sdAllTaxa,posTaxaAll.i,taxaAll.i)

  # extract x data
  xData=data[,predNames]
  rm(data,predNames)

  # identify taxa with xtreme sparse reads
  xtremeSparsTaxa=taxaNames[(Matrix::colSums(w[,-nTaxa]>0)<taxaSparsThresh)]

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
        A[[i]]=zero.m
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

  # extract data driven covariance matrix
  if(independence & identity){
    results$cov=diag(nNorm)
  }

  rm(w)

  results$xtremeSparsTaxa=xtremeSparsTaxa
  rm(xtremeSparsTaxa)
  results$taxaLowSd=taxaLowSd
  rm(taxaLowSd)
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
  results$nTaxa=nTaxa
  results$nNorm=nNorm
  results$nSub=nSub
  results$nPredics=nPredics

  return(results)
}
