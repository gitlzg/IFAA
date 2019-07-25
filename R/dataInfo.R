
##' @export


dataInfo=function(
  data,
  Mprefix,
  covsPrefix,
  binPredInd,
  refReadsThresh=0.1,
  SDThresh=1,
  SDquantilThresh=0.05,
  balanceCut=0.1
){
  results=list()

  # get the original sample size
  nSub=nrow(data)
  MVarNamLength=nchar(Mprefix)
  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  nTaxa = length(microPositions)
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  # get the biggest index number of the taxa name
  taxaNameNum=rep(NA,nTaxa)

  for (i in 1:nTaxa){
    taxa.i.NamLength=nchar(taxaNames[i])
    taxaNameNum[i]=as.numeric(substr(taxaNames[i],(MVarNamLength+1),taxa.i.NamLength))
  }
  maxTaxaNameNum=max(taxaNameNum)
  rm(taxaNameNum)

  w=data[,taxaNames]

  # check zero taxa and subjects with zero taxa reads
  # print("percentage of taxa present for at least 25% subjects:")
  # print(sum((Matrix::colSums(w>0)/nSub)>=0.25))

  results$non0perct=Matrix::colSums(w>0)/nSub

  taxaOverThresh=taxaNames[(Matrix::colSums(w>0)>=nSub*refReadsThresh)]
  if(length(taxaOverThresh)==0){
    stop(paste("There are no taxa with presence over the threshold:",refReadsThresh,
               ". Try lower the reference taxa reads threshold"))
  }

  # check the sd threshold
  sdTaxaOverThresh=rep(0,length(taxaOverThresh))
  for (i in 1:length(taxaOverThresh)){
    taxa.i=w[,taxaOverThresh[i]]
    if(sum(taxa.i>0)>1){
      sdTaxaOverThresh[i]=sd(taxa.i[(taxa.i>0)])
    }
  }

  results$sdTaxa=sdTaxaOverThresh

  TaxaOverSdThresh=taxaOverThresh[(sdTaxaOverThresh>=SDThresh)]
  if(length(TaxaOverSdThresh)==0){
    stop(paste("There are no taxa with SD over the SD threshold:",SDThresh,
               ". Try lower the SD threshold"))
  }
  rm(taxa.i,taxaOverThresh)

  # check the sd quantile threshold
  sdAllTaxa=rep(0,nTaxa)
  for (i in 1:nTaxa){
    taxaAll.i=w[,taxaNames[i]]
    posTaxaAll.i=taxaAll.i[(taxaAll.i>0)]
    if(length(posTaxaAll.i)>1){sdAllTaxa[i]=sd(posTaxaAll.i)}
  }
  goodRefTaxaCandi=TaxaOverSdThresh[(TaxaOverSdThresh>=quantile(sdAllTaxa,probs=SDquantilThresh))]
  rm(sdAllTaxa,posTaxaAll.i,TaxaOverSdThresh)

  if(length(goodRefTaxaCandi)==0){
    stop(paste("There are no taxa with SD over the SD quantile threshold:",SDquantilThresh,
               ". Try lower the SD quantile threshold"))
  }

  nSubNoReads=length(which(Matrix::rowSums(w>0)==0))
  rm(w)

  # get predictor data
  xVarNamLength=nchar(covsPrefix)

  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  predNames=colnames(data)[predPositions]
  EName=predNames[1]
  nPredics=length(predNames)
  rm(predics,predPositions)

  # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
  if(length(binPredInd)>0){
    firstBinPredNam=paste0(covsPrefix,binPredInd)
    binPredStart=which(predNames%in%firstBinPredNam)
    binPredStop=length(predNames)
    allBinPred=predNames[binPredStart:binPredStop]
    nBinPred=length(allBinPred)
    rm(predNames)

    taxaAndBinIndexNoInt=vector()
    taxaNoBin=c()
    taxaBalanceBin=c()

    for(i in 1:nTaxa){
      for(j in 1:nBinPred){
        twoColumns.ij=data[,c(taxaNames[i],allBinPred[j])]
        nNonZero=length(which(twoColumns.ij[,1]>0))
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        rm(twoColumns.ij)
        if(sumOfBin%in%c(0,nNonZero)){
          taxaNoBin=c(taxaNoBin,taxaNames[i])
          index.ij=(i-1)*nPredics+binPredInd+j-1
          taxaAndBinIndexNoInt=c(taxaAndBinIndexNoInt,index.ij)
        }
        if(min(sumOfBin,(nNonZero-sumOfBin))>=balanceCut*nNonZero){
          taxaBalanceBin=c(taxaBalanceBin,taxaNames[i])
        }
      }
    }
    results$taxaAndBinIndexNoInt=taxaAndBinIndexNoInt
    rm(taxaAndBinIndexNoInt,allBinPred,data)


    # remove unbalanced taxa across binary variables
    goodRefTaxaCandi=goodRefTaxaCandi[!(goodRefTaxaCandi%in%taxaNoBin)]

    # keep balanced taxa
    goodRefTaxaCandi=goodRefTaxaCandi[(goodRefTaxaCandi%in%taxaBalanceBin)]

    rm(taxaNoBin)
  }
  # return
  results$taxaNames=taxaNames
  rm(taxaNames)
  results$goodRefTaxaCandi=goodRefTaxaCandi
  rm(goodRefTaxaCandi)
  results$nTaxa=nTaxa
  results$nSub=nSub
  results$nSubNoReads=nSubNoReads
  results$nPredics=nPredics
  results$maxTaxaNameNum=maxTaxaNameNum
  results$EName=EName
  return(results)
}
