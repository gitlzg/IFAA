
##' @export


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
      logwi=as(as.numeric(log(w[i,taxa.nonzero])),"sparseVector")
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

##' @export

allUserFunc=function(){
  return(c("dataRecovTrans", "AIcalcu","metaData",
           "dataSparsCheck","dataInfo",
           "cvPicasso","runPicasso","groupBetaToFullBeta"))
}
##' @export



bootResuHDCI=function(
  data,
  refTaxa,
  bootB,
  bootLassoAlpha,
  binPredInd,
  covsPrefix,
  Mprefix,
  standardize,
  seed
){
  
  results=list()
  # load data info
  basicInfo=dataInfo(data=data,binPredInd=binPredInd,
                     covsPrefix=covsPrefix,Mprefix=Mprefix)
  taxaNames=basicInfo$taxaNames
  ii=which(basicInfo$taxaNames%in%refTaxa)
  
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # inital Lasso OLS estimate
  dataForEst=dataRecovTrans(data=data,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix)
  
  x=as(dataForEst$xTildalong,"sparseMatrix")
  y=as(dataForEst$UtildaLong,"sparseVector")
  rm(dataForEst)
  
  penal=runBootLassoHDCI(x=x,y=y,nPredics=nPredics,nTaxa=nTaxa,
                         refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,
                         bootB=bootB,standardize=standardize,
                         seed=seed)
  rm(x,y)
  
  results$finalBetaEst=penal$beta
  results$CIvecLow=penal$betaCIlow
  results$CIvecUp=penal$betaCIhi
  
  results$finalBetaEst.LPR=penal$beta.LPR
  results$CIvecLow.LPR=penal$betaCIlow.LPR
  results$CIvecUp.LPR=penal$betaCIhi.LPR
  
  rm(penal)
  return(results)
}



##' @export


cvPicasso=function(
  x,
  y,
  family,
  lamMax,
  lambda,
  nfolds,
  zeroSDCut,
  lambda.min.ratio,
  nLam,
  method,
  standardize,
  seed,
  seedi
){
  
  results=list()
  
  nObsAll=length(y)
  nBeta=ncol(x)
  
  # obtain lambda range
  if(length(lambda)==0){
    #lamMin=lambda.min.ratio*lamMax
    lamMin=0
    lamList=seq(lamMax,lamMin,length=(nLam+1))[1:nLam]
  }
  if(length(lambda)>0){
    lamList=lambda
  }
  
  # partition the data randomly into nfolds partitions
  parSize=floor(nObsAll/nfolds)
  
  if(length(seed)>0)set.seed(as.numeric(seed)+10^7+seedi)
  
  randomShuf=sample(nObsAll, nObsAll)
  sampleInd=list()
  for (i in 1:nfolds){
    if(i<nfolds){
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):(i*parSize)]
    }else{
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):nObsAll]
    }
  }
  rm(randomShuf)
  
  # cross validation
  #sink(paste("picassoPredict",method,".txt",sep="")) # to avoid output from the predict() function
  
  cvPara=matrix(NA,nrow=nLam,ncol=nfolds)
  
  for(i in 1:nfolds){
    # check if there is zero-variance x in the partitions
    startT.i=proc.time()[3]
    xi=x[-sampleInd[[i]],]
    nObs.i=nrow(xi)
    sdX.i=apply(xi,2,sd)
    xWithNearZeroSd.i=which(sdX.i<=zeroSDCut)
    rm(sdX.i)
    
    # remove near constant columns in x
    if(length(xWithNearZeroSd.i)>0) {
      xi=xi[,-xWithNearZeroSd.i]
    }
    nearZeroSd.i=length(xWithNearZeroSd.i)
    
    yi=as(y[-sampleInd[[i]]],"sparseVector")
    
    cv.i=picasso(X=xi,Y=yi,lambda=lamList,method=method,family=family,standardize=standardize)
    rm(xi,yi)
    
    if(length(xWithNearZeroSd.i)>0) {
      missPositions=xWithNearZeroSd.i
      for(iLam in 2:nLam){
        missPositions=c(missPositions,xWithNearZeroSd.i+(iLam-1)*nBeta)
      }
      rm(xWithNearZeroSd.i)
      gc()
      betaa=as.vector(cv.i$beta)
      betaTrans=groupBetaToFullBeta(nTaxa=(nBeta*nLam),nPredics=1,
                                    unSelectList=sort(missPositions),newBetaNoInt=betaa)
      rm(betaa)
      gc()
      betaMatrix=Matrix(betaTrans$finalBeta,nrow=nBeta,ncol=nLam)
      rm(betaTrans)
    } else {betaMatrix=cv.i$beta}
    
    # training is done, start testing
    if(sum(is.na(cv.i$intercept))>0){
      cv.i$intercept[is.na(cv.i$intercept)]=mean(cv.i$intercept[!is.na(cv.i$intercept)])
    }
    if(sum(is.na(cv.i$intercept))>nLam*0.15){
      warnings("Many intercept estimates are missing in picasso estimates")
    }
    intcep=rep(1,length(sampleInd[[i]]))%*%t(as.matrix(cv.i$intercept))
    rm(cv.i)
    yHatTest.i=x[sampleInd[[i]],]%*%betaMatrix+intcep
    rm(betaMatrix,intcep)
    gc()
    
    resiVecs.i=(yHatTest.i-y[sampleInd[[i]]])
    rm(yHatTest.i)
    cvPara[,i]=apply(resiVecs.i,2,function(x)sum(x^2))
    rm(resiVecs.i)
  }
  rm(sampleInd,x,y)
  
  SSE=Matrix::rowSums(cvPara)
  rm(cvPara)
  
  optiLamInd=tail(which(SSE==min(SSE)),n=1)
  optiLam=lamList[optiLamInd]
  
  # return results
  results$SSE=SSE
  rm(SSE)
  
  results$lamList=lamList
  rm(lamList)
  
  results$optiLam=optiLam
  rm(optiLam)
  
  results$optiLamInd=optiLamInd
  rm(optiLamInd)
  
  results$nLam=nLam
  rm(nLam)
  return(results)
}



##' @export




dataInfo=function(
  data,
  Mprefix,
  covsPrefix,
  qualifyRefTax=F,
  binPredInd,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut
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
  
  if(qualifyRefTax){
    qualifyData=data[rowSums(data[,taxaNames]>0)>=2,,drop=F]
    w=qualifyData[,taxaNames]
    nSubQualif=nrow(qualifyData)
    taxaOverThresh=taxaNames[(Matrix::colSums(w>0)>=(nSubQualif*refReadsThresh))]
    if(length(taxaOverThresh)==0){
      cat("There are no taxa with presence over the threshold:",refReadsThresh,
          ". Try lower the reference taxa reads threshold.","\n")
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
      cat("There are no taxa with SD over the SD threshold:",SDThresh,
          ". Try lower the SD threshold","\n")
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
      cat("There are no taxa with SD over the SD quantile threshold:",SDquantilThresh,
          ". Try lower the SD quantile threshold","\n")
    }
    rm(w)
  }
  
  # get predictor data
  xVarNamLength=nchar(covsPrefix)
  
  predics = sapply(substr(colnames(data),1,xVarNamLength), function(x) {grep(covsPrefix, x)})
  predPositions=which(predics == 1)
  predNames=colnames(data)[predPositions]
  nPredics=length(predNames)
  rm(predics,predPositions,data)
  
  if(qualifyRefTax){
    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    if(length(binPredInd)>0){
      allBinPred=predNames[binPredInd:nPredics]
      nBinPred=length(allBinPred)
      
      taxaBalanceBin=c()
      
      for(i in 1:nTaxa){
        for(j in 1:nBinPred){
          twoColumns.ij=qualifyData[,c(taxaNames[i],allBinPred[j])]
          nNonZero=length(which(twoColumns.ij[,1]>0))
          sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
          rm(twoColumns.ij)
          if(min(sumOfBin,(nNonZero-sumOfBin))>=floor(balanceCut*nSubQualif)){
            taxaBalanceBin=c(taxaBalanceBin,taxaNames[i])
          }
        }
      }
      
      rm(allBinPred,qualifyData)
      
      # keep balanced taxa
      goodRefTaxaCandi=goodRefTaxaCandi[(goodRefTaxaCandi%in%taxaBalanceBin)]
    }
    results$goodRefTaxaCandi=goodRefTaxaCandi
    rm(goodRefTaxaCandi)
  }
  # return 
  results$taxaNames=taxaNames
  rm(taxaNames)
  results$predNames=predNames
  rm(predNames)
  results$nTaxa=nTaxa
  results$nSub=nSub
  results$nPredics=nPredics
  return(results)
}


##' @export



dataRecovTrans=function(
  data,
  ref,
  Mprefix,
  covsPrefix,
  xOnly=F,
  yOnly=F
){
  results=list()
  
  # load A and log-ratio transformed RA
  data.and.init=AIcalcu(data=data,ref=ref,Mprefix=Mprefix,
                        covsPrefix=covsPrefix)
  rm(data)
  
  taxaNames=data.and.init$taxaNames
  A=data.and.init$A
  logRatiow=data.and.init$logRatiow
  nSub=data.and.init$nSub
  nTaxa=data.and.init$nTaxa
  xData=data.and.init$xData
  nPredics=data.and.init$nPredics
  twoList=data.and.init$twoList
  lLast=data.and.init$lLast
  L=data.and.init$l
  lengthTwoList=data.and.init$lengthTwoList
  rm(data.and.init) 
  
  nNorm=nTaxa-1
  xDimension=nPredics+1 # predictors+intercept
  
  # create omegaRoot
  omegaRoot=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if(lLast[i]==nTaxa){
      #omega[[i]]=Diagonal(L[i]-1)
      omegaRoot[[i]]=Diagonal(L[i]-1)
    } else{ if(L[i]==2){
      #omega[[i]]=1/2
      omegaRoot[[i]]=sqrt(0.5)
    }else{
      dim=L[i]-1
      a=(1+(dim-2)/2)/(1/2*(1+(dim-1)/2))
      b=-1/(1+(dim-1)/2)
      #omega[[i]]=1/2*((a-b)*Diagonal(dim)+b*matrix(1,nrow=dim,ncol=dim))
      
      # calculate the square root of omega assuming it is exchangeable
      aStar=dim^2/((dim-1)^2)
      bStar=b*(dim-2)/(dim-1)-a*(dim^2-2*dim+2)/((dim-1)^2)
      cStar=(0.5*b-0.5*a*(dim-2)/(dim-1))^2
      cSquare=(-bStar+sqrt(bStar^2-4*aStar*cStar))/(2*aStar)
      if(cSquare<0) {stop("no solution for square root of omega")}
      d=sqrt((0.5*a-cSquare)/(dim-1))
      if (is.na(d)){stop("no solution for off-diagnal elements for square root of omega")}
      c=(0.5*b-(dim-2)*(d^2))/(2*d)
      omegaRoot[[i]]=-((c-d)*Diagonal(dim)+d*matrix(1,nrow=dim,ncol=dim))
    }
    }
  }
  rm(L,lLast)
  
  if(xOnly){
    # create X_i in the regression equation using Kronecker product
    xDataWithInter=as.matrix(cbind(rep(1,nSub),xData))
    rm(xData)
    
    xInRegres=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      xInRegres[[i]]=as(t(as.matrix(kronecker(Diagonal(nNorm),xDataWithInter[i,]))),"sparseMatrix")
    } 
    rm(xDataWithInter)
    
    xDataTilda=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      xDataTilda[[i]]=omegaRoot[[i]]%*%A[[i]]%*%xInRegres[[i]]
    }
    rm(omegaRoot,logRatiow,xInRegres)
    
    # stack xTilda
    for (j in 1:lengthTwoList){
      i=twoList[j]
      if (j==1) {xTildalong=xDataTilda[[i]]
      } else {
        xTildalong=rbind(xTildalong,xDataTilda[[i]])
      }
    }
    
    results$xTildalong=xTildalong
    rm(xTildalong)
    return(results)
  }
  
  if(yOnly){
    Utilda=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      Utilda[[i]]=omegaRoot[[i]]%*%logRatiow[[i]]
    }
    rm(omegaRoot,logRatiow)
    
    # stack Utilda
    for (j in 1:lengthTwoList){
      i=twoList[j]
      if (j==1) {UtildaLong=Utilda[[i]]
      } else {
        UtildaLong=rbind(UtildaLong,Utilda[[i]])
      }
    }
    
    rm(Utilda)
    results$UtildaLong=UtildaLong
    rm(UtildaLong)
    return(results)
  }
  
  # create X_i in the regression equation using Kronecker product
  xDataWithInter=data.matrix(cbind(rep(1,nSub),xData))
  colnames(xDataWithInter)[1]="Inter"
  rm(xData)
  
  xInRegres=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    xInRegres[[i]]=as(t(as.matrix(kronecker(Diagonal(nNorm),xDataWithInter[i,]))),"sparseMatrix")
  } 
  
  rm(xDataWithInter)
  
  # transfer U and xInRegres using omega root
  xDataTilda=list()
  Utilda=list()
  
  for (j in 1:lengthTwoList){
    i=twoList[j]
    Utilda[[i]]=omegaRoot[[i]]%*%logRatiow[[i]]
    xDataTilda[[i]]=omegaRoot[[i]]%*%A[[i]]%*%xInRegres[[i]]
  }
  rm(omegaRoot,logRatiow,xInRegres)
  
  # stack Utilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {UtildaLong=Utilda[[i]]
    } else {
      UtildaLong=rbind(UtildaLong,Utilda[[i]])
    }
  }
  rm(Utilda)
  
  # stack xTilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {xTildalong=xDataTilda[[i]]
    } else {
      xTildalong=rbind(xTildalong,xDataTilda[[i]])
    }
  }
  rm(xDataTilda,lengthTwoList)
  
  # return objects 
  results$UtildaLong=UtildaLong
  rm(UtildaLong)
  results$xTildalong=xTildalong
  rm(xTildalong)
  
  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}


##' @export



dataSparsCheck=function(
  data,
  Mprefix
){
  results=list()
  
  # get the original sample size
  nSub=nrow(data)
  MVarNamLength=nchar(Mprefix)
  
  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)
  
  taxaNames=colnames(data)[microPositions]
  rm(microPositions)
  
  w=data[,taxaNames,drop=F]
  rm(data,taxaNames)
  overallSparsity=round(100*sum(w==0)/(nrow(w)*ncol(w)),2)
  cat(overallSparsity,"percent of microbiome sequencing reads are zero","\n")
  
  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads=length(which(Matrix::colSums(w)==0))
  if(numTaxaNoReads>0){
    cat("There are",numTaxaNoReads,"taxa without any sequencing reads and 
        excluded from the analysis","\n")
  }
  rm(numTaxaNoReads)
  
  numSubNoReads=length(which(Matrix::rowSums(w)==0))
  if(numSubNoReads>0){
    cat("There are",numSubNoReads,"subjects without any sequencing reads and 
        excluded from the analysis","\n")
  }
  rm(numSubNoReads,w)
}

##' @export



getScrResu=function(
  method,
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nPermu,
  doPermut=T,
  x1permut,
  nRef,
  paraJobs,
  refTaxa=NULL,
  fwerRate,
  goodIndeCutPerc=0.33,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed
){
  results=list()
  
  # run permutation 
  scrParal=runScrParal(data=data,testCovInd=testCovInd,
                       testCovInOrder=testCovInOrder,
                       testCovInNewNam=testCovInNewNam,nRef=nRef,
                       nPermu=nPermu,doPermut=doPermut,
                       x1permut=x1permut,paraJobs=paraJobs,
                       refTaxa=refTaxa,method=method,
                       standardize=standardize,
                       sequentialRun=sequentialRun,
                       allFunc=allFunc,
                       refReadsThresh=refReadsThresh,SDThresh=SDThresh,
                       SDquantilThresh=SDquantilThresh,balanceCut=balanceCut,
                       Mprefix=Mprefix,covsPrefix=covsPrefix,
                       binPredInd=binPredInd,seed=seed)
  
  selecCountOverall=scrParal$countOfSelecForAPred
  selecCountMatIndv=scrParal$testCovCountMat
  taxaNames=scrParal$taxaNames
  goodRefTaxaCandi=scrParal$goodRefTaxaCandi
  
  nTaxa=scrParal$nTaxa
  nPredics=scrParal$nPredics
  nTestCov=scrParal$nTestCov
  binomPar=scrParal$binomPar
  maxVec=scrParal$maxVec
  results$maxVec=maxVec
  MaxMatTestCovByPermu=scrParal$MaxMatTestCovByPermu
  results$refTaxa=scrParal$refTaxa
  rm(scrParal)
  gc()
  
  if(doPermut){
    # control family wise error rate
    originFwerCut=quantile(maxVec,probs=(1-fwerRate))
    fwerCut=originFwerCut
    
    results$twoMeanUsed=0
    
    if(fwerCut<=min(selecCountOverall)){
      # use two mean clustering
      twoMeanClusResu=twoMeanClus(matrix=selecCountOverall,nRef=nRef)
      fwerCut=twoMeanClusResu$twoMeanCut
      rm(twoMeanClusResu)
      if(fwerCut==0)fwerCut=0.1
      results$twoMeanUsed=1
    }
    results$fwerCut=fwerCut
    
    results$selecTaxaFWER=as(as.vector((selecCountOverall>=fwerCut)+0),"sparseVector")
    rm(maxVec)
    
    if(nTestCov==1){
      results$selecMatIndv=matrix(results$selecTaxaFWER,nrow=1)
      results$fwerCutIndv=fwerCut
    }
    
    if(nTestCov>1){
      fwerCutIndv=vector()
      for(i in 1:nTestCov){
        originFwerCut.i=quantile(MaxMatTestCovByPermu[i,],probs=(1-fwerRate))
        fwerCutIndv[i]=originFwerCut.i
        if(originFwerCut.i<=min(selecCountMatIndv[i,])){
          # use two mean clustering
          twoMean=twoMeanClus(matrix=selecCountMatIndv[i,,drop=F],nRef=nRef)
          fwerCutIndv[i]=twoMean$twoMeanCut
          if(twoMean$twoMeanCut==0)fwerCutIndv[i]=0.1
          rm(twoMean,originFwerCut.i)
        }
      }
      results$selecMatIndv=as((selecCountMatIndv>=fwerCutIndv)+0,"sparseMatrix")
      selecIndvInOverall_trans=t(as.matrix(results$selecMatIndv))*as.vector(results$selecTaxaFWER)
      results$selecIndvInOverall=as(t(as.matrix(selecIndvInOverall_trans)),"sparseMatrix")
      rm(selecIndvInOverall_trans)
      
      results$selecCountMatIndv=as(selecCountMatIndv,"sparseMatrix")
      results$fwerCutIndv=as(fwerCutIndv,"sparseVector")
      results$MaxMatTestCovByPermu=as(MaxMatTestCovByPermu,"sparseMatrix")
      rm(selecCountMatIndv,fwerCutIndv,MaxMatTestCovByPermu)
    }
    
    goodIndpCut=quantile(selecCountOverall[1,(selecCountOverall[1,]<=fwerCut)],prob=goodIndeCutPerc)
    taxaLessGoodCut=selecCountOverall[1,(selecCountOverall[1,]<=goodIndpCut)]
    taxaLessFWERCut=selecCountOverall[1,(selecCountOverall[1,]<fwerCut)]
    results$taxaNames=taxaNames
    results$selecCountOverall=selecCountOverall
    rm(taxaNames,selecCountOverall)
    
    goodIndpRefTaxWithCount=taxaLessGoodCut[(names(taxaLessGoodCut)%in%goodRefTaxaCandi)]
    if(length(goodIndpRefTaxWithCount)==0){
      results$goodIndpRefTaxLeastCount=NULL
    } else {
      results$goodIndpRefTaxLeastCount=names(tail(goodIndpRefTaxWithCount[goodIndpRefTaxWithCount==min(goodIndpRefTaxWithCount)],n=1))
    }
    
    results$goodIndpRefTaxWithCount=goodIndpRefTaxWithCount
    goodIndpRefTaxFWERcut=taxaLessFWERCut[(names(taxaLessFWERCut)%in%goodRefTaxaCandi)]
    results$goodIndpRefTaxFWERcut=goodIndpRefTaxFWERcut
    if(length(goodIndpRefTaxFWERcut)==0){
      results$goodIndpRefTaxFWERcutLeastCount=NULL
    } else {
      results$goodIndpRefTaxFWERcutLeastCount=names(tail(goodIndpRefTaxFWERcut[goodIndpRefTaxFWERcut==min(goodIndpRefTaxFWERcut)],n=1))
    }
    results$goodIndpRefTaxFWERcut=goodIndpRefTaxFWERcut
    results$goodRefTaxaCandi=goodRefTaxaCandi
    rm(goodRefTaxaCandi)
    
    minPosGoodCut=NULL
    if(sum(taxaLessGoodCut>0)>0){
      minPosGoodCut=min(taxaLessGoodCut[taxaLessGoodCut>0])
    }
    
    minPosFWERCut=NULL
    if(sum(taxaLessFWERCut>0)>0){
      minPosFWERCut=min(taxaLessFWERCut[taxaLessFWERCut>0])
    }
    results$taxaLessGoodCut=taxaLessGoodCut
    results$taxaLessFWERCut=taxaLessFWERCut
    
    results$minPosGoodCut=minPosGoodCut
    results$minPosFWERCut=minPosFWERCut
    
    results$refTaxonQualified=1
    if(length(results$goodIndpRefTaxLeastCount)==1){
      results$finalIndpRefTax=results$goodIndpRefTaxLeastCount
    } else if (length(results$goodIndpRefTaxFWERcutLeastCount)==1){
      results$finalIndpRefTax=results$goodIndpRefTaxFWERcutLeastCount
    } else if(length(minPosGoodCut)>0){
      results$finalIndpRefTax=names(tail(taxaLessGoodCut[taxaLessGoodCut==minPosGoodCut],n=1))
      results$refTaxonQualified=0
    } else if(length(minPosFWERCut)>0){
      results$finalIndpRefTax=names(tail(taxaLessFWERCut[taxaLessFWERCut==minPosFWERCut],n=1))
      results$refTaxonQualified=0
    }else{
      taxaLessGoodCutLeastCount=taxaLessGoodCut[taxaLessGoodCut==min(taxaLessGoodCut)]
      results$finalIndpRefTax=names(taxaLessGoodCutLeastCount)[1]
      results$refTaxonQualified=0
    }
    rm(taxaLessGoodCut,taxaLessFWERCut)
  }
  
  if(!doPermut){results$finalIndpRefTax=refTaxa}
  return(results)
}

##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of covariates of microbiome
##'
##' The `IFAA()` uses a novel approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem.
##'
##' @param MicrobData Microbiome data matrix containing microbiome abundance with each row per sample and each column per taxon/OTU/ASV. It should contain an `"id"` variable to correspond to the `"id"` variable in the covariates data: `CovData`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row per sample and each column per variable. It should also contain an `"id"` variable to correspond to the `"id"` variable in the microbiome data: `MicrobData`.
##' @param linkIDname Variable name of the `"id"` variable in both `MicrobData` and `CovData`. The two data sets will be merged by this `"id"` variable.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `ctrlCov` is not `NULL`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nPermu The number of permutation used in phase 1. Default number is `40`.
##' @param x1permut This takes a logical value `TRUE` or `FALSE`. If true, it will permute the variables in testCov. If false, it will use residual-permutation proposed by Freedman and Lane (1983). 
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param reguMethod regularization approach used in phase 1 of the algorithm. Default is `"mcp"`. Other methods are under development.
##' @param fwerRate The family wise error rate for identifying taxa/OTU/ASV associated with `testCov` in phase 1. Default is `0.25`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Sometimes parallel jobs can not be successfully run for unknown reasons. For example, socket related errors may pop up or some slave cores return simple error instead of numerical results. In those scenarios, setting `sequentialRun = TRUE` may help, but it will take more time to run. Default is `FALSE`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. Default is `8`. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, all design matrix X in phase 1 and phase 2 will be standardized in the analyses. Default is `FALSE`.
##' @param nRefMaxForEsti The maximum number of reference taxa used in phase 2. The default is `1`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2. The default is `500`.
##' @param bootLassoAlpha The significance level in phase 2. Default is `0.05`.
##' @param refReadsThresh The threshold of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param SDThresh The threshold of standard deviations of sequencing reads for choosing the reference taxon in phase 2. The default is `0.5` which means the standard deviation of sequencing reads should be at least `0.5`.
##' @param balanceCut The threshold of non-zero sequencing reads in each group of a binary variable for choosing the reference taxon in phase 2. The default number is `0.2` which means at least 20% sequencing reads are non-zero in each group.
##' @param seed Random seed for reproducibility. Default is `1`.
##' @return A list containing the estimation results.
##'
##' - `analysisResults$estByCovList`: A list containing estimating results for all the variables in `testCov`. See details.
##'
##' - `covariatesData`: A dataset containing covariates and confounders used in the analyses.
##'
##' @examples
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' results <- IFAA(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"), nRef = 4,
##'                 nPermu = 4,
##'                 fwerRate = 0.25,
##'                 bootB = 5)
##'
##' @references Li et al.(2020) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. arXiv:1909.10101v3
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Freedman and Lane (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics. 1(4):292-298.

##' @export
##' @md





IFAA=function(
  MicrobData,
  CovData,
  linkIDname,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=T,
  ctrlMany=F,
  nRef=40,
  nRefMaxForEsti=1,
  nPermu=40,
  x1permut=T,
  refTaxa=NULL,
  reguMethod=c("mcp"), 
  fwerRate=0.25,
  paraJobs=NULL,
  bootB=500,
  bootLassoAlpha=0.05,
  standardize=F,
  sequentialRun=F,
  allFunc=allUserFunc(),
  refReadsThresh=0.2,
  SDThresh=0.05,
  SDquantilThresh=0,
  balanceCut=0.2,
  seed=1
){
  results=list()
  start.time = proc.time()[3] 
  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=testCov,
                   ctrlCov=ctrlCov,testMany=testMany,ctrlMany=ctrlMany)
  data=runMeta$data
  results$covariatesData=runMeta$covariatesData
  binaryInd=runMeta$binaryInd
  covsPrefix=runMeta$covsPrefix
  Mprefix=runMeta$Mprefix
  testCovInd=runMeta$testCovInd
  testCovInOrder=runMeta$testCovInOrder
  testCovInNewNam=runMeta$testCovInNewNam
  ctrlCov=runMeta$ctrlCov
  microbName=runMeta$microbName
  results$covriateNames=runMeta$xNames
  rm(runMeta)
  
  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa have no sequencing reads 
      or are not in the data set. Double check the names of the reference taxa and their 
           sparsity levels.")
    }
  }
  
  if(nRef>(length(microbName))){
    stop("Error: number of random reference taxa can not be larger than the total number
           of taxa in the data. Try lower nRef")
  }
  results$analysisResults=Regulariz(data=data,testCovInd=testCovInd,
                                    testCovInOrder=testCovInOrder,
                                    testCovInNewNam=testCovInNewNam,
                                    microbName=microbName,nRef=nRef,
                                    nRefMaxForEsti=nRefMaxForEsti,
                                    nPermu=nPermu,binaryInd=binaryInd,
                                    x1permut=x1permut,
                                    covsPrefix=covsPrefix,Mprefix=Mprefix,
                                    refTaxa=refTaxa,paraJobs=paraJobs,
                                    reguMethod=reguMethod,fwerRate=fwerRate,
                                    bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                    standardize=standardize,
                                    sequentialRun=sequentialRun,
                                    allFunc=allFunc,refReadsThresh=refReadsThresh,
                                    SDThresh=SDThresh,
                                    SDquantilThresh=SDquantilThresh,
                                    balanceCut=balanceCut,seed=seed
  )
  rm(data)
  
  results$testCov=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  results$bootB=bootB
  results$bootLassoAlpha=bootLassoAlpha
  results$refReadsThresh=refReadsThresh
  results$balanceCut=balanceCut
  results$SDThresh=SDThresh
  results$paraJobs=paraJobs
  results$SDquantilThresh=SDquantilThresh
  results$nRef=nRef
  results$nPermu=nPermu
  results$x1permut=x1permut
  
  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }
  
  rm(testCovInOrder,ctrlCov,microbName)
  
  totalTimeMins = (proc.time()[3] - start.time)/60
  cat("The entire analysis took",totalTimeMins, "minutes","\n")
  
  results$totalTimeMins=totalTimeMins
  
  return(results)
}


##' @export


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

##' Conditional regression for microbiome analysis based on multivariate zero-inflated logistic normal model
##'
##' Make inference on the associations of microbiome with covariates given a user-specified reference taxon/OTU/ASV.
##'
##' The `MZILN()` function can implement the Multivariate Zero-Inflated Logistic Normal model. It estimate and test the association given a user-specified reference taxon/OTU/ASV, whereas the `IFAA()` does not require any user-specified reference taxa.
##'
##' @param MicrobData Microbiome data matrix containing microbiome abundance with each row per sample and each column per taxon/OTU/ASV. It should contain an `"id"` variable to correspond to the `"id"` variable in the covariates data: `CovData`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row per sample and each column per variable. It should also contain an `"id"` variable to correspond to the `"id"` variable in the microbiome data: `MicrobData`.
##' @param linkIDname Variable name of the `"id"` variable in both `MicrobData` and `CovData`. The two data sets will be merged by this `"id"` variable.
##' @param allCov All covariates of interest (including confounders) for estimating and testing their associations with microbiome. Default is all covariates in covData are of interest.
##' @param refTaxa Reference taxa specified by the user and will be used as the reference taxa. 
##' @param reguMethod regularization approach used in phase 1 of the algorithm. Default is `"mcp"`. Other methods are under development.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Sometimes parallel jobs can not be successfully run for unknown reasons. For example, socket related errors may pop up or some slave cores return simple error instead of numerical results. In those scenarios, setting `sequentialRun = TRUE` may help, but it will take more time to run. Default is `TRUE`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. Default is `8`. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, all design matrix X in phase 1 and phase 2 will be standardized in the analyses. Default is `FALSE`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2. The default is `500`.
##' @param bootLassoAlpha The significance level in phase 2. Default is `0.05`.
##' @param seed Random seed for reproducibility. Default is `1`.
##' @return A list containing the estimation results.
##'
##' - `analysisResults$estByRefTaxaList`: A list containing estimating results for all reference taxa and all the variables in 'allCov'. See details.
##'
##' - `covariatesData`: A dataset containing all covariates used in the analyses.
##'
##' @examples
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' results <- MZILN(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 allCov=c("v1","v2","v3"),
##'                 refTaxa=c("rawCount11"))
##'
##' @references Li et al.(2018) Conditional Regression Based on a Multivariate Zero-Inflated Logistic-Normal Model for Microbiome Relative Abundance Data. Statistics in Biosciences 10(3): 587-608 
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.

##' @export
##' @md




MZILN=function(
  MicrobData,
  CovData,
  linkIDname,
  allCov=NULL,
  refTaxa,
  reguMethod=c("mcp"), 
  paraJobs=NULL,
  bootB=500,
  bootLassoAlpha=0.05,
  standardize=F,
  sequentialRun=T,
  allFunc=allUserFunc(),
  seed=1
){
  results=list()
  start.time = proc.time()[3] 
  MZILN=T
  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=allCov,MZILN=MZILN)
  
  data=runMeta$data
  results$covariatesData=runMeta$covariatesData
  binaryInd=runMeta$binaryInd
  covsPrefix=runMeta$covsPrefix
  Mprefix=runMeta$Mprefix
  testCovInd=runMeta$testCovInd
  testCovInOrder=runMeta$testCovInOrder
  testCovInNewNam=runMeta$testCovInNewNam
  ctrlCov=runMeta$ctrlCov
  microbName=runMeta$microbName
  results$covriateNames=runMeta$xNames
  rm(runMeta)
  
  nRef=length(refTaxa)
  
  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa have no sequencing reads 
      or are not in the data set. Double check the names of the reference taxa and their 
           sparsity levels.")
    }
  }
  
  results$analysisResults=Regulariz_MZILN(data=data,testCovInd=testCovInd,
                                          testCovInOrder=testCovInOrder,microbName=microbName,
                                          binaryInd=binaryInd,covsPrefix=covsPrefix,Mprefix=Mprefix,
                                          refTaxa=refTaxa,paraJobs=paraJobs,reguMethod=reguMethod,
                                          bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                          standardize=standardize,sequentialRun=sequentialRun,
                                          allFunc=allFunc,seed=seed
  )
  rm(data)
  
  results$testCov=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  results$bootB=bootB
  results$bootLassoAlpha=bootLassoAlpha
  results$paraJobs=paraJobs
  results$nRef=nRef
  
  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }
  
  rm(testCovInOrder,ctrlCov,microbName)
  
  totalTimeMins = (proc.time()[3] - start.time)/60
  cat("The entire analysis took",totalTimeMins, "minutes","\n")
  
  results$totalTimeMins=totalTimeMins
  
  return(results)
}

##' @export

originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
  standardize,
  sequentialRun,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  
  startT=proc.time()[3]
  cat("start Original screen","\n")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  
  c1<-snow::makeCluster(paraJobs)
  
  if(!sequentialRun){
    cat(paraJobs, "parallel jobs are registered for analyzing", nRef, "reference taxa in Phase 1a.","\n")
  }
  
  snow::clusterExport(c1, allFunc)
  doSNOW::registerDoSNOW(c1)
  
  if(sequentialRun){foreach::registerDoSEQ()}
  
  startT1=proc.time()[3]
  cat("OriginDataScreen parallel setup took",startT1-startT,"seconds","\n")
  # start parallel computing
  scr1Resu=foreach(i=1:nRef,.multicombine=T,
                   .packages=c("picasso","expm","doSNOW","snow","foreach","Matrix"),
                   .errorhandling="pass") %dopar% {
                     
                     ii=which(taxaNames==refTaxa[i])
                     dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                               covsPrefix=covsPrefix)
                     xTildLongTild.i=dataForEst$xTildalong
                     yTildLongTild.i=dataForEst$UtildaLong
                     rm(dataForEst)
                     gc()
                     
                     if(method=="mcp") {
                       Penal.i=runPicasso(x=xTildLongTild.i,y=yTildLongTild.i,
                                          lambda=lambda,nPredics=nPredics,
                                          method="mcp",permutY=F,
                                          standardize=standardize,
                                          seed=seed,seedi=i)
                     }
                     rm(xTildLongTild.i)
                     
                     BetaNoInt.i=as(Penal.i$betaNoInt,"sparseVector")
                     rm(Penal.i)
                     gc()
                     selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                     if (ii==1){
                       selection.i[-seq(1,nPredics)]=as(BetaNoInt.i!=0,"sparseVector")
                     }
                     if (ii==nTaxa) {
                       selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=as(BetaNoInt.i!=0,"sparseVector")
                     }
                     if ((ii>1) & (ii<nTaxa)) {
                       selection.i[1:(nPredics*(ii-1))]=as(BetaNoInt.i[1:(nPredics*(ii-1))]!=0,"sparseVector")
                       selection.i[(nPredics*ii+1):nAlphaSelec]=as(BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]!=0,"sparseVector")
                     }
                     rm(BetaNoInt.i)
                     time4=proc.time()[3]
                     # create return vector
                     recturnlist=list()
                     recturnlist[[1]]=selection.i
                     recturnlist[[2]]=yTildLongTild.i
                     
                     rm(selection.i,yTildLongTild.i)
                     return(recturnlist)
                   }
  snow::stopCluster(c1)
  rm(data)
  gc()
  
  endT=proc.time()[3]
  
  cat("Original screen done and took",(endT-startT1)/60,"minutes","\n")
  
  selecList=list()
  for(i in 1:nRef){
    selecList[[i]]=scr1Resu[[i]][[1]]
  }
  
  results$yTildLongList=list()
  for(i in 1:nRef){
    results$yTildLongList[[i]]=scr1Resu[[i]][[2]]
  }
  rm(scr1Resu)
  
  selecList<- lapply(selecList, as, "sparseMatrix")
  scr1ResuSelec=do.call(cbind, selecList)
  rm(selecList)
  
  results$scr1ResuSelec=scr1ResuSelec
  
  # create count of selection for individual testCov
  countOfSelecForAllPred=as(matrix(Matrix::rowSums(scr1ResuSelec),nrow=nPredics),"sparseMatrix")
  testCovCountMat=countOfSelecForAllPred[testCovInd,,drop=F]
  rm(testCovInd,countOfSelecForAllPred)
  
  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred=as(matrix(rep(0,nTaxa),nrow=1),"sparseMatrix")
  for (tax in 1:nTaxa){
    countMatForTaxni=scr1ResuSelec[(1+(tax-1)*nPredics):(tax*nPredics),,drop=F]
    totCountVecForTaxoni=Matrix::colSums(countMatForTaxni)
    countOfSelecForAPred[1,tax]=sum(totCountVecForTaxoni>0)
  }
  rm(tax,scr1ResuSelec,countMatForTaxni,totCountVecForTaxoni)
  gc()
  
  colnames(countOfSelecForAPred)=taxaNames
  rm(taxaNames)
  
  # return results
  results$testCovCountMat=testCovCountMat
  rm(testCovCountMat)
  results$countOfSelecForAPred=countOfSelecForAPred
  rm(countOfSelecForAPred)
  return(results)
}




##' @export




Regulariz=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  microbName,
  nRef,
  nRefMaxForEsti,
  nPermu,
  x1permut,
  refTaxa,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  reguMethod,
  fwerRate,
  bootB,
  bootLassoAlpha,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  seed
){
  results=list()
  regul.start.time = proc.time()[3]
  
  nTestCov=length(testCovInd)
  
  dataSparsCheck(data=data,Mprefix=Mprefix)
  
  # load abundance data info
  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)
  
  results$x1permut=x1permut
  
  regul.start.time = proc.time()[3]
  cat("Start Phase 1 association identification","\n")
  selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                           testCovInOrder=testCovInOrder,
                           testCovInNewNam=testCovInNewNam,nRef=nRef,
                           nPermu=nPermu,paraJobs=paraJobs,
                           x1permut=x1permut,
                           refTaxa=NULL,fwerRate=fwerRate,
                           method=reguMethod,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           allFunc=allFunc,
                           refReadsThresh=refReadsThresh,
                           SDThresh=SDThresh,
                           SDquantilThresh=SDquantilThresh,
                           balanceCut=balanceCut,
                           Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binaryInd,seed=seed)
  
  results$selecTaxaFWER=selectRegroup$selecTaxaFWER
  results$selecCountOverall=selectRegroup$selecCountOverall
  results$selecMatIndv=selectRegroup$selecMatIndv
  results$selecIndvInOverall=selectRegroup$selecIndvInOverall
  results$selecCountMatIndv=selectRegroup$selecCountMatIndv
  results$maxVec=selectRegroup$maxVec
  results$MaxMatTestCovByPermu=selectRegroup$MaxMatTestCovByPermu
  
  finalIndpRefTax=microbName[taxaNames%in%(selectRegroup$finalIndpRefTax)]
  results$finalRefTaxonQualified=selectRegroup$refTaxonQualified
  
  results$goodIndpRefTaxLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxLeastCount)]
  results$goodIndpRefTaxWithCount=selectRegroup$goodIndpRefTaxWithCount
  names(results$goodIndpRefTaxWithCount)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxWithCount)]
  
  results$goodIndpRefTaxFWERcutLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxFWERcutLeastCount)]
  results$goodIndpRefTaxFWERcut=selectRegroup$goodIndpRefTaxFWERcut
  names(results$goodIndpRefTaxFWERcut)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxFWERcut)]
  
  results$goodRefTaxaCandi=microbName[taxaNames%in%(selectRegroup$goodRefTaxaCandi)]
  
  results$taxaLessGoodCut=selectRegroup$taxaLessGoodCut
  names(results$taxaLessGoodCut)=microbName[taxaNames%in%names(selectRegroup$taxaLessGoodCut)]
  
  
  results$taxaLessFWERCut=selectRegroup$taxaLessFWERCut
  names(results$taxaLessFWERCut)=microbName[taxaNames%in%names(selectRegroup$taxaLessFWERCut)]
  
  results$randomRefTaxa=microbName[taxaNames%in%(selectRegroup$refTaxa)]
  
  results$fwerCut=selectRegroup$fwerCut
  results$twoMeanUsed=selectRegroup$twoMeanUsed
  results$fwerCutIndv=selectRegroup$fwerCutIndv
  rm(selectRegroup)
  
  MCPExecuTime = (proc.time()[3] - regul.start.time)/60
  results$MCPExecuTime=MCPExecuTime
  cat("Phase 1 Associaiton identification is done and used", MCPExecuTime,"minutes","\n")
  
  results$finalizedBootRefTaxon=finalIndpRefTax
  
  startT=proc.time()[3]
  cat("Start Phase 2 parameter estimation","\n")
  
  unestimableTaxa=c()
  qualifyData=data
  if(length(binaryInd)>0){
    qualifyData=data[rowSums(data[,results$taxaNames]>0)>=2,,drop=F]
    firstBinPredNam=paste0(covsPrefix,binaryInd)
    allBinPred=paste0(covsPrefix,binaryInd:nPredics)
    nBinPred=length(allBinPred)
    
    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax=taxaNames[!taxaNames%in%(results$finalizedBootRefTaxon)]
    unbalanceTaxa=c()
    for(i in AllTaxaNamesNoRefTax){
      for(j in allBinPred){
        twoColumns.ij=qualifyData[,c(i,j)]
        nNonZero=sum(twoColumns.ij[,1]>0)
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        if(sumOfBin%in%c(0,1,(nNonZero-1),nNonZero)){
          unbalanceTaxa=c(unbalanceTaxa,i)
        }
      }
    }
    unestimableTaxa=unique(unbalanceTaxa)
    rm(allBinPred,unbalanceTaxa)
  }
  
  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(qualifyData[,taxaNames])==0)
  rm(qualifyData)
  unestimableTaxa=unique(c(unestimableTaxa,taxaNames[TaxaNoReads]))
  results$unEstTaxaPos=which(taxaNames%in%unestimableTaxa)
  rm(TaxaNoReads,unestimableTaxa)
  
  goodIndpRefTax.ascend=sort(results$goodIndpRefTaxWithCount)
  #print("goodIndpRefTax.ascend:")
  #print(goodIndpRefTax.ascend)
  
  goodIndpRefTaxNam=names(goodIndpRefTax.ascend)
  allRefTaxNam=unique(c(results$finalizedBootRefTaxon,goodIndpRefTaxNam))
  nGoodIndpRef=length(allRefTaxNam)
  results$allRefTaxNam=allRefTaxNam
  
  #print("allRefTaxNam:")
  #print(allRefTaxNam)
  
  results$nRefUsedForEsti=min(nGoodIndpRef,nRefMaxForEsti)
  
  cat("Final Reference Taxa are:",allRefTaxNam[seq(results$nRefUsedForEsti)],"\n")
  
  #selecResu=paste0("../../../ufrc/selecResul_",sample,".RData")
  #selecResu=paste0("./sim_results/selecResul_",sample,".RData")
  
  #save(results,file=selecResu)
  
  results$estiList=list()
  for(iii in 1:(results$nRefUsedForEsti)){
    cat("Start estimation for the", iii,"th final reference taxon:",allRefTaxNam[iii],"\n")
    time11=proc.time()[3]
    newRefTaxNam=taxaNames[microbName%in%(allRefTaxNam[iii])]
    results$estiList[[allRefTaxNam[iii]]]=bootResuHDCI(data=data,
                                                       refTaxa=newRefTaxNam,
                                                       bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                                       binPredInd=binaryInd,covsPrefix=covsPrefix,
                                                       Mprefix=Mprefix,
                                                       standardize=standardize,
                                                       seed=seed)
    time12=proc.time()[3]
    cat("Estimation done for the", iii,"th final reference taxon:",allRefTaxNam[iii],
        "and it took",(time12-time11)/60,"minutes","\n")
  }
  estiResults=results$estiList[[results$finalizedBootRefTaxon]]
  
  endT=proc.time()[3]
  
  cat("Phase 2 parameter estimation done and took",(endT-startT)/60,"minutes.","\n")
  
  
  results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
  results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
  results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")
  
  results$betaMat.LPR=as(matrix(estiResults$finalBetaEst.LPR,nrow=nPredics),"sparseMatrix")
  results$CILowMat.LPR=as(matrix(estiResults$CIvecLow.LPR,nrow=nPredics),"sparseMatrix")
  results$CIUpMat.LPR=as(matrix(estiResults$CIvecUp.LPR,nrow=nPredics),"sparseMatrix")
  rm(estiResults)
  
  estByCovList=list()
  
  for(i in 1:nTestCov){
    sigTaxaPosition=which(results$selecIndvInOverall[i,]!=0)
    nrow=length(sigTaxaPosition)
    if(nrow==0)next
    
    ncol=3
    estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)
    
    for(j in 1:nrow){
      for(k in 1:ncol){
        if(k==1)estByCovMat[j,k]=results$betaMat.LPR[i,sigTaxaPosition[j]]
        if(k==2)estByCovMat[j,k]=results$CILowMat.LPR[i,sigTaxaPosition[j]]
        if(k==3)estByCovMat[j,k]=results$CIUpMat.LPR[i,sigTaxaPosition[j]]
      }
    }
    
    rownames(estByCovMat)=microbName[results$selecIndvInOverall[i,]!=0]
    colnames(estByCovMat)=c("Beta.LPR","LowB95%CI.LPR","UpB95%CI.LPR")  
    
    estByCovList[[testCovInOrder[i]]]=estByCovMat
    rm(estByCovMat)
  }
  
  if(length(estByCovList)==0){
    results$estByCovList="No significant assoication is identified. You may consider a higher FDR level."
  }else{
    results$estByCovList=estByCovList
  }
  rm(estByCovList)
  
  SigCovByTaxaList=list()
  for(i in 1:nTaxa){
    sigCov=which(results$selecIndvInOverall[,i]!=0)
    if(length(sigCov)==0)next
    SigCovByTaxaList[[microbName[i]]]=testCovInOrder[sigCov]
  }
  if(length(SigCovByTaxaList)==0){
    results$SigCovByTaxaList="No significant assoication is identified. You may consider a higher FDR level."
  }else{
    results$SigCovByTaxaList=SigCovByTaxaList
  }
  rm(SigCovByTaxaList,microbName)
  
  results$reguMethod=reguMethod
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  return(results)
  
  rm(data)
  
  # return results
  results$fwerRate=fwerRate
  results$nPermu=nPermu
  results$nRef=nRef
  return(results)
}



##' @export



Regulariz_MZILN=function(
  data,
  testCovInd,
  testCovInOrder,
  microbName,
  refTaxa,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  reguMethod,
  fwerRate,
  bootB,
  bootLassoAlpha,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  seed
){
  results=list()
  regul.start.time = proc.time()[3]
  
  nTestCov=length(testCovInd)
  
  dataSparsCheck(data=data,Mprefix=Mprefix)
  
  # load abundance data info
  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)
  
  newRefTaxNam=taxaNames[microbName%in%refTaxa]
  nRef=length(refTaxa)
  refTaxa_reOrder=microbName[microbName%in%refTaxa]
  
  reguResu=originDataScreen(method=reguMethod,data=data,testCovInd=testCovInd,
                            paraJobs=paraJobs,lambda=NULL,refTaxa=newRefTaxNam,standardize=standardize,
                            sequentialRun=sequentialRun,allFunc=allFunc,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,binPredInd=binaryInd,seed=seed)
  
  results$SelecAllRefTaxaPred=reguResu$scr1ResuSelec
  rm(reguResu)
  
  cat("Reference taxa are:",refTaxa,"\n")
  
  betaMatList=list()
  CILowMatList=list()
  CIUpMatList=list()
  results$estiList=list()
  for(iii in 1:nRef){
    time11=proc.time()[3]
    results$estiList[[refTaxa_reOrder[iii]]]=bootResuHDCI(data=data,
                                                          refTaxa=newRefTaxNam[iii],bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                                          binPredInd=binaryInd,covsPrefix=covsPrefix,Mprefix=Mprefix,
                                                          standardize=standardize,seed=seed)
    
    time12=proc.time()[3]
    cat("Estimation done for the", iii,"th reference taxon:",refTaxa[iii],
        "and it took",(time12-time11)/60,"minutes","\n")
    
    estiResults=results$estiList[[refTaxa_reOrder[iii]]]
    
    betaMatList[[iii]]=as(matrix(estiResults$finalBetaEst.LPR,nrow=nPredics),"sparseMatrix")
    CILowMatList[[iii]]=as(matrix(estiResults$CIvecLow.LPR,nrow=nPredics),"sparseMatrix")
    CIUpMatList[[iii]]=as(matrix(estiResults$CIvecUp.LPR,nrow=nPredics),"sparseMatrix")
  }
  rm(estiResults)
  
  estByRefTaxaList=list()
  estByCovList=list()
  
  for(iii in 1:nRef){
    SelecRefTaxa.i=matrix(results$SelecAllRefTaxaPred[,iii],nrow=nPredics)
    for(i in 1:nTestCov){
      sigTaxaPosition=which(SelecRefTaxa.i[i,]!=0)
      nrow=length(sigTaxaPosition)
      if(nrow==0)next
      
      ncol=3
      estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)
      
      for(j in 1:nrow){
        for(k in 1:ncol){
          if(k==1)estByCovMat[j,k]=betaMatList[[iii]][i,sigTaxaPosition[j]]
          if(k==2)estByCovMat[j,k]=CILowMatList[[iii]][i,sigTaxaPosition[j]]
          if(k==3)estByCovMat[j,k]=CIUpMatList[[iii]][i,sigTaxaPosition[j]]
        }
      }
      
      rownames(estByCovMat)=microbName[SelecRefTaxa.i[i,]!=0]
      colnames(estByCovMat)=c("Beta.LPR","LowB95%CI.LPR","UpB95%CI.LPR")  
      
      estByCovList[[testCovInOrder[i]]]=estByCovMat
      rm(estByCovMat)
    }
    
    if(length(estByCovList)==0){
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["estByCovList"]]="No significant assoication is identified."
    }else{
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["estByCovList"]]=estByCovList
    }
  }
  
  
  rm(estByCovList)
  
  SigCovByTaxaList=list()
  
  for(iii in 1:nRef){
    SelecRefTaxa.i=matrix(results$SelecAllRefTaxaPred[,iii],nrow=nPredics)
    for(i in 1:nTaxa){
      sigCov=which(SelecRefTaxa.i[,i]!=0)
      if(length(sigCov)==0)next
      SigCovByTaxaList[[microbName[i]]]=testCovInOrder[sigCov]
    }
    if(length(SigCovByTaxaList)==0){
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["SigCovByTaxaList"]]="No significant assoication is identified."
    }else{
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["SigCovByTaxaList"]]=SigCovByTaxaList
    }
  }
  
  rm(SigCovByTaxaList,microbName)
  
  results$reguMethod=reguMethod
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  rm(data)
  # return results
  return(results)
}


##' @export



runBootLassoHDCI=function(
  x,
  y,
  nPredics,
  nTaxa,
  nfolds=10,
  lambdaOPT=NULL,
  refTaxaPosition,
  #zeroSDCut=10^(-6),
  zeroSDCut=0,
  correCut=1.1,
  standardize,
  bootB,
  bootLassoAlpha,
  seed
){
  
  results=list()
  
  # check if there are x variables with zero variance
  nObsAll=length(y)
  nBeta=ncol(x)
  
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  
  # write.csv(cbind(as.matrix(x[,-xWithNearZeroSd]),as.vector(y)),file="xy.csv",row.names = F)
  
  df.cor=suppressWarnings(cor(as.matrix(x)))
  df.cor[is.na(df.cor)]=0
  df.cor[!lower.tri(df.cor)]=0
  excluCorColumns=which(apply(df.cor, 2, function(x) any(abs(x)>=correCut)))
  
  xWithNearZeroSd=sort(unique(c(xWithNearZeroSd,excluCorColumns)))
  rm(excluCorColumns)
  
  # remove near constant columns
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd]
  }
  rm(sdX)
  
  nearZeroSd=length(xWithNearZeroSd)
  
  # print(paste("Number of x with near-zero sd: ", nearZeroSd))
  
  cvStartTime= proc.time()[3]
  
  availCores=availableCores()
  if(is.numeric(availCores)){
    ncores.boot=max(1,availableCores())
    #ncores.boot=max(1,availableCores()-1)
  }else{
    ncores.boot=1
  }
  
  cat(ncores.boot, "parallel jobs are registered for bootstrp in Phase 2.","\n")
  
  c3 <- parallel::makeCluster(ncores.boot)
  doParallel::registerDoParallel(c3)
  
  if(length(seed)>0){
    set.seed(as.numeric(seed)+10^2)
    parallel::clusterSetRNGStream(cl=c3,(as.numeric(seed)+10^3))
  }
  
  bootResu=bootLOPR(x=x,y=as.vector(y),B=bootB,nfolds=nfolds,
                    standardize=standardize,parallel.boot=T,
                    ncores.boot=ncores.boot,alpha=bootLassoAlpha)
  parallel::stopCluster(c3)
  
  # convert to a sparse vector format from sparse matrix format
  beta=as(bootResu$Beta,"sparseVector")
  betaCI=as(bootResu$interval,"sparseMatrix")
  
  beta.LPR=as(bootResu$Beta.LPR,"sparseVector")
  betaCI.LPR=as(bootResu$interval.LPR,"sparseMatrix")
  rm(bootResu)
  # transform vector back
  if(length(xWithNearZeroSd)>0){
    betaTransLasso=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                       unSelectList=sort(xWithNearZeroSd),newBetaNoInt=beta)
    beta=betaTransLasso$finalBeta
    
    betaTransLasso.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                         unSelectList=sort(xWithNearZeroSd),newBetaNoInt=beta.LPR)
    beta.LPR=betaTransLasso.L$finalBeta
    
    betaTransCIlow=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                       unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI[1,])
    betaCIlow=betaTransCIlow$finalBeta
    
    betaTransCIlow.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                         unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI.LPR[1,])
    betaCIlow.LPR=betaTransCIlow.L$finalBeta
    
    betaTransCIhi=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                      unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI[2,])
    betaCIhi=betaTransCIhi$finalBeta
    
    betaTransCIhi.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                        unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI.LPR[2,])
    betaCIhi.LPR=betaTransCIhi.L$finalBeta
    
    rm(betaTransLasso,betaTransLasso.L,betaTransCIlow,
       betaTransCIlow.L,betaTransCIhi,betaTransCIhi.L)
  } else{
    betaCIlow=betaCI[1,]
    betaCIhi=betaCI[2,]
    betaCIlow.LPR=betaCI.LPR[1,]
    betaCIhi.LPR=betaCI.LPR[2,]
  }
  
  rm(xWithNearZeroSd)
  
  BetaNoInt.i=as(beta[-seq(1,length(beta),by=(nPredics+1))],"sparseVector")
  betaCIlow.i=as(betaCIlow[-seq(1,length(betaCIlow),by=(nPredics+1))],"sparseVector")
  betaCIhi.i=as(betaCIhi[-seq(1,length(betaCIhi),by=(nPredics+1))],"sparseVector")
  
  BetaNoInt.LPR.i=as(beta.LPR[-seq(1,length(beta.LPR),by=(nPredics+1))],"sparseVector")
  betaCIlow.LPR.i=as(betaCIlow.LPR[-seq(1,length(betaCIlow.LPR),by=(nPredics+1))],"sparseVector")
  betaCIhi.LPR.i=as(betaCIhi.LPR[-seq(1,length(betaCIhi.LPR),by=(nPredics+1))],"sparseVector")
  
  nAlphaSelec=nPredics*nTaxa
  nAlphaNoInt=nPredics*(nTaxa-1)
  
  BetaNoIntWithRef=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIlow=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIhi=as(rep(0,nAlphaSelec),"sparseVector")
  
  BetaNoIntWithRef.LPR=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIlow.LPR=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIhi.LPR=as(rep(0,nAlphaSelec),"sparseVector")
  
  ii=refTaxaPosition
  
  if (ii==1){
    BetaNoIntWithRef[-seq(1,nPredics)]=BetaNoInt.i
    betaCIlow[-seq(1,nPredics)]=betaCIlow.i
    betaCIhi[-seq(1,nPredics)]=betaCIhi.i
    
    BetaNoIntWithRef.LPR[-seq(1,nPredics)]=BetaNoInt.LPR.i
    betaCIlow.LPR[-seq(1,nPredics)]=betaCIlow.LPR.i
    betaCIhi.LPR[-seq(1,nPredics)]=betaCIhi.LPR.i
  }
  if (ii==nTaxa){
    BetaNoIntWithRef[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
    betaCIlow[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIlow.i
    betaCIhi[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIhi.i
    
    BetaNoIntWithRef.LPR[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.LPR.i
    betaCIlow.LPR[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIlow.LPR.i
    betaCIhi.LPR[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIhi.LPR.i
  }
  if ((ii>1) & (ii<nTaxa)) {
    BetaNoIntWithRef[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
    BetaNoIntWithRef[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
    betaCIlow[1:(nPredics*(ii-1))]=betaCIlow.i[1:(nPredics*(ii-1))]
    betaCIlow[(nPredics*ii+1):nAlphaSelec]=betaCIlow.i[(nPredics*(ii-1)+1):nAlphaNoInt]
    betaCIhi[1:(nPredics*(ii-1))]=betaCIhi.i[1:(nPredics*(ii-1))]
    betaCIhi[(nPredics*ii+1):nAlphaSelec]=betaCIhi.i[(nPredics*(ii-1)+1):nAlphaNoInt]
    
    BetaNoIntWithRef.LPR[1:(nPredics*(ii-1))]=BetaNoInt.LPR.i[1:(nPredics*(ii-1))]
    BetaNoIntWithRef.LPR[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.LPR.i[(nPredics*(ii-1)+1):nAlphaNoInt]
    betaCIlow.LPR[1:(nPredics*(ii-1))]=betaCIlow.LPR.i[1:(nPredics*(ii-1))]
    betaCIlow.LPR[(nPredics*ii+1):nAlphaSelec]=betaCIlow.LPR.i[(nPredics*(ii-1)+1):nAlphaNoInt]
    betaCIhi.LPR[1:(nPredics*(ii-1))]=betaCIhi.LPR.i[1:(nPredics*(ii-1))]
    betaCIhi.LPR[(nPredics*ii+1):nAlphaSelec]=betaCIhi.LPR.i[(nPredics*(ii-1)+1):nAlphaNoInt]
  }
  
  rm(BetaNoInt.i,betaCIlow.i,betaCIhi.i,
     BetaNoInt.LPR.i,betaCIlow.LPR.i,betaCIhi.LPR.i)
  
  results$beta=BetaNoIntWithRef
  results$betaCIlow=betaCIlow
  results$betaCIhi=betaCIhi
  
  results$beta.LPR=BetaNoIntWithRef.LPR
  results$betaCIlow.LPR=betaCIlow.LPR
  results$betaCIhi.LPR=betaCIhi.LPR
  
  rm(BetaNoIntWithRef,betaCIlow,betaCIhi,
     BetaNoIntWithRef.LPR,betaCIlow.LPR,betaCIhi.LPR)
  
  return(results)
}


##' @export



runPicasso=function(
  x,
  y,
  nPredics,
  family="gaussian",
  nfolds=10,
  nCV=1,
  cvRestartPosition=0.9,
  zeroSDCut=0,
  lambda=NULL,
  lambda.min.ratio=0.05,
  nLam=100,
  method="mcp",
  standardize,
  permutY,
  allFunc=allFunc,
  seed,
  seedi
){
  
  results=list()
  
  # check if there are x variables with zero variance
  nObsAll=length(y)
  nBeta=ncol(x)
  
  if(permutY){
    # permute y
    y=sample(y,size=nObsAll)
  }
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  
  # remove near constant columns
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd]
  }
  rm(sdX)
  nearZeroSd=length(xWithNearZeroSd)
  
  # calculate lambda max
  lamMax=max(abs(Matrix::colSums(as.matrix(x)*as.vector(y))))/nObsAll
  
  if(length(lambda)!=1){
    for(i in 1:nCV){
      cvResul=cvPicasso(x=x,y=y,lamMax=lamMax,lambda=lambda,nfolds=nfolds,zeroSDCut=zeroSDCut,
                        lambda.min.ratio=lambda.min.ratio,nLam=nLam,
                        method=method,family=family,standardize=standardize,seed=seed,seedi=seedi)
    }
    
    results$lamList=cvResul$lamList
    results$SSE=cvResul$SSE
    lambda=cvResul$optiLam
  }
  
  rm(cvResul)
  
  finalRun=picasso(X=x,Y=y,lambda=lambda,method=method,family=family,standardize=standardize)
  endT=proc.time()[3]
  rm(x,y)
  
  # convert to a sparse vector format from sparse matrix format
  
  OverallIntercp=finalRun$intercept[1]
  finalRunBeta=as(finalRun$beta[,1],"sparseVector")
  if(length(xWithNearZeroSd)>0){
    betaTrans=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                  unSelectList=sort(xWithNearZeroSd),newBetaNoInt=finalRunBeta)
    beta=betaTrans$finalBeta
    rm(betaTrans)
  } else {
    beta=finalRunBeta
  }
  lambda=finalRun$lambda
  
  rm(finalRun,finalRunBeta)
  
  rm(xWithNearZeroSd)
  
  # return results
  results$betaInt=beta
  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  rm(beta,nPredics)
  results$lambda=lambda
  results$overalIntercp=OverallIntercp
  rm(lambda)
  
  results$successCV=1
  results$nCV=nCV
  rm(nCV)
  return(results)
}


##' @export


runScrParal=function(
  method=c("mcp"),
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  doPermut,
  permutY=F,
  x1permut,
  nPermu,
  refTaxa=NULL,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,qualifyRefTax=T,
                     refReadsThresh=refReadsThresh,
                     SDThresh=SDThresh,SDquantilThresh=SDquantilThresh,
                     balanceCut=balanceCut,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  
  results$goodRefTaxaCandi=basicInfo$goodRefTaxaCandi
  nSub=basicInfo$nSub
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  # make reference taxa list
  if(length(refTaxa)<nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa_extra=sample((taxaNames[!(taxaNames%in%refTaxa)]),(nRef-length(refTaxa)))
    refTaxa=c(refTaxa,refTaxa_extra)
    results$refTaxa=refTaxa
  }
  
  if(length(refTaxa)>=nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(refTaxa,nRef)
    results$refTaxa=refTaxa
  }
  #
  ## run original data screen
  #
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           seed=seed)
  
  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  rm(screen1)
  gc()
  
  if(!x1permut){
    # generate Xbeta and residues for reduced model permutation
    XbetResi=XbetaAndResidu(data=data,testCovInd=testCovInd,
                            nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                            method=method,allFunc=allFunc,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,standardize=standardize,
                            binPredInd=binPredInd,
                            sequentialRun=sequentialRun,
                            seed=seed)
    xBetaList=XbetResi$xBetaList
    residuList=XbetResi$residuList
    xTildLong=XbetResi$xTildLong
    rm(XbetResi)
  }
  #
  ## start to run permutation screen
  #
  if(doPermut){
    startT=proc.time()[3]
    cat("start to run permutation","\n")
    
    # permut the exposure variable
    if(length(seed)>0)set.seed(as.numeric(seed)+10^6)
    
    permutOrder=lapply(rep(nSub,nPermu),sample)
    if(!x1permut){
      residuPermuOrder=lapply(rep(length(residuList[[1]]),nPermu),sample)
    }
    
    screenStartTime = proc.time()[3]
    
    EName=testCovInNewNam
    EVar=data[,EName,drop=F]
    
    totNumOfLoops=nRef*nPermu
    
    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
      if(!is.numeric(availCores))paraJobs=1
    }
    
    c2 <- snow::makeCluster(paraJobs)
    
    if(!sequentialRun){
      cat(paraJobs, "parallel jobs are registered for the permutation analysis in Phase 1b.","\n")
    }
    
    snow::clusterExport(c2, allFunc)
    doSNOW::registerDoSNOW(c2)
    
    if(sequentialRun){foreach::registerDoSEQ()}
    
    refResu=foreach (i=1:totNumOfLoops,.multicombine=T,
                     .packages=c("picasso","expm","doSNOW","snow","foreach","Matrix"),
                     .errorhandling="pass") %dopar% {
                       #for(j in 1:nRef){
                       
                       permut.i=1+(i-1)%%nPermu
                       ref.i=1+floor((i-1)/nPermu)
                       ii=which(taxaNames==refTaxa[ref.i])
                       
                       if(x1permut){
                         permutX1=EVar[permutOrder[[permut.i]],,drop=F]
                         newData=data
                         newData[,EName]=permutX1
                         rm(permutX1)
                         
                         xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=T)
                         xLongTild.i=xLong.i$xTildalong
                         rm(newData,xLong.i)
                         gc()
                         
                         if(method=="mcp") {
                           Penal.i=runPicasso(x=xLongTild.i,y=yTildLongList[[ref.i]],nPredics=nPredics,
                                              method="mcp",permutY=permutY,
                                              standardize=standardize,
                                              seed=seed,seedi=i)
                         }
                         rm(xLongTild.i)
                       }else{
                         resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                         yTildLong.i=xBetaList[[ref.i]]+resid
                         rm(resid)
                         
                         if(method=="mcp") {
                           Penal.i=runPicasso(x=xTildLong,y=yTildLong.i,nPredics=nPredics,
                                              method="mcp",permutY=permutY,
                                              standardize=standardize,
                                              seed=seed,seedi=i)
                         }
                         rm(yTildLong.i)
                       }
                       BetaNoInt.i=as(Penal.i$betaNoInt,"sparseVector")
                       rm(Penal.i)
                       gc()
                       selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                       if (ii==1){
                         selection.i[-seq(1,nPredics)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if (ii==nTaxa) {
                         selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if ((ii>1) & (ii<nTaxa)) {
                         selection.i[1:(nPredics*(ii-1))]=as(BetaNoInt.i[1:(nPredics*(ii-1))]!=0,"sparseVector")
                         selection.i[(nPredics*ii+1):nAlphaSelec]=as(BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]!=0,"sparseVector")
                       }
                       rm(BetaNoInt.i)
                       # create return vector
                       recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                       recturnVec=selection.i
                       rm(selection.i)
                       return(recturnVec)
                     }
    rm(yTildLongList)
    snow::stopCluster(c2)
    gc()
    
    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)
    
    rm(data,permutOrder,EVar)
    
    endT=proc.time()[3]
    
    cat("Permutation done and took", (endT-startT)/60,"minutes","\n")
    
    # obtain the maximum vector
    permuColInd=1+seq(0,totNumOfLoops-1)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)
    
    for(i in 1:nPermu){
      matrix.i.permu=refResu[,(permuColInd==i),drop=F]
      vec.i=as(rep(0,nTaxa),"sparseVector")
      
      for(j in 1:nRef){
        vector.j.ref=matrix.i.permu[,j]
        matrix.i.j=as(matrix(vector.j.ref,nrow=nPredics),"sparseMatrix")
        testCovVec.j.ref=as((Matrix::colSums(matrix.i.j[testCovInd,,drop=F])>0)+0,"sparseVector")
        vec.i=vec.i+testCovVec.j.ref
      }
      
      permutResuMat[,i]=as.vector(vec.i)
    }
    rm(vec.i,matrix.i.permu,vector.j.ref,matrix.i.j,testCovVec.j.ref)
    
    if(nTestCov>1){
      for(i in 1:nPermu){
        matrix.i.permu=refResu[,(permuColInd==i),drop=F]
        allCovCountMat.i=as(matrix(Matrix::rowSums(matrix.i.permu),nrow=nPredics),"sparseMatrix")
        permutTestCovList[[i]]=allCovCountMat.i[testCovInd,,drop=F]
      }
      rm(matrix.i.permu,allCovCountMat.i)
    }
    
    maxVec=rep(NA,nPermu)
    for(i in 1:nPermu){
      maxVec[i]=max(permutResuMat[,i])
    }
    rm(permutResuMat)
    
    results$nTestCov=nTestCov
    
    if(nTestCov>1){
      MaxMatTestCovByPermu=matrix(NA,nrow=nTestCov,ncol=nPermu)
      for(k in 1:nTestCov){
        for(i in 1:nPermu){
          MaxMatTestCovByPermu[k,i]=max(permutTestCovList[[i]][k,])
        }
      }
      rm(permutTestCovList)
      
      results$MaxMatTestCovByPermu=MaxMatTestCovByPermu
      rm(MaxMatTestCovByPermu)
    }
    
    # obtain the null binomial distribution for each taxa
    totSeleCount=Matrix::rowSums(refResu)
    rm(refResu)
    binomPar=totSeleCount/totNumOfLoops
    rm(totSeleCount,totNumOfLoops)
    
    results$maxVec=maxVec
    rm(maxVec)
    
    results$binomPar=binomPar
    rm(binomPar)
    results$nTaxa=nTaxa
    results$nPredics=nPredics
  }
  
  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}


##' @export



twoMeanClus=function(
  matrix,
  nRef
){
  results=list()
  if(sum(matrix<0)>0) stop("Count matrix contains negative values")
  
  matrix=as.matrix(t(as.matrix(matrix)))
  
  n.col=ncol(matrix)
  n.row=nrow(matrix)
  finalValVec=rep(NA,n.col)
  twoMeanSucVec=rep(NA,n.col)
  finalSelec=matrix(NA,nrow=n.row,ncol=n.col)
  for (i in 1:n.col){
    vecI=matrix[,i]
    cluster2=try(kmeans(vecI,centers=c(0,nRef)),silent=T)
    twoMeanSucVec[i]=1-inherits(cluster2,"try-error")
    
    if(inherits(cluster2,"try-error")){
      finalSelec[1,i]=1
      finalSelec[2:n.row,i]=0
    }
    else{
      hiClus=which(cluster2$centers==max(cluster2$centers))
      finalSelec[,i]=(cluster2$cluster==hiClus)+0
    }
    finalValVec[i]=min(vecI[(finalSelec[,i]==1)])
    rm(vecI)
  }
  
  results$finalSelec=t(as.matrix(finalSelec))
  
  rm(finalSelec)
  results$twoMeanCut=finalValVec
  rm(finalValVec)
  results$twoMeanSucVec=twoMeanSucVec
  rm(twoMeanSucVec)
  return(results)
}

##' @export




groupBetaToFullBeta=function(
  nTaxa,
  nPredics,
  unSelectList,
  newBetaNoInt
){
  
  results=list()
  unSelectList=unique(sort(unSelectList))
  nUnSelec=length(unSelectList)
  nAlphaSelec=nTaxa*nPredics
  nNewBetaNoInt=length(newBetaNoInt)
  
  if(nNewBetaNoInt!=((nTaxa-nUnSelec)*nPredics)) {
    stop("Error: Beta dimension from grouped analyis does not match the expected number")
  }  
  
  if(nTaxa<max(unSelectList) | 1>min(unSelectList)){
    stop("Error: unSelectList out of range")
  }
  finalBeta=newBetaNoInt
  for (i in unSelectList){
    finalBetaTemp=rep(NA,(length(finalBeta)+nPredics))
    lengthTemp=length(finalBetaTemp)
    
    if(i==1){
      finalBetaTemp[1:nPredics]=0
      finalBetaTemp[(nPredics+1):lengthTemp]=finalBeta
      finalBeta=finalBetaTemp
    }
    
    if(i>1){
      if((i*nPredics)<=(length(finalBeta))){
        finalBetaTemp[1:((i-1)*nPredics)]=finalBeta[1:((i-1)*nPredics)]
        finalBetaTemp[((i-1)*nPredics+1):(i*nPredics)]=0
        finalBetaTemp[(i*nPredics+1):lengthTemp]=finalBeta[((i-1)*nPredics+1):(length(finalBeta))]
      }else{
        finalBetaTemp[1:((i-1)*nPredics)]=finalBeta
        finalBetaTemp[((i-1)*nPredics+1):lengthTemp]=0
      }
      finalBeta=finalBetaTemp
    }
  }
  rm(finalBetaTemp)
  results$finalBeta=finalBeta
  rm(finalBeta)
  return(results)
}

##' @export


XbetaAndResidu=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
  sequentialRun,
  standardize,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # generate x matrix for the permutation
  dataForEst1=dataRecovTrans(data=data,ref=refTaxa[1],Mprefix=Mprefix,
                             covsPrefix=covsPrefix)
  results$xTildLong=dataForEst1$xTildalong
  rm(dataForEst1)
  
  # generate reduced data
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  predNames=basicInfo$predNames
  testCovNam=predNames[testCovInd]
  colnames(data)%in%testCovNam
  reducedData=data[,!(colnames(data)%in%testCovNam)]
  rm(basicInfo)
  
  # load reduced data info
  reducedBasicInfo=dataInfo(data=reducedData,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,
                            binPredInd=binPredInd)
  
  taxaNames=reducedBasicInfo$taxaNames
  nTaxa=reducedBasicInfo$nTaxa
  nPredics=reducedBasicInfo$nPredics
  rm(reducedBasicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  
  startT=proc.time()[3]
  cat("start generating residues for permutation","\n")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  
  c4<-snow::makeCluster(paraJobs)
  
  if(!sequentialRun){
    cat(paraJobs, "parallel jobs are registered for generate residues in Phase 1a.","\n")
  }
  
  snow::clusterExport(c4, allFunc)
  doSNOW::registerDoSNOW(c4)
  
  if(sequentialRun){foreach::registerDoSEQ()}
  
  startT1=proc.time()[3]
  # start parallel computing
  residu1Resu=foreach(i=1:nRef,.multicombine=T,
                      .packages=c("picasso","expm","doSNOW","snow","foreach","Matrix"),
                      .errorhandling="pass") %dopar% {
                        
                        ii=which(taxaNames==refTaxa[i])
                        dataForEst=dataRecovTrans(data=reducedData,ref=refTaxa[i],Mprefix=Mprefix,
                                                  covsPrefix=covsPrefix)
                        xTildLongTild.i=dataForEst$xTildalong
                        yTildLongTild.i=dataForEst$UtildaLong
                        rm(dataForEst)
                        gc()
                        
                        if(method=="mcp") {
                          Penal.i=runPicasso(x=xTildLongTild.i,y=yTildLongTild.i,
                                             lambda=lambda,nPredics=nPredics,
                                             standardize=standardize,
                                             method="mcp",permutY=F,seed=seed,seedi=i)
                        }
                        
                        betaInt=Penal.i$betaInt
                        overalIntercp=Penal.i$overalIntercp
                        rm(Penal.i)
                        gc()
                        
                        xBeta=xTildLongTild.i%*%betaInt+overalIntercp
                        rm(xTildLongTild.i)
                        
                        residu=yTildLongTild.i-xBeta
                        rm(yTildLongTild.i)
                        
                        # create return vector
                        recturnlist=list()
                        recturnlist[[1]]=xBeta
                        recturnlist[[2]]=residu
                        return(recturnlist)
                      }
  snow::stopCluster(c4)
  rm(data)
  gc()
  
  endT=proc.time()[3]
  
  cat("Generating residu is done and took",(endT-startT1)/60,"minutes","\n")
  
  xBetaList=list()
  for(i in 1:nRef){
    xBetaList[[i]]=residu1Resu[[i]][[1]]
  }
  
  residuList=list()
  for(i in 1:nRef){
    residuList[[i]]=residu1Resu[[i]][[2]]
  }
  rm(residu1Resu)
  
  # return results
  results$xBetaList=xBetaList
  rm(xBetaList)
  results$residuList=residuList
  rm(residuList)
  return(results)
}

