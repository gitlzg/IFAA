##' @export


originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
  independence=T,
  identity=T,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd){

  results=list()

  # load data info
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics

  nSub=basicInfo$nSub
  nSubNoReads=basicInfo$nSubNoReads
  maxTaxaNameNum=basicInfo$maxTaxaNameNum
  rm(basicInfo)
  gc()

  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa

  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()

  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  cat("start Original screen","\n")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-1)
    if(!is.numeric(availCores))paraJobs=1
  }

  c1<-makeCluster(paraJobs)

  clusterExport(c1, varlist=allFunc)

  clusterSetupRNGstream(cl=c1,seed=1)

  registerDoSNOW(c1)

  # start parallel computing
  scr1Resu=foreach(i=1:nRef,.multicombine=T,
                   .packages=c("picasso","glmnet","expm","doSNOW","snow","foreach","Matrix"),
                   .errorhandling="pass") %dopar% {
                      ii=which(taxaNames==refTaxa[i])

                     dataForEst=dataRecovTrans(data=data,independence=independence,
                                               identity=identity,ref=refTaxa[i],
                                               Mprefix=Mprefix,covsPrefix=covsPrefix,
                                               binPredInd=binPredInd)

                     x=as(dataForEst$xTildalong,"sparseMatrix")
                     y=as(dataForEst$UtildaLong,"sparseVector")
                     rm(dataForEst)
                     if(method=="lasso") {
                       Penal.i=runGlmnet(x=x,y=y,nPredics=nPredics)
                       MCPoriginSuccCV=0
                       nCV=0
                     }

                     if(method=="mcp") {
                       Penal.i=runPicasso(x=x,y=y,lambda=lambda,nPredics=nPredics,method="mcp",
                                          permutY=F,allFunc=allFunc)
                       MCPoriginSuccCV=Penal.i$successCV
                       nCV=Penal.i$nCV
                     }

                     rm(x,y)
                     gc()

                     if(is.null(Penal.i)) {
                       Penal.i$betaNoInt=as(rep(0,nNorm),"sparseVector")
                     }
                     BetaNoInt.i=as(Penal.i$betaNoInt,"sparseVector")
                     lambda=Penal.i$lambda
                     rm(Penal.i)
                     gc()

                     selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                     BetaNoIntWithRef.i=as(rep(0,nAlphaSelec),"sparseVector")

                     if (ii==1){
                       selection.i[-seq(1,nPredics)]=as(BetaNoInt.i!=0,"sparseVector")
                       BetaNoIntWithRef.i[-seq(1,nPredics)]=BetaNoInt.i
                     }
                     if (ii==nTaxa) {
                       selection.i[-seq((nAlphaSelec-nPredics-1),nAlphaSelec)]=as(BetaNoInt.i!=0,"sparseVector")
                       BetaNoIntWithRef.i[-seq((nAlphaSelec-nPredics-1),nAlphaSelec)]=BetaNoInt.i
                     }
                     if ((ii>1) & (ii<nTaxa)) {
                       selection.i[1:(nPredics*(ii-1))]=as(BetaNoInt.i[1:(nPredics*(ii-1))]!=0,"sparseVector")
                       selection.i[(nPredics*ii+1):nAlphaSelec]=as(BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]!=0,"sparseVector")
                       BetaNoIntWithRef.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                       BetaNoIntWithRef.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                     }
                     rm(BetaNoInt.i)
                     gc()

                     # create return vector
                     recturnVec=as(rep(0,(2*nAlphaSelec+3)),"sparseVector")
                     recturnVec[1]=nCV
                     recturnVec[2]=MCPoriginSuccCV
                     recturnVec[3]=lambda
                     recturnVec[4:(nAlphaSelec+3)]=selection.i
                     recturnVec[(nAlphaSelec+4):(2*nAlphaSelec+3)]=BetaNoIntWithRef.i
                     rm(selection.i,BetaNoIntWithRef.i)
                      return((recturnVec))
                   }
  registerDoSEQ()
  stopCluster(c1)
  gc()

  cat("Original parallel screen done","\n")

  scr1Resu<- lapply(scr1Resu, as, "sparseMatrix")
  scr1Resu=do.call(cbind, scr1Resu)

  # create count of for each predictor,each row is the count
  # of selection for a predictor

  countOfSelecForAPred=matrix(Matrix::rowSums(scr1Resu[4:(nAlphaSelec+3),,drop=F]),nrow=nPredics)
  testCovCountMat=countOfSelecForAPred[testCovInd,,drop=F]
  countOfSelecForAPred=as(matrix(Matrix::colSums(testCovCountMat),nrow=1),"sparseMatrix")
  betaAveOverRef=matrix(Matrix::rowSums(scr1Resu[(nAlphaSelec+4):(2*nAlphaSelec+3),,drop=F]),nrow=nPredics)

  rm(testCovInd)

  results$MCPoriginSuccCV=scr1Resu[2,]
  results$nCV=scr1Resu[1,1]
  results$lambda=scr1Resu[3,]

  rm(scr1Resu)
  gc()

  colnames(countOfSelecForAPred)=taxaNames
  rm(taxaNames)

  # return results
  results$testCovCountMat=testCovCountMat
  rm(testCovCountMat)
  results$countOfSelecForAPred=countOfSelecForAPred
  rm(countOfSelecForAPred)
  results$betaAveOverRef=betaAveOverRef
  rm(betaAveOverRef)
  return(results)
}

