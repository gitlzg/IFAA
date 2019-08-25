##' @export


runBootLassoHDCI=function(
  x=x,
  y=y,
  nPredics,
  nTaxa,
  nfolds=10,
  lambdaOPT=NULL,
  refTaxaPosition,
  zeroSDCut=0,
  standardize=F,
  bootB,
  bootLassoAlpha
){

  results=list()

  # check if there are x variables with zero variance
  nObsAll=length(y)
  nBeta=ncol(x)

  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)

  # remove near constant columns
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd]
  }
  rm(sdX)
  nearZeroSd=length(xWithNearZeroSd)

  # print(paste("Number of x with near-zero sd: ", nearZeroSd))

  cvStartTime= proc.time()[3]

  availCores=availableCores()
  if(is.numeric(availCores))ncores.boot=max(1,availableCores()-2)
  if(!is.numeric(availCores))ncores.boot=1

  c3 <- parallel::makeCluster(ncores.boot)
  doParallel::registerDoParallel(c3)

  bootResu=bootLOPR(x=x,y=as.vector(y),B=bootB,nfolds=nfolds,
                    standardize=standardize,parallel.boot=T,
                    ncores.boot=ncores.boot,alpha=bootLassoAlpha)
  parallel::stopCluster(c3)

  # convert to a sparse vector format from sparse matrix format
  beta=as(bootResu$Beta.LPR,"sparseVector")
  betaCI=as(bootResu$interval.LOPR,"sparseMatrix")
  rm(bootResu)

  # filling zeros
  if(length(xWithNearZeroSd)>0){
    betaTransLasso=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                       unSelectList=sort(xWithNearZeroSd),newBetaNoInt=beta)
    beta=betaTransLasso$finalBeta

    betaTransCIlow=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                       unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI[1,])
    betaCIlow=betaTransCIlow$finalBeta

    betaTransCIhi=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                      unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI[2,])
    betaCIhi=betaTransCIhi$finalBeta

    rm(betaTransLasso,betaTransCIlow,betaTransCIhi)
  } else{
    betaCIlow=betaCI[1,]
    betaCIhi=betaCI[2,]
  }

  rm(xWithNearZeroSd)

  BetaNoInt.i=as(beta[-seq(1,length(beta),by=(nPredics+1))],"sparseVector")
  betaCIlow.i=as(betaCIlow[-seq(1,length(betaCIlow),by=(nPredics+1))],"sparseVector")
  betaCIhi.i=as(betaCIhi[-seq(1,length(betaCIhi),by=(nPredics+1))],"sparseVector")

  nAlphaSelec=nPredics*nTaxa
  nAlphaNoInt=nPredics*(nTaxa-1)

  BetaNoIntWithRef=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIlow=as(rep(0,nAlphaSelec),"sparseVector")
  betaCIhi=as(rep(0,nAlphaSelec),"sparseVector")

  ii=refTaxaPosition

  if (ii==1){
    BetaNoIntWithRef[-seq(1,nPredics)]=BetaNoInt.i
    betaCIlow[-seq(1,nPredics)]=betaCIlow.i
    betaCIhi[-seq(1,nPredics)]=betaCIhi.i
  }
  if (ii==nTaxa){
    BetaNoIntWithRef[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
    betaCIlow[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIlow.i
    betaCIhi[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=betaCIhi.i
  }
  if ((ii>1) & (ii<nTaxa)) {
    BetaNoIntWithRef[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
    BetaNoIntWithRef[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]

    betaCIlow[1:(nPredics*(ii-1))]=betaCIlow.i[1:(nPredics*(ii-1))]
    betaCIlow[(nPredics*ii+1):nAlphaSelec]=betaCIlow.i[(nPredics*(ii-1)+1):nAlphaNoInt]

    betaCIhi[1:(nPredics*(ii-1))]=betaCIhi.i[1:(nPredics*(ii-1))]
    betaCIhi[(nPredics*ii+1):nAlphaSelec]=betaCIhi.i[(nPredics*(ii-1)+1):nAlphaNoInt]
  }

  rm(BetaNoInt.i,betaCIlow.i,betaCIhi.i)


  results$beta=BetaNoIntWithRef
  results$betaCIlow=betaCIlow
  results$betaCIhi=betaCIhi

  rm(BetaNoIntWithRef,betaCIlow,betaCIhi)

  return(results)
}
