

runBootLassoHDCI=function(
  x,
  y,
  nPredics,
  nTaxa,
  nfolds=10,
  lambdaOPT=NULL,
  refTaxaPosition,
  zeroSDCut=0,
  correCut=0.996,
  standardize,
  paraJobs,
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

  df.cor=suppressWarnings(corSparse(x))
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

  # message(paste("Number of x with near-zero sd: ", nearZeroSd))

  cvStartTime= proc.time()[3]

  c3 <- parallel::makeCluster(paraJobs)
  doParallel::registerDoParallel(c3)

  if(length(seed)>0){
    set.seed(as.numeric(seed)+10^2)
    parallel::clusterSetRNGStream(cl=c3,(as.numeric(seed)+10^3))
  }

  bootResu=bootLOPR(x=as.matrix(x),y=as.vector(y),B=bootB,nfolds=nfolds,
                    standardize=standardize,parallel.boot=TRUE,
                    ncores.boot=paraJobs,alpha=bootLassoAlpha)
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
