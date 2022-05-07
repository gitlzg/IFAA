

runBootLassoHDCI=function(
  x,
  y,
  nfolds=10,
  lambdaOPT=NULL,
  zeroSDCut=10^(-20),
  correCut=0.996,
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
                    standardize=FALSE,parallel.boot=TRUE,
                    ncores.boot=paraJobs,alpha=bootLassoAlpha)
  parallel::stopCluster(c3)

  # convert to a sparse vector format from sparse matrix format

  beta.LPR=bootResu$Beta.LPR
  betaCI.LPR=bootResu$interval.LPR
  rm(bootResu)
  # transform vector back
  if(length(xWithNearZeroSd)>0){

    betaTransLasso.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                         unSelectList=sort(xWithNearZeroSd),newBetaNoInt=beta.LPR)
    beta.LPR=betaTransLasso.L$finalBeta

    betaTransCIlow.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                         unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI.LPR[1,])
    betaCIlow.LPR=betaTransCIlow.L$finalBeta

    betaTransCIhi.L=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                        unSelectList=sort(xWithNearZeroSd),newBetaNoInt=betaCI.LPR[2,])
    betaCIhi.LPR=betaTransCIhi.L$finalBeta

    rm(betaTransLasso.L,betaTransCIlow.L,betaTransCIhi.L)
  } else{
    betaCIlow.LPR=betaCI.LPR[1,]
    betaCIhi.LPR=betaCI.LPR[2,]
  }

  rm(xWithNearZeroSd)


  results$Beta.LPR=beta.LPR
  results$interval.LPR=rbind(betaCIlow.LPR,betaCIhi.LPR)



  return(results)
}
