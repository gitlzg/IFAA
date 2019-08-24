
##' @export


IFAA=function(
  Microbdata,
  CovData,
  linkIDname,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=F,
  ctrlMany=F,
  nRef,
  nPermu,
  refTaxa=NULL,
  reguMethod=c("mcp"), ## "lasso" or "mcp"
  fwerRate=0.30, #0.30
  boot=T,
  paraJobs=NULL,
  bootB=5,
  bootLassoAlpha=0.05,
  allFunc=allUserFunc(),
  refReadsThresh=0.3,
  SDThresh=0.0001,
  SDquantilThresh=0.0001,
  balanceCut=0.3,
  seed=1
){
  results=list()
  start.time = proc.time()[3]
  runMeta=metaData(Microbdata=Microbdata,CovData=CovData,
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

  if(length(refTaxa)>1){
    stop("Too many reference taxa specified. Only one at a time is allowed.")
  }
  if(length(refTaxa)==1 ){
    if(!refTaxa%in%microbName)
      stop("The specified reference taxon is not in the data set,
           double check the name of the reference taxon.")
  }
  results$analysisResults=Regulariz(data=data,testCovInd=testCovInd,
                                    testCovInOrder=testCovInOrder,
                                    testCovInNewNam=testCovInNewNam,
                                    microbName=microbName,nRef=nRef,
                                    nPermu=nPermu,binaryInd=binaryInd,
                                    covsPrefix=covsPrefix,Mprefix=Mprefix,
                                    refTaxa=refTaxa,paraJobs=paraJobs,
                                    reguMethod=reguMethod,fwerRate=fwerRate,
                                    boot=boot,bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                    allFunc=allFunc,refReadsThresh=refReadsThresh,
                                    SDThresh=SDThresh,
                                    SDquantilThresh=SDquantilThresh,
                                    balanceCut=balanceCut,seed=seed
  )
  rm(data)

  results$testCovInOrder=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  rm(testCovInOrder,ctrlCov,microbName)

  totalTimeMins = (proc.time()[3] - start.time)/60
  cat("The entire analysis took",totalTimeMins, "minutes","\n")

  return(results)
  }
