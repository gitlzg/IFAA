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
##' @param reguMethod regularization approach used in phase 1 of the algorithm. Take value `"mcp"` or `"lasso"`, default is `"mcp"`.
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