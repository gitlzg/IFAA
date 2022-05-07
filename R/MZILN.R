##' Conditional regression for microbiome analysis based on multivariate zero-inflated logistic normal model
##'
##' For estimating and testing the associations of abundance ratios with covariates.
##' \loadmathjax
##'
##' Most of the time, users just need to feed the first six inputs to the function: `MicrobData`, `CovData`,  `linkIDname`, `targetTaxa`, `refTaxa` and `allCov`. All other inputs can just take their default values.
##' The regression model for `MZILN()` can be expressed as follows:
##' \mjdeqn{\log\bigg(\frac{\mathcal{Y}_i^k}{\mathcal{Y}_i^{K+1}}\bigg)|\mathcal{Y}_i^k>0,\mathcal{Y}_i^{K+1}>0=\alpha^{0k}+\mathcal{X}_i^T\alpha^k+\epsilon_i^k,\hspace{0.2cm}k=1,...,K}{}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{} in subject \mjeqn{i}{} in the entire
##' ecosystem.
##' - \mjeqn{\mathcal{Y}_i^{K+1}}{} is the reference taxon (specified by user).
##' - \mjeqn{\mathcal{X}_i}{} is the covariate matrix for all covariates including confounders.
##' - \mjeqn{\alpha^k}{} is the regression coefficients along with their 95% confidence intervals that will be estimated by the `MZILN()` function.
##'
##' High-dimensional \mjeqn{X_i}{} is handled by regularization.
##'
##' @param MicrobData Microbiome data matrix containing microbiome absolute abundance or relative abundance with each row
##' per sample and each column per taxon/OTU/ASV (or any other unit). It should contain an id variable to
##' be linked with the id variable in the covariates data: `CovData`. This argument can
##' take directory path. For example, `MicrobData="C://...//microbiomeData.tsv"`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row
##' per sample and each column per variable. Any categorical variable should be converted into dummy variables in this data matrix unless it can be treated as a continuous variable.
##' It should also contain an id variable to
##' be linked with the id variable in the microbiome data: `MicrobData`. This argument can
##' take directory path. For example, `CovData="C://...//covariatesData.tsv"`.
##'
##'
##' @param linkIDname The common variable name of the id variable in both `MicrobData` and `CovData`. The two data sets will be merged by this id variable.
##' @param targetTaxa The numerator taxa names specified by users for the targeted ratios. Default is NULL in which case all taxa are numerator taxa (except the taxa in the argument 'refTaxa').
##' @param refTaxa Denominator taxa names specified by the user for the targeted ratios. This could be a vector of names.
##' @param allCov All covariates of interest (including confounders) for estimating and testing their associations with the targeted ratios. Default is 'NULL' meaning that all covariates in covData are of interest.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `allCov`. Default is `0.15`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression. The default is `500`.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `TRUE`. It can be set to be "FALSE" to increase speed if there are multiple taxa in the argument 'refTaxa'.
##' @param seed Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding.
##' @return A list containing the estimation results.
##'
##'##' - `targettaxa_result_list`: A list containing estimating results for the targeted ratios. Only available when targetTaxa is non-empty.
##' - `sig_results`: A list containing estimating results for all significant ratios
##' - `covariatesData`: A dataset containing all covariates used in the analyses.
##'
##' @examples
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' \donttest{
##' results <- MZILN(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 targetTaxa = "rawCount6",
##'                 refTaxa=c("rawCount11"),
##'                 allCov=c("v1","v2","v3"),
##'                 fdrRate=0.15)
##'
##' }
##'
##'
##'
##' @references Li et al.(2018) Conditional Regression Based on a Multivariate Zero-Inflated Logistic-Normal Model for Microbiome Relative Abundance Data. Statistics in Biosciences 10(3): 587-608
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Liu et al.(2020) A bootstrap lasso + partial ridge method to construct confidence intervals for parameters in high-dimensional sparse linear models. Statistica Sinica

##' @export
##' @md



MZILN=function(
  MicrobData,
  CovData,
  linkIDname,
  targetTaxa=NULL,
  refTaxa,
  allCov=NULL,
  adjust_method="BY",
  fdrRate=0.15,
  paraJobs=NULL,
  bootB=500,
  taxDropThresh=0,
  standardize=FALSE,
  sequentialRun=TRUE,
  seed=1
){
  allFunc=allUserFunc()

  results=list()

  start.time = proc.time()[3]
  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=allCov,
                   taxDropThresh=taxDropThresh,standardize=standardize)

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
  newMicrobNames=runMeta$newMicrobNames
  results$covriateNames=runMeta$xNames
  rm(runMeta)

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]

  nRef=length(refTaxa)
  refTaxa_newNam=newMicrobNames[microbName%in%refTaxa]

  if ((length(targetTaxa)>0)) {
    if (sum(targetTaxa%in%microbName)!=length(targetTaxa)) {
      stop("Error: One or more of the target taxa have no sequencing reads or are not in the data set.
      Double check the names of the target taxa and their
           sparsity levels.")
    }
  }

  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels.")
    }
  }

  results$analysisResults=Regulariz_MZILN(data=data,nRef=nRef,testCovInd=testCovInd,
                                          testCovInOrder=testCovInOrder,testCovInNewNam=testCovInNewNam,
                                          microbName=microbName,
                                          binaryInd=binaryInd,
                                          binaryInd_test=binaryInd_test,
                                          covsPrefix=covsPrefix,Mprefix=Mprefix,
                                          targetTaxa=targetTaxa,
                                          refTaxa=refTaxa_newNam,adjust_method=adjust_method,
                                          paraJobs=paraJobs,fdrRate=fdrRate,
                                          bootB=bootB,
                                          sequentialRun=sequentialRun,
                                          allFunc=allFunc,seed=seed
  )
  rm(data)

  results$sig_results<-results$analysisResults$sig_results
  results$targettaxa_result_list<-results$analysisResults$targettaxa_result_list
  results$testCov=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  results$bootB=bootB
  results$paraJobs=paraJobs
  results$nRef=nRef

  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }

  rm(testCovInOrder,ctrlCov,microbName)

  totalTimeMins = (proc.time()[3] - start.time)/60
  message("The entire analysis took ",round(totalTimeMins,2), " minutes")

  results$totalTimeMins=totalTimeMins

  return(results)
}
