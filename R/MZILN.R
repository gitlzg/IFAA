##' Conditional regression for microbiome analysis based on multivariate zero-inflated logistic normal model
##'
##' For estimating and testing the associations of abundance ratios with covariates.
##' \loadmathjax
##'
##' Most of the time, users just need to feed the first three inputs to the function: `experiment_dat`, `refTaxa` and `allCov`. All other inputs can just take their default values.
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
##'
##' @param experiment_dat A SummarizedExperiment object containing microbiome data and covarites (see example on how to create a SummarizedExperiment object). The microbiome data can be 
##' absolute abundance or relative abundance with each column per sample and each row per taxon/OTU/ASV (or any other unit). No imputation is needed for zero-valued data points. 
##' The covarites data contains covariates and confounders with each row per sample and each column per variable. The covarites data has to be numeric or binary.
##' @param refTaxa Denominator taxa names specified by the user for the targeted ratios. This could be a vector of names.
##' @param allCov All covariates of interest (including confounders) for estimating and testing their associations with the targeted ratios. Default is 'NULL' meaning that all covariates in covData are of interest.
##' @param sampleID Name of the ID variable. Could be set to NULL. 
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `allCov`. Default is `0.15`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression. The default is `500`.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `TRUE`. It can be set to be "FALSE" to increase speed if there are multiple taxa in the argument 'refTaxa'.
##' @param verbose Whether the process message is printed out to the console. The default is TRUE.  
##' @param seed Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding.
##' @return A list with two elements. 
##'
##' - `full_results`: The main results for MZILN containing the estimation and testing results for all associations between all taxa ratios with refTaxan being the denominator and all covariates in `allCov`. It is a dataframe with each row representing an association, and ten columns named as "ref_tax", "taxon", "cov", "estimate", "SE.est", "CI.low", "CI.up", "adj.p.value", "unadj.p.value", and "sig_ind". The columns correspond to the denominator taxon, numerator taxon, covariate name, association estimates, standard error estimates, lower bound and upper bound of the 95% confidence interval, adjusted p value, and the indicator showing whether the association is significant after multiple testing adjustment.
##'
##' - `metadata`: The metadata is a list containing total time used in minutes, random seed used, FDR rate, and multiple testing adjustment method used.
##'
##' @examples
##' library(IFAA)
##' library(SummarizedExperiment)
##' 
##' ## load the example microbiome data. This could be relative abundance or absolute 
##' ## abundance data. If you have a csv or tsv file for the microbiome data, you 
##' ## can use read.csv() function or read.table() function in R to read the 
##' ## data file into R.
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' 
##' ## load the example covariates data. If you have a csv or tsv file for the 
##' ## covariates data, you can use read.csv() function or read.table() function 
##' ## in R to read the data file into R.
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' 
##' ## Merge microbiome data and covariate data by id, to avoid unmatching observations. 
##' data_merged<-merge(dataM,dataC,by="id",all=FALSE)
##' 
##' ## Seperate microbiome data and covariate data, drop id variable from the microbiome data
##' dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]
##' dataC_sub<-data_merged[,colnames(dataC)]
##' 
##' ## Create SummarizedExperiment object for inputs 
##' test_dat<-SummarizedExperiment(assays=list(counts=t(dataM_sub)), colData=dataC_sub)
##' 
##' ### Run MZILN function
##' results <- MZILN(experiment_dat = test_dat,
##'                 refTaxa=c("rawCount11"),
##'                 allCov=c("v1","v2","v3"),
##'                 sampleID=c("id"),
##'                 fdrRate=0.05)
##' ## to extract the results for all ratios with rawCount11 as the denominator:
##' summary_res<-results$full_results
##' ## to extract results for the ratio of a specific taxon (e.g., rawCount45) over rawCount11:
##' target_ratio=summary_res[summary_res$taxon=="rawCount45",]
##' ## to extract all of the ratios having significant associations:
##' sig_ratios=subset(summary_res,sig_ind==TRUE)


##'
##'
##' @references Li et al.(2018) Conditional Regression Based on a Multivariate Zero-Inflated Logistic-Normal Model for Microbiome Relative Abundance Data. Statistics in Biosciences 10(3): 587-608
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Liu et al.(2020) A bootstrap lasso + partial ridge method to construct confidence intervals for parameters in high-dimensional sparse linear models. Statistica Sinica

##' @export
##' @md

MZILN = function(experiment_dat,
                 refTaxa,
                 allCov = NULL,
                 sampleID = NULL,
                 adjust_method = "BY",
                 fdrRate = 0.15,
                 paraJobs = NULL,
                 bootB = 500,
                 taxDropThresh = 0,
                 standardize = FALSE,
                 sequentialRun = TRUE,
                 verbose = TRUE,
                 seed = 1) {
  allFunc = allUserFunc()
  
  results = list()
  
  start.time = proc.time()[3]
  
  MicrobData <- data.frame(t(assays(experiment_dat)$counts))
  MicrobData$ID_char_1234 <- rownames(MicrobData)
  CovData <- data.frame(colData(experiment_dat))
  CovData$ID_char_1234 <- rownames(CovData)
  linkIDname <- "ID_char_1234"
  
  if (verbose) {
    runMeta = metaData(
      MicrobData = MicrobData,
      CovData = CovData,
      linkIDname = linkIDname,
      testCov = allCov,
      taxDropThresh = taxDropThresh,
      standardize = standardize
    )
  } else {
    runMeta = suppressMessages(
      metaData(
        MicrobData = MicrobData,
        CovData = CovData,
        linkIDname = linkIDname,
        testCov = allCov,
        taxDropThresh = taxDropThresh,
        standardize = standardize
      )
    )
  }
  
  
  data = runMeta$data
  covariatesData = runMeta$covariatesData
  binaryInd = runMeta$binaryInd
  covsPrefix = runMeta$covsPrefix
  Mprefix = runMeta$Mprefix
  testCovInd = runMeta$testCovInd
  testCovInOrder = runMeta$testCovInOrder
  testCovInNewNam = runMeta$testCovInNewNam
  ctrlCov = runMeta$ctrlCov
  microbName = runMeta$microbName
  newMicrobNames = runMeta$newMicrobNames
  results$covriateNames = runMeta$xNames
  rm(runMeta)
  
  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]
  
  nRef = length(refTaxa)
  refTaxa_newNam = newMicrobNames[microbName %in% refTaxa]
  
  
  
  if (length(refTaxa) > 0) {
    if (sum(refTaxa %in% microbName) != length(refTaxa)) {
      stop(
        "Error: One or more of the specified reference taxa have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels."
      )
    }
  }
  
  if (verbose) {
    results$analysisResults = Regulariz_MZILN(
      data = data,
      nRef = nRef,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      adjust_method = adjust_method,
      paraJobs = paraJobs,
      fdrRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      seed = seed
    )
  } else {
    results$analysisResults = suppressMessages(
      Regulariz_MZILN(
        data = data,
        nRef = nRef,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        adjust_method = adjust_method,
        paraJobs = paraJobs,
        fdrRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc,
        seed = seed
      )
    )
  }
  
  rm(data)
  
  if (length(sampleID)>0) {
    covariatesData <- merge(CovData[,c(sampleID,linkIDname)],covariatesData,all=FALSE)
  } 
  covariatesData <- covariatesData[,!colnames(covariatesData) %in% c(linkIDname)]
  
  
  rm(testCovInOrder, ctrlCov, microbName)
  
  totalTimeMins = (proc.time()[3] - start.time) / 60
  message("The entire analysis took ", round(totalTimeMins, 2), " minutes")
  
  
  output_se_obj <-
    list(
      full_results  = results$analysisResults$full_results,
      metadata = list(
        covariatesData = covariatesData,
        totalTimeMins = totalTimeMins,
        seed = seed,
        fdrRate = fdrRate,
        adjust_method = adjust_method
      )
    )
  
  return(output_se_obj)
}
