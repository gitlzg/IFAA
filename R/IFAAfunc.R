##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of microbiome with covariates
##'
##' Most of the time, users just need to feed the first three inputs to the function: `experiment_dat`, `testCov` and `ctrlCov`. All other inputs can just take their default values.
##' To model the association, the following equation is used:
##'
##' \loadmathjax
##'
#' \mjdeqn{\log(\mathcal{Y}_i^k)|\mathcal{Y}_i^k>0=\beta^{0k}+X_i^T\beta^k+W_i^T\gamma^k+Z_i^Tb_i+\epsilon_i^k,\hspace{0.2cm}k=1,...,K+1}{}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{} in subject \mjeqn{i}{} in the entire ecosystem.
##' - \mjeqn{X_i}{} is the covariate matrix.
##' - \mjeqn{W_i}{} is the confounder matrix.
##' - \mjeqn{Z_i}{} is the design matrix for random effects.
##' - \mjeqn{\beta^k}{} is the regression coefficients that will be estimated and tested with the `IFAA()` function.
##'
##' The challenge in microbiome analysis is that \mjeqn{\mathcal{Y}_i^k}{} can not be observed. What is observed is its small proportion:
##' \mjeqn{Y_i^k=C_i\mathcal{Y}^k_i}{}, where \mjeqn{C_i}{} is an unknown number between 0
##' and 1 that denote the observed proportion.
##'
##' The IFAA method can successfully addressed this challenge. The `IFAA()` will estimate the parameter
##' \mjeqn{\beta^k}{} and their 95% confidence intervals. High-dimensional \mjeqn{X_i}{} is handled by
##' regularization.
##'
##' @param experiment_dat A SummarizedExperiment object containing microbiome data and covarites (see example on how to create a SummarizedExperiment object). The microbiome data can be absolute abundance or relative abundance 
##' with each column per sample and each row per taxon/OTU/ASV (or any other unit). No imputation is needed for zero-valued data points. The covarites data contains covariates and confounders with each row per sample and each 
##' column per variable. The covarites data has to be numeric or binary.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param sampleID Name of the sample ID variable in the data. This is not required. Default is NULL. 
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `FALSE`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nRefMaxForEsti The maximum number of final reference taxa used in phase 2. The default is `2`.
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `testCov`. Default is `0.15`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. The default is `500`.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `FALSE`. This argument could be useful for debug.
##' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param SDThresh The threshold of standard deviations of sequencing reads for been chosen as the reference taxon in phase 2. The default is `0.05` which means the standard deviation of sequencing reads should be at least `0.05` in order to be chosen as reference taxon.
##' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as reference taxon. The default is `0`.
##' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default number is `0.2` which means at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
##' @param verbose Whether the process message is printed out to the console. The default is TRUE.  
##' @param seed Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding.
##' @return A list containing 2 elements
##' \itemize{
##' \item {`full_results`: The main results for IFAA containing the estimation and testing results for all associations between all taxa and all test covariates in `testCov`. It is a dataframe with each row 
##' representing an association, and eight columns named as "taxon", "cov", "estimate", "SE.est", "CI.low", "CI.up", "adj.p.value", and "sig_ind". The columns correspond to taxon name, covariate name, association estimates, 
##' standard error estimates, lower bound and upper bound of the 95% confidence interval, adjusted p value, and the indicator showing whether the association is significant after multiple testing adjustment.
##' }
##' \item {`metadata`: The metadata is a list. 
##' \itemize{
##' \item {`covariatesData`: A dataset containing covariates and confounders used in the analyses.}
##' \item {`final_ref_taxon`: The final 2 reference taxon used for analysis.}
##' \item {`ref_taxon_count`: The counts of selection for the associations of all taxa with test covariates in Phase 1.} 
##' \item {`totalTimeMins`: The average magnitude estimates for the associations of all taxa with test covariates in Phase 1.} 
##' \item {`ref_taxon_est`: Total time used for the entire analysis.} 
##' \item {`seed`: The seed used for the analysis for reproducibility.} 
##' \item {`fdrRate`: FDR rate used for the analysis.} 
##' \item {`adjust_method`: Multiple testing adjust method used for the analysis.} 

##' }
##' }
##' }
##' @examples
##' library(IFAA)
##' library(SummarizedExperiment)
##' 
##' ## If you already have a SummarizedExperiment format data, you can ignore the data processing steps below.
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
##' ## Seperate microbiome data and covariate data, drop id variable from microbiome data
##' dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]
##' dataC_sub<-data_merged[,colnames(dataC)]
##' 
##' ## Create SummarizedExperiment object 
##' test_dat<-SummarizedExperiment(assays=list(counts=t(dataM_sub)), colData=dataC_sub)
##' 
##' ## If you already have a SummarizedExperiment format data, you can ignore the above steps.
##' 
##' ## run IFAA function
##' results <- IFAA(experiment_dat = test_dat,
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"),
##'                 sampleID = c("id"),
##'                 fdrRate = 0.05)
##'
##' ## to extract all results:
##' summary_res<-results$full_results
##' ## to extract significant results:
##' sig_results=subset(summary_res,sig_ind==TRUE)
##'
##'
##'
##' @references Li et al.(2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Liu et al.(2020) A bootstrap lasso + partial ridge method to construct confidence intervals for parameters in high-dimensional sparse linear models. Statistica Sinica
##' @importFrom foreach `%dopar%` foreach
##' @importFrom methods as
##' @importFrom parallelly availableCores
##' @importFrom Matrix Diagonal Matrix sparseVector
##' @importFrom HDCI bootLOPR
##' @importFrom qlcMatrix corSparse
##' @import mathjaxr
##' @import glmnet
##' @import stats
##' @import utils
##' @importFrom SummarizedExperiment assays colData SummarizedExperiment
##' @importFrom stringr str_order
##' @importFrom S4Vectors DataFrame
##' @importFrom speedglm speedlm
##' @import DescTools
##' @import MatrixExtra
##' @export
##' @md

## arg deleted: nperm, x1permt, reguMethod, method, bootLassoAlpha


IFAA = function(experiment_dat,
                testCov = NULL,
                ctrlCov = NULL,
                sampleID = NULL,
                testMany = TRUE,
                ctrlMany = FALSE,
                nRef = 40,
                nRefMaxForEsti = 2,
                refTaxa = NULL,
                adjust_method = "BY",
                fdrRate = 0.15,
                paraJobs = NULL,
                bootB = 500,
                standardize = FALSE,
                sequentialRun = FALSE,
                refReadsThresh = 0.2,
                taxDropThresh = 0,
                SDThresh = 0.05,
                SDquantilThresh = 0,
                balanceCut = 0.2,
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
      testCov = testCov,
      ctrlCov = ctrlCov,
      testMany = testMany,
      ctrlMany = ctrlMany,
      taxDropThresh = taxDropThresh,
      standardize = standardize
    )
  } else {
    runMeta = suppressMessages(
      metaData(
        MicrobData = MicrobData,
        CovData = CovData,
        linkIDname = linkIDname,
        testCov = testCov,
        ctrlCov = ctrlCov,
        testMany = testMany,
        ctrlMany = ctrlMany,
        taxDropThresh = taxDropThresh,
        standardize = standardize
      )
    )
  }
  
  
  data = runMeta$data
  covariatesData<-runMeta$covariatesData
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
  
  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]
  
  
  rm(runMeta)
  
  if(length(refTaxa)>0){
    if(length(unique(refTaxa))!=length(refTaxa)){
      message("Duplicated names in refTaxa are removed")
      refTaxa=unique(refTaxa)
    }
    
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      refTaxa<-refTaxa[refTaxa%in%microbName]
      message("One or more of the specified reference taxa in phase 1 have no sequencing reads
      or are not in the data set.")
    }
  }
  
  
  if (nRefMaxForEsti < 2) {
    nRefMaxForEsti <- 2
    warning("Warning: Needs at least two final reference taxon for estimation.")
  }
  
  
  if(nRef>(length(microbName))){
    nRef=length(microbName)
    message("The number of reference taxa in Phase 1 is set to be equal to the total number
           of taxa in the data because it cannot exceed that.")
  }
  
  refTaxa_newNam = newMicrobNames[microbName %in% refTaxa]
  
  if (verbose) {
    results$analysisResults = Regulariz(
      data = data,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      nRef = nRef,
      nRefMaxForEsti = nRefMaxForEsti,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      paraJobs = paraJobs,
      adjust_method = adjust_method,
      fwerRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      refReadsThresh = refReadsThresh,
      SDThresh = SDThresh,
      SDquantilThresh = SDquantilThresh,
      balanceCut = balanceCut,
      seed = seed
    )
  } else {
    results$analysisResults = suppressMessages(
      Regulariz(
        data = data,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        nRef = nRef,
        nRefMaxForEsti = nRefMaxForEsti,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        paraJobs = paraJobs,
        adjust_method = adjust_method,
        fwerRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc,
        refReadsThresh = refReadsThresh,
        SDThresh = SDThresh,
        SDquantilThresh = SDquantilThresh,
        balanceCut = balanceCut,
        seed = seed
      )
    )
  }
  
  
  
  rm(data)
  
  
  
  
  
  totalTimeMins = (proc.time()[3] - start.time) / 60
  message("The entire analysis took ", round(totalTimeMins, 2), " minutes")
  
  if (length(seed) == 1) {
    results$seed = seed
  } else{
    results$seed = "No seed used."
  }
  
  if (length(sampleID)>0) {
    covariatesData <- merge(CovData[,c(sampleID,linkIDname)],covariatesData,by=linkIDname,all=FALSE)
  } 
  covariatesData <- covariatesData[,!colnames(covariatesData) %in% c(linkIDname)]
  
  
  
  output_se_obj <-
    list(
      full_results  = results$analysisResults$full_results,
      metadata = list(
        covariatesData = covariatesData,
        final_ref_taxon = results$analysisResults$fin_ref_taxon_name,
        ref_taxon_count =
          results$analysisResults$goodIndpRefTaxWithCount,
        ref_taxon_est = results$analysisResults$goodIndpRefTaxWithEst,
        totalTimeMins = totalTimeMins,
        seed = seed,
        fdrRate = fdrRate,
        adjust_method = adjust_method
      )
    )
  
  
  
  return(output_se_obj)
}


