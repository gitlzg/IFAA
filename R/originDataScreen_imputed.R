#---------------------------------------------------------------------------------------
## Original screen function to select independent reference taxa with imputed dataset
#---------------------------------------------------------------------------------------
originDataScreen_imputed = function(data,
                                    testCovInd,
                                    nRef,
                                    paraJobs,
                                    refTaxa,
                                    maxDimensionScr = 434 * 5 * 10 ^ 5,
                                    sequentialRun,
                                    allFunc,
                                    Mprefix,
                                    covsPrefix,
                                    binPredInd,
                                    adjust_method = "fdr",
                                    fwerRate = 0.25,
                                    seed) {
  results = list()
  
  # load data info
  basicInfo = dataInfo(
    data = data,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd
  )
  
  taxaNames = basicInfo$taxaNames
  nTaxa = basicInfo$nTaxa
  nPredics = basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm = nTaxa - 1
  nAlphaNoInt = nPredics * nNorm
  nAlphaSelec = nPredics * nTaxa
  
  countOfSelec = rep(0, nAlphaSelec)
  
  MVarNamLength = nchar(Mprefix)
  
  micros = lapply(substr(colnames(data), 1, MVarNamLength), function(x) {
    grep(Mprefix, x)
  })
  microPositions = which(micros == 1)
  
  xVarNamLength = nchar(covsPrefix)
  predics = lapply(substr(colnames(data), 1, xVarNamLength), function(x) {
    grep(covsPrefix, x)
  })
  predPositions = which(predics == 1)
  predNames = colnames(data)[predPositions]
  nPredics = length(predNames)
  
  # extract x data
  xData = data[, predNames]
  
  # overwrite nRef if the reference taxon is specified
  nRef = length(refTaxa)
  
  startT1 = proc.time()[3]
  if (length(paraJobs) == 0) {
    availCores = parallelly::availableCores()
    if (is.numeric(availCores))
      paraJobs = max(1, parallelly::availableCores() - 2)
    if (!is.numeric(availCores))
      paraJobs = 1
  }
  
  if (!sequentialRun) {
    message(
      paraJobs,
      " parallel jobs are registered for analyzing ",
      nRef,
      " reference taxa in Phase 1"
    )
  }
  
  cl <- parallel::makeCluster(paraJobs)
  
  parallel::clusterExport(cl = cl,
                          varlist = allFunc,
                          envir = environment())
  
  doParallel::registerDoParallel(cl)
  
  if (length(seed) > 0) {
    set.seed(as.numeric(seed) + 10 ^ 4)
    parallel::clusterSetRNGStream(cl = cl, (as.numeric(seed) + 10 ^ 5))
  }
  
  if (sequentialRun) {
    foreach::registerDoSEQ()
  }
  
  #start parallel computing
  scr1Resu = foreach(
    i = seq_len(nRef),
    .multicombine = TRUE,
    .packages = c("glmnet", "Matrix"),
    .errorhandling = "pass"
  ) %dopar% {
    ref = refTaxa[i]
    
    otherTaxaNames = taxaNames[(taxaNames != ref)]
    microb_ref_dat <- log(data[, otherTaxaNames] / data[, ref])
    
    
    
    est_temp_save <- numeric(length(taxaNames) * nPredics)
    p_value_temp_save <- numeric(length(taxaNames) * nPredics)
    
    k <- 0
    for (j in taxaNames) {
      k <- k + 1
      if (j == ref) {
        est_temp_save[(k * 3 - 2):(k * 3)] <- 0
        p_value_temp_save[(k * 3 - 2):(k * 3)] <- NA
      } else {
        lm_res <- lm(microb_ref_dat[, j] ~ xData)
        
        full_name_coef <- names(lm_res$coefficients)
        valid_coef <- summary(lm_res)$coefficients
        bootResu <- matrix(nrow = length(full_name_coef), ncol = 4)
        rownames(bootResu) <- full_name_coef
        bootResu[rownames(bootResu) %in% rownames(valid_coef), ] <-
          valid_coef
        bootResu <- bootResu[-1,]
        
        est_temp_save[(k * 3 - 2):(k * 3)] <- abs(bootResu[, 1])
        p_value_temp_save[(k * 3 - 2):(k * 3)] <- bootResu[, 4]
      }
      
    }
    
    p_value_adj_temp_save <-
      p.adjust(p_value_temp_save, adjust_method,
               n = sum(!is.na(p_value_temp_save)))
    p_value_adj_temp_save[is.na(p_value_adj_temp_save)] <- 1
    
    CovCountall <- p_value_adj_temp_save < fwerRate
    
    countOfSelecForAllPred <-
      matrix(
        CovCountall,
        nrow = nPredics,
        ncol = length(taxaNames),
        byrow = FALSE
      )
    colnames(countOfSelecForAllPred) <- taxaNames
    EstOfAllPred <-
      matrix(
        est_temp_save,
        nrow = nPredics,
        ncol = length(taxaNames),
        byrow = FALSE
      )
    colnames(EstOfAllPred) <- taxaNames
    
    
    
    selection.i <- countOfSelecForAllPred[testCovInd,]
    coef.i <- EstOfAllPred[testCovInd,]
    
    recturnlist = list()
    recturnlist[[1]] = selection.i
    recturnlist[[2]] <- coef.i
    return(recturnlist)
  }
  
  parallel::stopCluster(cl)
  gc()
  rm(data)
  endT = proc.time()[3]
  
  
  selecList = list()
  for (i in seq_len(nRef)) {
    selecList[[i]] = scr1Resu[[i]][[1]]
  }
  
  
  
  
  estList = list()
  for (i in seq_len(nRef)) {
    estList[[i]] = scr1Resu[[i]][[2]]
  }
  
  
  
  testCovCountMat <- base::Reduce("+", selecList)
  testEstMat <- base::Reduce("+", estList)
  
  rm(selecList, estList)
  
  
  
  # create count of selection for individual testCov
  
  
  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred <-
    matrix(Matrix::colSums(testCovCountMat), nrow =
             1)
  estOfSelectForAPred <-
    matrix(Matrix::colSums(testEstMat), nrow = 1)
  gc()
  
  colnames(countOfSelecForAPred) = taxaNames
  colnames(estOfSelectForAPred) <- taxaNames
  rm(taxaNames)
  
  # return results
  results$testCovCountMat = testCovCountMat
  results$testEstMat <- testEstMat
  rm(testCovCountMat, testEstMat)
  results$countOfSelecForAPred = countOfSelecForAPred
  results$estOfSelectForAPred <- estOfSelectForAPred
  rm(countOfSelecForAPred, estOfSelectForAPred)
  return(results)
}
