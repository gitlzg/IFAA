#---------------------------------------------------------------------------------------
## Original screen function to select independent reference taxa
#---------------------------------------------------------------------------------------
originDataScreen <- function(data,
                             testCovInd,
                             nRef,
                             paraJobs,
                             refTaxa,
                             maxDimensionScr = 434 * 5 * 10^5,
                             sequentialRun,
                             allFunc,
                             Mprefix,
                             covsPrefix,
                             binPredInd,
                             adjust_method,
                             subsamp_cut = 3,
                             dfmax_div = 3,
                             nlambda_num = 10) {
  results <- list()

  # load data info
  basicInfo <- dataInfo(
    data = data,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd
  )

  taxaNames <- basicInfo$taxaNames
  nTaxa <- basicInfo$nTaxa
  nPredics <- basicInfo$nPredics
  rm(basicInfo)
  gc()

  nNorm <- nTaxa - 1
  nAlphaNoInt <- nPredics * nNorm
  nAlphaSelec <- nPredics * nTaxa

  countOfSelec <- rep(0, nAlphaSelec)

  # overwrite nRef if the reference taxon is specified
  nRef <- length(refTaxa)

  startT1 <- proc.time()[3]
  if (length(paraJobs) == 0) {
    availCores <- parallelly::availableCores()
    if (is.numeric(availCores)) {
      paraJobs <- max(1, parallelly::availableCores() - 2)
    }
    if (!is.numeric(availCores)) {
      paraJobs <- 1
    }
  }

  if (!sequentialRun) {
    message(
      paraJobs,
      " parallel jobs are registered for the analysis."
    )
  }

  cl <- parallel::makeCluster(paraJobs)

  parallel::clusterExport(
    cl = cl,
    varlist = allFunc[c(1, 2)],
    envir = environment()
  )

  doParallel::registerDoParallel(cl)

  if (sequentialRun) {
    foreach::registerDoSEQ()
  }

  # start parallel computing


  scr1Resu <- foreach(
    i = seq_len(nRef),
    .multicombine = TRUE,
    .errorhandling = "pass"
  ) %dorng% {

    # for(i in seq_len(nRef)){

    ii <- which(taxaNames == refTaxa[i])

    dataForEst <- dataRecovTrans(
      data = data,
      ref = refTaxa[i],
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd,
      contCovStd = TRUE
    )

    xTildLongTild.i <- dataForEst$xTildalong
    yTildLongTild.i <- dataForEst$UtildaLong
    rm(dataForEst)
    gc()

    maxSubSamplSiz <- min(50000, floor(maxDimensionScr /
      ncol(xTildLongTild.i)))
    nToSamplFrom <- nrow(xTildLongTild.i)

    subSamplK <- ceiling(nToSamplFrom / maxSubSamplSiz)
    if (subSamplK == 1) {
      maxSubSamplSiz <- nToSamplFrom
    }

    nRuns <- ceiling(subSamplK / subsamp_cut)

    for (k in seq_len(nRuns)) {
      rowToKeep <- sample(nToSamplFrom, maxSubSamplSiz)
      x <- xTildLongTild.i[rowToKeep, ]
      y <- yTildLongTild.i[rowToKeep]

      dfmax <-
        ceiling(min((nPredics + 1) * nNorm / dfmax_div, length(rowToKeep) /
          dfmax_div))
      
      pmax <- ceiling(min((nPredics + 1) * nNorm, length(rowToKeep)))
      
      lasso.est <- glmnet::glmnet(
        x = x,
        y = as.vector(y),
        dfmax = dfmax,
        pmax = pmax,
        family = "gaussian",
        intercept = TRUE,
        standardize = FALSE,
        nlambda = nlambda_num
      )

      rm(x, y)
      gc()

      betaMat <- as.matrix(lasso.est$beta)
      rm(lasso.est)
      gc()

      lasso_beta <- betaMat[, ncol(betaMat)]
      rm(betaMat)
      gc()

      betaNoInt <-
        lasso_beta[-seq(1, length(lasso_beta), by = (nPredics + 1))]
      rm(lasso_beta)

      BetaNoInt.k <- (0 + (betaNoInt != 0))
      EstNoInt.k <- abs(betaNoInt)

      if (k == 1) {
        BetaNoInt.i <- BetaNoInt.k
        EstNoInt.i <- EstNoInt.k
      }
      if (k > 1) {
        BetaNoInt.i <- BetaNoInt.i + BetaNoInt.k
        EstNoInt.i <- EstNoInt.i + EstNoInt.k
      }
      rm(BetaNoInt.k, EstNoInt.k)
    }
    rm(k, xTildLongTild.i)
    gc()
    BetaNoInt.i <- BetaNoInt.i / nRuns
    EstNoInt.i <- EstNoInt.i / nRuns

    selection.i <- rep(0, nAlphaSelec)
    coef.i <- rep(0, nAlphaSelec)
    if (ii == 1) {
      selection.i[-seq(1, nPredics)] <- BetaNoInt.i
      coef.i[-seq(1, nPredics)] <- EstNoInt.i
    }
    if (ii == nTaxa) {
      selection.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <-
        BetaNoInt.i
      coef.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <-
        EstNoInt.i
    }
    if ((ii > 1) & (ii < nTaxa)) {
      selection.i[seq_len((nPredics * (ii - 1)))] <-
        BetaNoInt.i[seq_len((nPredics *
          (ii - 1)))]
      selection.i[(nPredics * ii + 1):nAlphaSelec] <-
        BetaNoInt.i[(nPredics *
          (ii - 1) + 1):nAlphaNoInt]
      coef.i[seq_len((nPredics * (ii - 1)))] <-
        EstNoInt.i[seq_len((nPredics *
          (ii - 1)))]
      coef.i[(nPredics * ii + 1):nAlphaSelec] <-
        EstNoInt.i[(nPredics *
          (ii - 1) + 1):nAlphaNoInt]
    }
    rm(BetaNoInt.i)
    # create return vector
    recturnlist <- list()
    recturnlist[[1]] <- selection.i
    recturnlist[[2]] <- coef.i
    rm(selection.i, yTildLongTild.i, coef.i)
    return(recturnlist)
  }
  parallel::stopCluster(cl)
  gc()
  rm(data)

  endT <- proc.time()[3]

  selecList <- list()
  for (i in seq_len(nRef)) {
    selecList[[i]] <- scr1Resu[[i]][[1]]
  }

  estList <- list()
  for (i in seq_len(nRef)) {
    estList[[i]] <- scr1Resu[[i]][[2]]
  }

  rm(scr1Resu)

  scr1ResuSelec <- DescTools::DoCall(cbind, selecList)
  scr1ResuEst <- DescTools::DoCall(cbind, estList)

  rm(selecList, estList)

  # create count of selection for individual testCov
  countOfSelecForAllPred <-
    matrix(rowSums(scr1ResuSelec), nrow = nPredics)
  EstOfAllPred <- matrix(rowMeans(scr1ResuEst), nrow = nPredics)

  testCovCountMat <-
    countOfSelecForAllPred[testCovInd, , drop = FALSE]
  testEstMat <- EstOfAllPred[testCovInd, , drop = FALSE]
  rm(
    scr1ResuSelec,
    testCovInd,
    countOfSelecForAllPred,
    EstOfAllPred
  )

  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred <- matrix(colSums(testCovCountMat), nrow = 1)
  estOfSelectForAPred <- matrix(colSums(testEstMat), nrow = 1)
  gc()

  colnames(countOfSelecForAPred) <- taxaNames
  colnames(estOfSelectForAPred) <- taxaNames
  rm(taxaNames)

  # return results
  results$testCovCountMat <- testCovCountMat
  results$testEstMat <- testEstMat
  rm(testCovCountMat, testEstMat)
  results$countOfSelecForAPred <- countOfSelecForAPred
  results$estOfSelectForAPred <- estOfSelectForAPred
  rm(countOfSelecForAPred, estOfSelectForAPred)
  return(results)
}
