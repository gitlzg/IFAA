getScrResu <- function(data,
                       testCovInd,
                       testCovInOrder,
                       testCovInNewNam,
                       nRef,
                       paraJobs,
                       refTaxa,
                       goodIndeCutPerc = 0.33,
                       sequentialRun,
                       allFunc = allFunc,
                       refReadsThresh,
                       SDThresh,
                       SDquantilThresh,
                       balanceCut,
                       Mprefix,
                       covsPrefix,
                       binPredInd,
                       adjust_method) {
  results <- list()
  # run permutation
  scrParal <- runScrParal(
    data = data,
    testCovInd = testCovInd,
    testCovInOrder = testCovInOrder,
    testCovInNewNam = testCovInNewNam,
    nRef = nRef,
    paraJobs = paraJobs,
    refTaxa = refTaxa,
    sequentialRun = sequentialRun,
    allFunc = allFunc,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd,
    adjust_method = adjust_method
  )

  selecCountOverall <- scrParal$countOfSelecForAPred
  selecEstOverall <- scrParal$estOfSelectForAPred

  selecCountMatIndv <- scrParal$testCovCountMat
  selecEstMatIndv <- scrParal$testEstMat

  taxaNames <- scrParal$taxaNames
  goodRefTaxaCandi <- scrParal$goodRefTaxaCandi

  nTaxa <- scrParal$nTaxa
  nPredics <- scrParal$nPredics
  nTestCov <- scrParal$nTestCov
  results$refTaxa <- scrParal$refTaxa
  rm(scrParal)
  gc()

  if (nTestCov == 1) {
    results$selecCountMatIndv <- selecCountOverall
    results$selecEstMatIndv <- selecEstOverall
  }

  if (nTestCov > 1) {
    results$selecCountMatIndv <- selecCountMatIndv
    results$selecEstMatIndv <- selecEstMatIndv

    rm(selecCountMatIndv)
  }

  goodIndpRefTaxWithCount <-
    selecCountOverall[1, (colnames(selecCountOverall) %in%
      goodRefTaxaCandi)]
  goodIndpRefTaxWithEst <-
    selecEstOverall[1, (colnames(selecEstOverall) %in%
      goodRefTaxaCandi)]

  restRefTaxWithCount <-
    selecCountOverall[1, !(colnames(selecCountOverall) %in%
      goodRefTaxaCandi)]
  restRefTaxWithEst <-
    selecEstOverall[1, !(colnames(selecEstOverall) %in% goodRefTaxaCandi)]

  sort_goodIndpRefTaxWithCount <-
    goodIndpRefTaxWithCount[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))]
  sort_goodIndpRefTaxWithEst <-
    abs(goodIndpRefTaxWithEst[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))])

  sort_restRefTaxWithCount <-
    restRefTaxWithCount[order(restRefTaxWithCount, abs(restRefTaxWithEst))]
  sort_restRefTaxWithEst <-
    abs(restRefTaxWithEst[order(restRefTaxWithCount, abs(restRefTaxWithEst))])


  results$finalIndpRefTax <-
    names(c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount))[seq_len(2)]

  results$selecCountOverall <- selecCountOverall

  results$goodIndpRefTaxWithCount <-
    c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount)
  results$goodIndpRefTaxWithEst <-
    c(sort_goodIndpRefTaxWithEst, sort_restRefTaxWithEst)

  results$goodRefTaxaCandi <- goodRefTaxaCandi
  rm(goodRefTaxaCandi)

  return(results)
}
