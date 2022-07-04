Regulariz_MZILN = function(data,
                           nRef,
                           testCovInd,
                           testCovInOrder,
                           testCovInNewNam,
                           microbName,
                           refTaxa,
                           adjust_method,
                           paraJobs,
                           binaryInd,
                           binaryInd_test,
                           covsPrefix,
                           Mprefix,
                           fdrRate,
                           bootB,
                           sequentialRun,
                           allFunc = allFunc,
                           seed) {
  results = list()
  regul.start.time = proc.time()[3]
  
  nTestCov = length(testCovInd)
  
  dataSparsCheck(data = data, Mprefix = Mprefix)
  
  # load abundance data info
  data.info = dataInfo(
    data = data,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd
  )
  nSub = data.info$nSub
  taxaNames = data.info$taxaNames
  nPredics = data.info$nPredics
  nTaxa = data.info$nTaxa
  rm(data.info)
  
  refTaxa_reOrder = microbName[taxaNames %in% refTaxa]
  
  qualifyData = data
  
  
  if (length(binaryInd_test) > 0) {
    qualifyData = data[rowSums(data[, taxaNames] > 0) >= 2, , drop = FALSE]
    allBinPred = paste0(covsPrefix, binaryInd_test)
    nBinPred = length(allBinPred)
    
    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax = taxaNames[!microbName %in% (results$finalizedBootRefTaxon)]
    unbalanceTaxa = c()
    unbalancePred <- c()
    for (i in AllTaxaNamesNoRefTax) {
      for (j in allBinPred) {
        twoColumns.ij = qualifyData[, c(i, j)]
        nNonZero = sum(twoColumns.ij[, 1] > 0)
        sumOfBin = sum(twoColumns.ij[(twoColumns.ij[, 1] > 0), 2])
        if (sumOfBin %in% c(0, 1, (nNonZero - 1), nNonZero)) {
          unbalanceTaxa = c(unbalanceTaxa, i)
          unbalancePred <- c(unbalancePred, j)
        }
      }
    }
    if (length(unbalanceTaxa) > 0) {
      unbalanceTaxa_ori_name <-
        microbName[unlist(lapply(unbalanceTaxa, function(x)
          which(taxaNames %in% x)))]
      unbalancePred_ori_name <-
        testCovInOrder[unlist(lapply(unbalancePred, function(x)
          which(testCovInNewNam %in% x)))]
    } else {
      unbalanceTaxa_ori_name <- NULL
      unbalancePred_ori_name <- NULL
    }
    # unestimableTaxa=unique(unbalanceTaxa)
    rm(allBinPred)
  } else {
    unbalanceTaxa_ori_name <- NULL
    unbalancePred_ori_name <- NULL
  }
  
  
  betaMatList = list()
  CILowMatList = list()
  CIUpMatList = list()
  results$estiList = list()
  for (iii in seq_len(nRef)) {
    time11 = proc.time()[3]
    originTaxNam = refTaxa_reOrder[iii]
    newRefTaxNam <- refTaxa[iii]
    bootLassoAlpha_bon <- fdrRate
    results$estiList[[originTaxNam]] = bootResuHDCI(
      data = data,
      refTaxa = newRefTaxNam,
      originRefTaxNam = originTaxNam,
      bootB = bootB,
      bootLassoAlpha = bootLassoAlpha_bon,
      binPredInd = binaryInd,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      unbalanceTaxa_ori_name =
        unbalanceTaxa_ori_name,
      unbalancePred_ori_name =
        unbalancePred_ori_name,
      testCovInOrder = testCovInOrder,
      adjust_method = adjust_method,
      microbName = microbName,
      fwerRate = fdrRate,
      paraJobs = paraJobs,
      seed = seed
    )
    
    time12 = proc.time()[3]
    message(
      "Estimation done for the ",
      iii,
      "th denominator taxon: ",
      refTaxa_reOrder[iii],
      " and it took ",
      round((time12 - time11) / 60, 2),
      " minutes"
    )
    
    
  }
  
  fin_ref_taxon_name <- names(results$estiList)
  all_cov_sig_list <- list()
  all_cov_list <- list()
  
  for (i in seq_len(length(fin_ref_taxon_name))) {
    all_cov_sig_list[[fin_ref_taxon_name[i]]] <-
      results$estiList[[i]]$sig_list_each
    save_list_temp <- results$estiList[[i]]$all_cov_list
    rearrage_res_list <- list()
    for (j in testCovInOrder) {
      est_res_save_all <-
        cbind(
          save_list_temp$est_save_mat[j,],
          save_list_temp$se_mat[j,],
          save_list_temp$CI_low_mat[j,],
          save_list_temp$CI_up_mat[j,],
          save_list_temp$p_value_save_mat[j,]
        )
      est_res_save_all <-
        data.frame(fin_ref_taxon_name[i],
                   rownames(est_res_save_all),
                   j,
                   est_res_save_all)
      
      colnames(est_res_save_all) <-
        c("ref_tax",
          "taxon",
          "cov",
          "estimate",
          "SE est",
          "CI low",
          "CI up",
          "adj p-value")
      est_res_save_all$unadj.p.value<-(1 - pnorm(abs(est_res_save_all$estimate / est_res_save_all$`SE est`))) * 2
      rearrage_res_list[[j]] <- est_res_save_all
    }
    
    unorder_long <- do.call("rbind", rearrage_res_list)
    all_cov_list[[fin_ref_taxon_name[i]]] <-
      data.frame(unorder_long[stringr::str_order(unorder_long[, c("taxon")], decreasing = FALSE, numeric = T), ], row.names = NULL)
  }
  
  all_cov_list <- do.call("rbind", all_cov_list)
  rownames(all_cov_list) <- NULL
  all_cov_list$sig_ind <- all_cov_list$adj.p.value < fdrRate
  
  
  
  
  
  results$full_results <- DataFrame(all_cov_list)
  results$nSub = nSub
  results$nTaxa = nTaxa
  results$nPredics = nPredics
  rm(data)
  # return results
  return(results)
}
