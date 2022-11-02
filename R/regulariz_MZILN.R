Regulariz_MZILN <- function(data,
                            nRef,
                            sub_taxa,
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
                            phase1_taxon_num = 200,
                            trans_x_col = 200,
                            spar_cutoff = 10) {
  results <- list()
  regul.start.time <- proc.time()[3]
  
  nTestCov <- length(testCovInd)
  
  dataSparsCheck(data = data, Mprefix = Mprefix)
  
  # load abundance data info
  data.info <- dataInfo(
    data = data,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd
  )
  nSub <- data.info$nSub
  taxaNames <- data.info$taxaNames
  nPredics <- data.info$nPredics
  nTaxa <- data.info$nTaxa
  rm(data.info)
  
  refTaxa_reOrder <- microbName[taxaNames %in% refTaxa]
  
  qualifyData <- data
  
  
  if (length(binaryInd_test) > 0) {
    qualifyData <-
      data[rowSums(data[, taxaNames] > 0) >= 2, , drop = FALSE]
    allBinPred <- paste0(covsPrefix, binaryInd_test)
    nBinPred <- length(allBinPred)
    
    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax <-
      taxaNames[!microbName %in% (results$finalizedBootRefTaxon)]
    unbalanceTaxa <- c()
    unbalancePred <- c()
    for (i in AllTaxaNamesNoRefTax) {
      for (j in allBinPred) {
        twoColumns.ij <- qualifyData[, c(i, j)]
        nNonZero <- sum(twoColumns.ij[, 1] > 0)
        sumOfBin <- sum(twoColumns.ij[(twoColumns.ij[, 1] > 0), 2])
        if (sumOfBin %in% c(0, 1, (nNonZero - 1), nNonZero)) {
          unbalanceTaxa <- c(unbalanceTaxa, i)
          unbalancePred <- c(unbalancePred, j)
        }
      }
    }
    if (length(unbalanceTaxa) > 0) {
      unbalanceTaxa_ori_name <-
        microbName[unlist(lapply(unbalanceTaxa, function(x) {
          which(taxaNames %in% x)
        }))]
      unbalancePred_ori_name <-
        testCovInOrder[unlist(lapply(unbalancePred, function(x) {
          which(testCovInNewNam %in% x)
        }))]
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
  
  ##### Partition in phase 2 ####
  num_taxa_each <- ceiling(trans_x_col / (nPredics + 1))
  
  rowSpars <- apply(data[, taxaNames], 1, function(x) {
    sum(x == 0) / length(x)
  })
  
  meadianRowSpars <- min(median(rowSpars), 0.999)
  
  num_taxa_each <-
    max(ceiling(spar_cutoff / (1 - meadianRowSpars)),
        num_taxa_each)
  num_taxa_each <- min(num_taxa_each, nTaxa)
  
  # shuffle_seq <- sample(seq_len(length(taxaNames)))
  taxaNames_shuff <- taxaNames
  
  spar_each_taxon <-
    apply(data[, taxaNames_shuff], 2, function(x) {
      sum(x == 0) / length(x)
    })
  high_spar_taxon <- spar_each_taxon[spar_each_taxon > 0.7]
  low_spar_taxon <- spar_each_taxon[spar_each_taxon <= 0.7]
  
  num_cut <- floor((nTaxa) / num_taxa_each)
  
  high_spar_cut_taxon <-
    suppressWarnings(split(high_spar_taxon, seq_len(num_cut)))
  low_spar_cut_taxon <-
    suppressWarnings(split(low_spar_taxon, seq_len(num_cut)))
  
  taxa_sepname_list <- list()
  for (i in seq_len(num_cut)) {
    taxa_sepname_list[[i]] <-
      c(names(low_spar_cut_taxon[[i]]), names(high_spar_cut_taxon[[num_cut + 1 -
                                                                     i]]))
  }
  
  
  betaMatList <- list()
  CILowMatList <- list()
  CIUpMatList <- list()
  results$estiList <- list()
  for (iii in seq_len(nRef)) {
    time11 <- proc.time()[3]
    originTaxNam <- refTaxa_reOrder[iii]
    newRefTaxNam <- refTaxa[iii]
    # bootLassoAlpha_bon <- fdrRate
    results$estiList[[originTaxNam]] <- bootResuHDCI(
      data = data,
      refTaxa = newRefTaxNam,
      originRefTaxNam = originTaxNam,
      bootB = bootB,
      binPredInd = binaryInd,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      allFunc = allFunc,
      unbalanceTaxa_ori_name =
        unbalanceTaxa_ori_name,
      unbalancePred_ori_name =
        unbalancePred_ori_name,
      testCovInOrder = testCovInOrder,
      adjust_method = adjust_method,
      microbName = microbName,
      fwerRate = fdrRate,
      paraJobs = paraJobs,
      sequentialRun = sequentialRun,
      taxa_sepname_list_arg = taxa_sepname_list
    )
    
    time12 <- proc.time()[3]
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
          "unadj p-value")
      
      if (any(sub_taxa!="all")) {
        est_res_save_all <-
          est_res_save_all[est_res_save_all$taxon %in% sub_taxa, , drop = FALSE]
      }
      
      est_res_save_all$adj.p.value <-
        p.adjust(est_res_save_all$`unadj p-value`, adjust_method)
      rearrage_res_list[[j]] <- est_res_save_all
    }
    
    unorder_long <- do.call("rbind", rearrage_res_list)
    all_cov_list[[fin_ref_taxon_name[i]]] <-
      data.frame(unorder_long[stringr::str_order(unorder_long[, c("taxon")], decreasing = FALSE, numeric = TRUE),], row.names = NULL)
  }
  
  all_cov_list <- do.call("rbind", all_cov_list)
  rownames(all_cov_list) <- NULL
  all_cov_list$sig_ind <- all_cov_list$adj.p.value < fdrRate
  all_cov_list$sig_ind[is.na(all_cov_list$sig_ind)] <- FALSE
  
  
  
  
  results$full_results <- S4Vectors::DataFrame(all_cov_list)
  results$nSub <- nSub
  results$nTaxa <- nTaxa
  results$nPredics <- nPredics
  rm(data)
  # return results
  return(results)
}
