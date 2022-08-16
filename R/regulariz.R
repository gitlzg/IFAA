Regulariz <- function(data,
                      testCovInd,
                      testCovInOrder,
                      testCovInNewNam,
                      microbName,
                      nRef,
                      nRefMaxForEsti,
                      refTaxa,
                      paraJobs,
                      binaryInd,
                      binaryInd_test,
                      covsPrefix,
                      Mprefix,
                      fwerRate,
                      bootB,
                      sequentialRun,
                      allFunc,
                      refReadsThresh,
                      SDThresh,
                      SDquantilThresh,
                      balanceCut,
                      adjust_method,
                      phase1_taxon_num = 200,
                      trans_x_col = 200,
                      spar_cutoff = 10) {
  results <- list()
  regul.start.time <- proc.time()[3]

  nTestCov <- length(testCovInd)

  dataSparsCheck(data = data, Mprefix = Mprefix)


  data.info <- dataInfo(
    data = data,
    qualifyRefTax = TRUE,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd
  )
  
  nSub <- data.info$nSub
  taxaNames <- data.info$taxaNames
  nPredics <- data.info$nPredics
  nTaxa <- data.info$nTaxa
  predNames <- data.info$predNames
  
  results$goodRefTaxaCandi <- data.info$goodRefTaxaCandi
  
  num_taxon_phase1 <- max(length(refTaxa), phase1_taxon_num)
  num_taxon_phase1 <- min(num_taxon_phase1, nTaxa)
  phase1_taxon_sample_pool <-
    results$goodRefTaxaCandi[!(results$goodRefTaxaCandi %in% refTaxa)]
  if (length(phase1_taxon_sample_pool) > (num_taxon_phase1 - length(refTaxa))) {
    phase1_taxon_sample1 <-
      c(
        sample(
          phase1_taxon_sample_pool,
          num_taxon_phase1 - length(refTaxa)
        ),
        refTaxa
      )
  } else {
    phase1_taxon_sample1 <- c(phase1_taxon_sample_pool, refTaxa)
  }
  phase1_taxon_sample <- phase1_taxon_sample1
  
  
  if (length(phase1_taxon_sample1) < num_taxon_phase1) {
    not_good_candi <-
      taxaNames[!((taxaNames %in% results$goodRefTaxaCandi) |
                    (taxaNames %in% refTaxa))]
    non_zero_per_taxon <- colSums(data[, not_good_candi, drop = FALSE] > 0)
    phase1_taxon_sample2 <-
      names(sort(non_zero_per_taxon, decreasing = TRUE)[seq_len(num_taxon_phase1 -
                                                                  length(phase1_taxon_sample1))])
    phase1_taxon_sample <-
      c(phase1_taxon_sample1, phase1_taxon_sample2)
  }
  
  data_sub_phase1 <- data[, c(phase1_taxon_sample, predNames)]
  
  rm(data.info)

  regul.start.time <- proc.time()[3]
  message("Start Phase 1 analysis")

  t1phase1 <- proc.time()[3]

  selectRegroup <- getScrResu(
    data = data_sub_phase1,
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
    binPredInd = binaryInd,
    adjust_method = adjust_method
  )
  nRef_smaller <- max(2, ceiling(nRef / 2))
  while_loop_ind <- FALSE
  loop_num <- 0

  t2phase1 <- proc.time()[3]
  minu <- round((t2phase1 - t1phase1) / 60, 2)

  message(
    "33 percent of phase 1 analysis has been done and it took ",
    minu,
    " minutes"
  )

  while (while_loop_ind == FALSE) {
    if (loop_num >= 2) {
      break
    }
    loop_num <- loop_num + 1

    refTaxa_smaller <-
      head(names(selectRegroup$goodIndpRefTaxWithCount),
        n = nRef_smaller
      )

    fin_ref_1 <- selectRegroup$finalIndpRefTax
    ref_taxa_1 <- selectRegroup$refTaxa
    selectRegroup <- getScrResu(
      data = data_sub_phase1,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      nRef = nRef_smaller,
      paraJobs = paraJobs,
      refTaxa = refTaxa_smaller,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      refReadsThresh = refReadsThresh,
      SDThresh = SDThresh,
      SDquantilThresh = SDquantilThresh,
      balanceCut = balanceCut,
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binaryInd,
      adjust_method = adjust_method
    )


    fin_ref_2 <- selectRegroup$finalIndpRefTax
    ref_taxa_2 <- selectRegroup$refTaxa
    while_loop_ind <-
      identical(fin_ref_1, fin_ref_2) ||
        identical(ref_taxa_1, ref_taxa_2)

    t3phase1 <- proc.time()[3]
    minu <- round((t3phase1 - t1phase1) / 60, 2)

    if (while_loop_ind == FALSE) {
      message(
        round(100 * (loop_num + 1) / 3, 0),
        " percent of phase 1 analysis has been done and it took ",
        minu,
        " minutes"
      )
    }
    if (while_loop_ind == TRUE) {
      message(
        "100 percent of phase 1 analysis has been done and it took ",
        minu,
        " minutes"
      )
    }
  }

  results$selecCountOverall <- selectRegroup$selecCountOverall
  microbName_ind <-
    vapply(colnames(results$selecCountOverall), function(x) {
      which(taxaNames %in% x)
    }, FUN.VALUE = 1)

  colnames(results$selecCountOverall) <- microbName[microbName_ind]

  results$selecCountMatIndv <- selectRegroup$selecCountMatIndv
  finalIndpRefTax <-
    microbName[taxaNames %in% (selectRegroup$finalIndpRefTax)]

  results$goodIndpRefTaxWithCount <-
    selectRegroup$goodIndpRefTaxWithCount
  names(results$goodIndpRefTaxWithCount) <-
    microbName[unlist(lapply(names(selectRegroup$goodIndpRefTaxWithCount), function(x) {
      which(taxaNames %in% x)
    }))]

  results$goodIndpRefTaxWithEst <-
    selectRegroup$goodIndpRefTaxWithEst
  names(results$goodIndpRefTaxWithEst) <-
    microbName[unlist(lapply(names(selectRegroup$goodIndpRefTaxWithEst), function(x) {
      which(taxaNames %in% x)
    }))]

  results$goodRefTaxaCandi <-
    microbName[taxaNames %in% (selectRegroup$goodRefTaxaCandi)]
  results$randomRefTaxa <-
    microbName[taxaNames %in% (selectRegroup$refTaxa)]
  rm(selectRegroup)

  goodIndpRefTaxNam <- names(results$goodIndpRefTaxWithCount)

  MCPExecuTime <- (proc.time()[3] - regul.start.time) / 60
  results$phase1ExecuTime <- MCPExecuTime

  results$finalIndpRefTax <- finalIndpRefTax

  startT <- proc.time()[3]
  message("Start Phase 2 parameter estimation")

  qualifyData <- data


  if (length(binaryInd_test) > 0) {
    qualifyData <-
      data[rowSums(data[, taxaNames] > 0) >= 2, , drop = FALSE]
    allBinPred <- paste0(covsPrefix, binaryInd_test)
    nBinPred <- length(allBinPred)

    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax <-
      taxaNames[!microbName %in% finalIndpRefTax]
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
    rm(allBinPred)
  } else {
    unbalanceTaxa_ori_name <- NULL
    unbalancePred_ori_name <- NULL
  }

  allRefTaxNam <- unique(c(finalIndpRefTax, goodIndpRefTaxNam))
  nGoodIndpRef <- length(allRefTaxNam)
  results$allRefTaxNam <- allRefTaxNam

  results$nRefUsedForEsti <- min(nGoodIndpRef, nRefMaxForEsti)
  
  ##### Partition in phase 2 ####
  num_taxa_each <- ceiling(trans_x_col / (nPredics + 1))
  
  rowSpars <- apply(data[, taxaNames], 1, function(x) {
    sum(x == 0) / length(x)
  })
  
  meadianRowSpars <- min(median(rowSpars), 0.999)
  
  num_taxa_each <-
    max(
      ceiling(spar_cutoff / (
        1 - meadianRowSpars
      )),
      num_taxa_each
    )
  num_taxa_each <- min(num_taxa_each, nTaxa)
  
  shuffle_seq <- sample(seq_len(length(taxaNames)))
  taxaNames_shuff <- taxaNames[shuffle_seq]
  
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
  
  

  results$estiList <- list()
  for (iii in seq_len(results$nRefUsedForEsti)) {
    time11 <- proc.time()[3]
    originTaxNam <- allRefTaxNam[iii]
    newRefTaxNam <- taxaNames[microbName %in% originTaxNam]
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
      fwerRate = fwerRate,
      paraJobs = paraJobs,
      sequentialRun = sequentialRun,
      taxa_sepname_list_arg = taxa_sepname_list
    )
    time12 <- proc.time()[3]
    if (iii < (results$nRefUsedForEsti)) {
      message(
        round((
          100 * iii / (results$nRefUsedForEsti)
        ), 2),
        " percent of Phase 2 is done and it took ",
        round((time12 - time11) / 60, 3),
        " minutes"
      )
    }
  }

  endT <- proc.time()[3]

  message(
    "Entire Phase 2 parameter estimation done and took ",
    round((endT - startT) / 60, 3),
    " minutes."
  )

  ### calculate mean ###
  fin_ref_taxon_name <- names(results$estiList)
  all_cov_list_sep <- list()


  for (i in seq_len(length(fin_ref_taxon_name))) {
    save_list_temp <- results$estiList[[i]]$all_cov_list
    rearrage_res_list <- list()
    for (j in testCovInOrder) {
      est_res_save_all <-
        cbind(
          save_list_temp$est_save_mat[j, ],
          save_list_temp$se_mat[j, ],
          save_list_temp$CI_low_mat[j, ],
          save_list_temp$CI_up_mat[j, ],
          save_list_temp$p_value_save_mat[j, ]
        )
      est_res_save_all <-
        data.frame(
          fin_ref_taxon_name[i],
          rownames(est_res_save_all),
          j,
          est_res_save_all
        )

      colnames(est_res_save_all) <-
        c(
          "ref_tax",
          "taxon",
          "cov",
          "estimate",
          "SE est",
          "CI low",
          "CI up",
          "adj p-value"
        )
      rearrage_res_list[[j]] <- est_res_save_all
    }

    unorder_long <- DescTools::DoCall("rbind", rearrage_res_list)
    all_cov_list_sep[[fin_ref_taxon_name[i]]] <-
      data.frame(unorder_long[stringr::str_order(unorder_long[, c("taxon")], decreasing = FALSE, numeric = TRUE), ], row.names = NULL)
  }

  all_cov_list_sep <- DescTools::DoCall("rbind", all_cov_list_sep)
  rownames(all_cov_list_sep) <- NULL
  all_cov_list_sep$sig_ind <-
    all_cov_list_sep$adj.p.value < fwerRate
  results$all_cov_list_sep <- S4Vectors::DataFrame(all_cov_list_sep)
  ##### Mean ####
  all_cov_list <- list()

  for (i in seq_len(length(fin_ref_taxon_name))) {
    all_cov_list[[fin_ref_taxon_name[i]]] <-
      results$estiList[[i]]$all_cov_list
  }
  results$all_cov_list <- all_cov_list


  ref_taxon_name <- names(all_cov_list)

  exclu_1 <-
    !colnames(all_cov_list[[1]]$est_save_mat) %in% ref_taxon_name[2]
  exclu_2 <-
    !colnames(all_cov_list[[2]]$est_save_mat) %in% ref_taxon_name[1]


  est_save_mat_mean <-
    (all_cov_list[[1]]$est_save_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$est_save_mat[, exclu_2,
      drop =
        FALSE
    ]) / 2
  se_mat_mean <-
    (all_cov_list[[1]]$se_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$se_mat[, exclu_2,
      drop =
        FALSE
    ]) / 2
  CI_low_mat_mean <-
    (all_cov_list[[1]]$CI_low_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$CI_low_mat[, exclu_2,
      drop =
        FALSE
    ]) / 2
  CI_up_mat_mean <-
    (all_cov_list[[1]]$CI_up_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$CI_up_mat[, exclu_2,
      drop =
        FALSE
    ]) / 2


  p_value_unadj_mean <-
    t(apply(est_save_mat_mean / se_mat_mean, 1, function(x) {
      (1 - pnorm(abs(
        x
      ))) * 2
    }))
  p_value_adj_mean <-
    t(apply(p_value_unadj_mean, 1, function(x) {
      p.adjust(x, method = adjust_method)
    }))
  colname_use <- colnames(est_save_mat_mean)

  sig_ind <-
    which(p_value_adj_mean <= fwerRate,
      arr.ind = TRUE,
      useNames = FALSE
    )

  full_results <- list()
  for (j in testCovInOrder) {
    est_res_save_all <-
      data.frame(
        colname_use,
        j,
        est_save_mat_mean[j, ],
        se_mat_mean[j, ],
        CI_low_mat_mean[j, ],
        CI_up_mat_mean[j, ],
        p_value_adj_mean[j, ],
        row.names = NULL
      )
    colnames(est_res_save_all) <-
      c(
        "taxon",
        "cov",
        "estimate",
        "SE est",
        "CI low",
        "CI up",
        "adj p-value"
      )
    fin_ref_taxon_dat <-
      data.frame(
        taxon = fin_ref_taxon_name,
        cov = j,
        estimate = 0,
        SE.est = NA,
        CI.low = NA,
        CI.up = NA,
        adj.p.value = 1
      )
    full_results[[j]] <-
      rbind(data.frame(est_res_save_all), fin_ref_taxon_dat)
  }
  full_results <- DescTools::DoCall("rbind", full_results)
  rownames(full_results) <- NULL
  full_results <-
    full_results[stringr::str_order(full_results$taxon,
      decreasing = FALSE,
      numeric = TRUE
    ), ]
  full_results$sig_ind <- full_results$adj.p.value < fwerRate
  full_results$sig_ind[is.na(full_results$sig_ind)] <- FALSE
  results$full_results <- S4Vectors::DataFrame(full_results)

  results$nTaxa <- nTaxa
  results$nPredics <- nPredics
  results$fin_ref_taxon_name <- fin_ref_taxon_name

  rm(data)

  # return results

  results$nRef <- nRef
  return(results)
}

#
