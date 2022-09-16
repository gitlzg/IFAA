#---------------------------------------------------------------------------------------
## bootstrap results function for Lasso OLS HDCI
#---------------------------------------------------------------------------------------

bootResuHDCI <- function(data,
                         refTaxa,
                         originRefTaxNam,
                         maxDimension = 434 * 5 * 10^5,
                         bootB,
                         bootLassoAlpha = 0.05,
                         binPredInd,
                         covsPrefix,
                         Mprefix,
                         allFunc,
                         unbalanceTaxa_ori_name,
                         unbalancePred_ori_name,
                         testCovInOrder,
                         adjust_method,
                         microbName,
                         fwerRate,
                         paraJobs,
                         sequentialRun,
                         taxa_sepname_list_arg,
                         subsamp_cut = 3) {
  results <- list()
  # load data info
  basicInfo <- dataInfo(
    data = data,
    binPredInd = binPredInd,
    covsPrefix = covsPrefix,
    Mprefix = Mprefix
  )
  taxaNames <- basicInfo$taxaNames
  ii <- which(basicInfo$taxaNames %in% refTaxa)
  predNames <- basicInfo$predNames
  nTaxa <- basicInfo$nTaxa
  nPredics <- basicInfo$nPredics
  rm(basicInfo)
  gc()

  nNorm <- nTaxa - 1
  nAlphaNoInt <- nPredics * nNorm
  nAlphaSelec <- nPredics * nTaxa

  countOfSelec <- rep(0, nAlphaSelec)
  resultsByRefTaxon <- list()

  taxa_sepname_list_noref <- lapply(taxa_sepname_list_arg, function(x) x[x!=refTaxa])
  taxa_sepname_list<-lapply(taxa_sepname_list_noref, function(x) c(x,refTaxa))
  
  taxa_sepname_list_noref_vec <- unlist(taxa_sepname_list_noref)
  taxa_sepname_numeric <-
    as.numeric(sub(
      pattern = Mprefix,
      replacement = "",
      x = taxa_sepname_list_noref_vec
    ))


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
      " parallel jobs are registered for analyzing reference taxa in Phase 2"
    )
  }

  cl <- parallel::makeCluster(paraJobs)

  #
  parallel::clusterExport(
    cl = cl,
    varlist = allFunc,
    envir = environment()
  )

  doParallel::registerDoParallel(cl)

  if (sequentialRun) {
    foreach::registerDoSEQ()
  }
  
  i <- numeric(0)
  
  phase2res <- foreach(
    i = seq_len(length(taxa_sepname_list)),
    .multicombine = TRUE,
    .errorhandling = "pass"
  ) %dorng% {
    data_sub_taxa <- data[, c(taxa_sepname_list[[i]], predNames)]
    dataForEst <- dataRecovTrans(
      data = data_sub_taxa,
      ref = refTaxa,
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd
    )
    x <- dataForEst$xTildalong
    y <- dataForEst$UtildaLong
    rm(dataForEst, data_sub_taxa)

    xCol <- ncol(x)

    maxSubSamplSiz <- min(100000, floor(maxDimension / xCol))
    nToSamplFrom <- length(y)
    subSamplK <- ceiling(nToSamplFrom / maxSubSamplSiz)
    if (subSamplK == 1) {
      maxSubSamplSiz <- nToSamplFrom
    }

    nRuns <- (ceiling(subSamplK / subsamp_cut))

    if (dim(x)[1] > (dim(x)[2])) {
      for (k in seq_len(nRuns)) {
        rowToKeep <- sample(nToSamplFrom, maxSubSamplSiz)
        xSub <- x[rowToKeep, ]
        ySub <- y[rowToKeep]
        bootResu <- lm_sparse(x = xSub, y = ySub)

        if (k == 1) {
          bootResu_k <- bootResu
        } else {
          bootResu_k <- bootResu_k + bootResu
        }
        gc()
      }
      rm(x, y)
      gc()


      fin_ref_taxon_name <- originRefTaxNam
      nTestcov <- length(testCovInOrder)


      boot_est <- bootResu_k[, 1] / nRuns
      se_est_all <- bootResu_k[, 2] / nRuns



      ref_taxon_name <- originRefTaxNam
      p_value_save_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      est_save_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_up_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_low_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      se_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))

      for (ii in seq_len(nTestcov)) {
        se_est <-
          se_est_all[seq(ii + 1, length(se_est_all), nPredics + 1)]
        boot_est_par <-
          boot_est[seq(ii + 1, length(boot_est), nPredics + 1)]

        p_value_unadj <- (1 - pnorm(abs(boot_est_par / se_est))) * 2

        boot_est_CI_low <- boot_est_par - 1.96 * se_est
        boot_est_CI_up <- boot_est_par + 1.96 * se_est

        p_value_adj <- p.adjust(p_value_unadj, adjust_method)
        p_value_save_mat[ii, ] <- p_value_unadj
        est_save_mat[ii, ] <- boot_est_par
        CI_low_mat[ii, ] <- boot_est_CI_low
        CI_up_mat[ii, ] <- boot_est_CI_up
        se_mat[ii, ] <- se_est
      }
    } else {
      for (k in seq_len(nRuns)) {
        rowToKeep <- sample(nToSamplFrom, maxSubSamplSiz)
        xSub <- x[rowToKeep, ]
        ySub <- y[rowToKeep]


        bootResu <- runBootLassoHDCI(
          x = as.matrix(xSub),
          y = as.vector(ySub),
          bootB = bootB,
          bootLassoAlpha = bootLassoAlpha
        )

        rm(xSub, ySub)
        gc()
        if (k == 1) {
          boot_est.k <- bootResu$Beta.LPR
          boot_CI.k <- bootResu$interval.LPR
        } else {
          boot_est.k <- boot_est.k + bootResu$Beta.LPR
          boot_CI.k <- boot_CI.k + bootResu$interval.LPR
        }
      }
      rm(x, y)
      gc()


      results$reg_res <- bootResu


      fin_ref_taxon_name <- originRefTaxNam
      nTestcov <- length(testCovInOrder)

      boot_est <- boot_est.k / nRuns
      boot_CI <- boot_CI.k / nRuns
      ref_taxon_name <- originRefTaxNam
      p_value_save_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      est_save_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_up_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_low_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      se_mat <-
        matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))

      calculate_se <- function(x) {
        se_est <- abs(x[2] - x[1]) / (2 * qnorm(1 - bootLassoAlpha / 2))
        return(se_est)
      }

      for (ii in seq_len(nTestcov)) {
        se_est <-
          apply(boot_CI[, seq(ii + 1, ncol(boot_CI), nPredics + 1)], 2, calculate_se)
        p_value_unadj <- numeric(length(se_est))
        boot_est_par <-
          boot_est[seq(ii + 1, ncol(boot_CI), nPredics + 1)]
        boot_est_CI_low <-
          boot_CI[, seq(ii + 1, ncol(boot_CI), nPredics + 1)][1, ]
        boot_est_CI_up <-
          boot_CI[, seq(ii + 1, ncol(boot_CI), nPredics + 1)][2, ]

        for (j in seq_len(length(boot_est_par))) {
          if (se_est[j] == 0) {
            p_value_unadj[j] <- 1
          } else {
            p_value_unadj[j] <-
              (1 - pnorm(abs(boot_est_par[j] / se_est[j]))) * 2
          }
        }
        p_value_adj <- p.adjust(p_value_unadj, adjust_method)
        p_value_save_mat[ii, ] <- p_value_unadj
        est_save_mat[ii, ] <- boot_est_par
        CI_low_mat[ii, ] <- boot_est_CI_low
        CI_up_mat[ii, ] <- boot_est_CI_up
        se_mat[ii, ] <- se_est
      }
    }
    return(
      list(
        p_value_save_mat = p_value_save_mat,
        est_save_mat = est_save_mat,
        CI_low_mat = CI_low_mat,
        CI_up_mat = CI_up_mat,
        se_mat = se_mat
      )
    )
  }
  parallel::stopCluster(cl)

  p_value_save_mat <-
    DescTools::DoCall(cbind, lapply(phase2res, function(x) {
      x[["p_value_save_mat"]]
    }))
  p_value_save_mat <-
    p_value_save_mat[, order(taxa_sepname_numeric, decreasing = FALSE), drop = FALSE]

  est_save_mat <-
    DescTools::DoCall(cbind, lapply(phase2res, function(x) {
      x[["est_save_mat"]]
    }))
  est_save_mat <-
    est_save_mat[, order(taxa_sepname_numeric, decreasing = FALSE), drop = FALSE]

  CI_low_mat <-
    DescTools::DoCall(cbind, lapply(phase2res, function(x) {
      x[["CI_low_mat"]]
    }))
  CI_low_mat <-
    CI_low_mat[, order(taxa_sepname_numeric, decreasing = FALSE), drop = FALSE]

  CI_up_mat <-
    DescTools::DoCall(cbind, lapply(phase2res, function(x) {
      x[["CI_up_mat"]]
    }))
  CI_up_mat <-
    CI_up_mat[, order(taxa_sepname_numeric, decreasing = FALSE), drop = FALSE]

  se_mat <- DescTools::DoCall(cbind, lapply(phase2res, function(x) {
    x[["se_mat"]]
  }))
  se_mat <- se_mat[, order(taxa_sepname_numeric, decreasing = FALSE), drop = FALSE]



  rownames(p_value_save_mat) <- testCovInOrder
  rownames(est_save_mat) <- testCovInOrder
  rownames(CI_low_mat) <- testCovInOrder
  rownames(CI_up_mat) <- testCovInOrder
  rownames(se_mat) <- testCovInOrder

  ref_taxon_name <- originRefTaxNam
  colname_use <- microbName[microbName != ref_taxon_name]
  colnames(p_value_save_mat) <- colname_use
  colnames(est_save_mat) <- colname_use
  colnames(CI_low_mat) <- colname_use
  colnames(CI_up_mat) <- colname_use
  colnames(se_mat) <- colname_use

  if (length(unbalanceTaxa_ori_name) > 0) {
    est_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
      NA
    CI_low_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
      NA
    CI_up_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
      NA
    se_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
      NA
    p_value_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
      NA
  }

  sig_ind <-
    which(p_value_save_mat < fwerRate,
      arr.ind = TRUE,
      useNames = FALSE
    )
  est_sig <- est_save_mat[sig_ind]
  CI_low_sig <- CI_low_mat[sig_ind]
  CI_up_sig <- CI_up_mat[sig_ind]
  p_adj_sig <- p_value_save_mat[sig_ind]
  se_sig <- se_mat[sig_ind]

  cov_sig_index <- sort(unique(sig_ind[, 1]))
  sig_list_each <- list()
  if (length(cov_sig_index) > 0) {
    for (iii in seq_len(length(cov_sig_index))) {
      sig_loc <- which(sig_ind[, 1] == cov_sig_index[iii])
      est_spe_cov <- est_sig[sig_loc]
      CI_low_spe_cov <- CI_low_sig[sig_loc]
      CI_up_spe_cov <- CI_up_sig[sig_loc]
      p_adj_spe_cov <- p_adj_sig[sig_loc]
      se_spe_cov <- se_sig[sig_loc]
      cov_sig_mat <- matrix(nrow = length(sig_loc), ncol = 5)
      colnames(cov_sig_mat) <-
        c("estimate", "SE est", "CI low", "CI up", "adj p-value")
      cov_sig_mat[, 1] <- est_spe_cov
      cov_sig_mat[, 2] <- se_spe_cov
      cov_sig_mat[, 3] <- CI_low_spe_cov
      cov_sig_mat[, 4] <- CI_up_spe_cov
      cov_sig_mat[, 5] <- p_adj_spe_cov
      rownames(cov_sig_mat) <- colname_use[sig_ind[sig_loc, 2]]
      sig_list_each[[testCovInOrder[cov_sig_index[iii]]]] <-
        cov_sig_mat
    }
  }
  all_cov_list <- list()


  all_cov_list$est_save_mat <- est_save_mat
  all_cov_list$p_value_save_mat <- p_value_save_mat
  all_cov_list$CI_low_mat <- CI_low_mat
  all_cov_list$CI_up_mat <- CI_up_mat
  all_cov_list$se_mat <- se_mat

  results$sig_list_each <- sig_list_each
  results$all_cov_list <- all_cov_list

  return(results)
}
