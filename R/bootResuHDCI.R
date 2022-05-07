#---------------------------------------------------------------------------------------
## bootstrap results function for Lasso OLS HDCI
#---------------------------------------------------------------------------------------

bootResuHDCI=function(
  data,
  refTaxa,
  originRefTaxNam,
  maxDimension=434*5*10^5,
  bootB,
  bootLassoAlpha=0.05,
  binPredInd,
  covsPrefix,
  Mprefix,
  unbalanceTaxa_ori_name,
  unbalancePred_ori_name,
  testCovInOrder,
  adjust_method,
  microbName,
  fwerRate,
  paraJobs,
  seed
){

  results=list()
  # load data info
  basicInfo=dataInfo(data=data,binPredInd=binPredInd,
                     covsPrefix=covsPrefix,Mprefix=Mprefix)
  taxaNames=basicInfo$taxaNames
  ii=which(basicInfo$taxaNames%in%refTaxa)

  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()

  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa

  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()

  # inital Lasso OLS estimate
  dataForEst=dataRecovTrans(data=data,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix)

  x=as(dataForEst$xTildalong,"sparseMatrix")
  y=dataForEst$UtildaLong
  rm(dataForEst)

  xCol=ncol(x)

  maxSubSamplSiz=min(50000,floor(maxDimension/xCol))
  nToSamplFrom=length(y)
  subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
  if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

  nRuns=(ceiling(subSamplK/3))


  if (dim(x)[1]>(dim(x)[2])) {
    for(k in 1:nRuns){
      rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
      xSub=as((x[rowToKeep,]),"sparseMatrix")
      ySub=(y[rowToKeep])


      lm_res<-lm(as.vector(ySub)~as.matrix(xSub)-1)



      rm(xSub,ySub)

      full_name_coef<-names(lm_res$coefficients)
      valid_coef<-summary(lm_res)$coefficients
      bootResu<-matrix(nrow = length(full_name_coef),ncol = 4)
      rownames(bootResu)<-full_name_coef
      bootResu[rownames(bootResu)%in%rownames(valid_coef),]<-valid_coef
      if (k==1) {
        bootResu_k<-bootResu
      } else {
        bootResu_k=bootResu_k+bootResu
      }
      gc()

    }
    rm(x,y)
    gc()


    fin_ref_taxon_name<-originRefTaxNam
    all_cov_list<-list()
    nTestcov<-length(testCovInOrder)


    boot_est<-bootResu_k[,1]/nRuns
    se_est_all<-bootResu_k[,2]/nRuns


    boot_CI<-confint(lm_res)


    ref_taxon_name<-originRefTaxNam
    p_value_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    est_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    CI_up_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    CI_low_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    se_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))

    for (ii in 1:nTestcov) {
      se_est<-se_est_all[seq(ii+1,length(full_name_coef),nPredics+1)]
      boot_est_par<-boot_est[seq(ii+1,length(full_name_coef),nPredics+1)]

      p_value_unadj<-(1-pnorm(abs(boot_est_par/se_est)))*2

      boot_est_CI_low<-boot_est_par-1.96*se_est
      boot_est_CI_up<-boot_est_par+1.96*se_est

      p_value_adj<-p.adjust(p_value_unadj,adjust_method)
      p_value_save_mat[ii,]<-p_value_adj
      est_save_mat[ii,]<-boot_est_par
      CI_low_mat[ii,]<-boot_est_CI_low
      CI_up_mat[ii,]<-boot_est_CI_up
      se_mat[ii,]<-se_est

    }
  } else {
    for(k in 1:nRuns){
      rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
      xSub=as((x[rowToKeep,]),"sparseMatrix")
      ySub=as((y[rowToKeep]),"sparseVector")


      c3 <- parallel::makeCluster(paraJobs)
      doParallel::registerDoParallel(c3)

      if(length(seed)>0){
        set.seed(as.numeric(seed)+10^2)
        parallel::clusterSetRNGStream(cl=c3,(as.numeric(seed)+10^3))
      }

      bootResu=runBootLassoHDCI(x=as.matrix(xSub),y=as.vector(ySub),bootB=bootB,
                        paraJobs=paraJobs,bootLassoAlpha=bootLassoAlpha,seed=seed)

      parallel::stopCluster(c3)

      rm(xSub,ySub)
      gc()
      if (k==1) {
        boot_est.k<-bootResu$Beta.LPR
        boot_CI.k<-bootResu$interval.LPR
      } else {
        boot_est.k<-boot_est.k+bootResu$Beta.LPR
        boot_CI.k<-boot_CI.k+bootResu$interval.LPR
      }
    }
    rm(x,y)
    gc()


    results$reg_res<-bootResu


    fin_ref_taxon_name<-originRefTaxNam
    all_cov_list<-list()
    nTestcov<-length(testCovInOrder)

    boot_est<-boot_est.k/nRuns
    boot_CI<-boot_CI.k/nRuns
    ref_taxon_name<-originRefTaxNam
    p_value_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    est_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    CI_up_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    CI_low_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
    se_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))

    calculate_se<-function(x) {
      se_est<-abs(x[2]-x[1])/(2*qnorm(1-bootLassoAlpha/2))
      return(se_est)
      }

    for (ii in 1:nTestcov) {
      se_est<-apply(boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)],2,calculate_se)
      p_value_unadj<-numeric(length(se_est))
      boot_est_par<-boot_est[seq(ii+1,ncol(boot_CI),nPredics+1)]
      boot_est_CI_low<-boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)][1,]
      boot_est_CI_up<-boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)][2,]

      for (j in 1:length(boot_est_par)) {
        if (se_est[j]==0) {
          p_value_unadj[j]<-1
        } else {
          p_value_unadj[j]<-(1-pnorm(abs(boot_est_par[j]/se_est[j])))*2
        }
      }
      p_value_adj<-p.adjust(p_value_unadj,adjust_method)
      p_value_save_mat[ii,]<-p_value_adj
      est_save_mat[ii,]<-boot_est_par
      CI_low_mat[ii,]<-boot_est_CI_low
      CI_up_mat[ii,]<-boot_est_CI_up
      se_mat[ii,]<-se_est

    }
  }


  rownames(p_value_save_mat)<-testCovInOrder
  rownames(est_save_mat)<-testCovInOrder
  rownames(CI_low_mat)<-testCovInOrder
  rownames(CI_up_mat)<-testCovInOrder
  rownames(se_mat)<-testCovInOrder

  colname_use<-microbName[microbName!=ref_taxon_name]
  colnames(p_value_save_mat)<-colname_use
  colnames(est_save_mat)<-colname_use
  colnames(CI_low_mat)<-colname_use
  colnames(CI_up_mat)<-colname_use
  colnames(se_mat)<-colname_use

  if (length(unbalanceTaxa_ori_name)>0) {
    est_save_mat[cbind(unbalancePred_ori_name,unbalanceTaxa_ori_name)]<-NA
    CI_low_mat[cbind(unbalancePred_ori_name,unbalanceTaxa_ori_name)]<-NA
    CI_up_mat[cbind(unbalancePred_ori_name,unbalanceTaxa_ori_name)]<-NA
    se_mat[cbind(unbalancePred_ori_name,unbalanceTaxa_ori_name)]<-NA
    p_value_save_mat[cbind(unbalancePred_ori_name,unbalanceTaxa_ori_name)]<-NA
  }

  sig_ind<-which(p_value_save_mat<fwerRate,arr.ind = TRUE,useNames = FALSE)
  est_sig<-est_save_mat[sig_ind]
  CI_low_sig<-CI_low_mat[sig_ind]
  CI_up_sig<-CI_up_mat[sig_ind]
  p_adj_sig<-p_value_save_mat[sig_ind]
  se_sig<-se_mat[sig_ind]

  cov_sig_index<-sort(unique(sig_ind[,1]))
  sig_list_each<-list()
  if (length(cov_sig_index)>0) {
    for (iii in 1:length(cov_sig_index)) {
      sig_loc<-which(sig_ind[,1]==cov_sig_index[iii])
      est_spe_cov<-est_sig[sig_loc]
      CI_low_spe_cov<-CI_low_sig[sig_loc]
      CI_up_spe_cov<-CI_up_sig[sig_loc]
      p_adj_spe_cov<-p_adj_sig[sig_loc]
      se_spe_cov<-se_sig[sig_loc]
      cov_sig_mat<-matrix(nrow=length(sig_loc),ncol = 5)
      colnames(cov_sig_mat)<-c("estimate","SE est","CI low","CI up","adj p-value")
      cov_sig_mat[,1]<-est_spe_cov
      cov_sig_mat[,2]<-se_spe_cov
      cov_sig_mat[,3]<-CI_low_spe_cov
      cov_sig_mat[,4]<-CI_up_spe_cov
      cov_sig_mat[,5]<-p_adj_spe_cov
      rownames(cov_sig_mat)<-colname_use[sig_ind[sig_loc,2]]
      sig_list_each[[testCovInOrder[cov_sig_index[iii]]]]<-cov_sig_mat
    }
  }

  all_cov_list$est_save_mat<-est_save_mat
  all_cov_list$p_value_save_mat<-p_value_save_mat
  all_cov_list$CI_low_mat<-CI_low_mat
  all_cov_list$CI_up_mat<-CI_up_mat
  all_cov_list$se_mat<-se_mat

  results$sig_list_each<-sig_list_each
  results$all_cov_list<-all_cov_list

  return(results)
}
