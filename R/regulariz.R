
Regulariz=function(
  data,
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
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  adjust_method,
  seed
){
  results=list()
  regul.start.time = proc.time()[3]

  nTestCov=length(testCovInd)

  dataSparsCheck(data=data,Mprefix=Mprefix)


  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)


  regul.start.time = proc.time()[3]
  message("Start Phase 1 analysis")

  selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                           testCovInOrder=testCovInOrder,
                           testCovInNewNam=testCovInNewNam,nRef=nRef,
                           paraJobs=paraJobs,
                           refTaxa=refTaxa,
                           sequentialRun=sequentialRun,
                           allFunc=allFunc,
                           refReadsThresh=refReadsThresh,
                           SDThresh=SDThresh,
                           SDquantilThresh=SDquantilThresh,
                           balanceCut=balanceCut,
                           Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binaryInd,
                           adjust_method=adjust_method,
                           seed=seed)
  nRef_smaller<-max(2,ceiling(nRef/2))
  while_loop_ind<-FALSE
  loop_num<-0
  message("33 percent of phase 1 analysis has been done")
  while (while_loop_ind==FALSE) {
    if (loop_num>=2) {
      break
    }
    loop_num<-loop_num+1

    refTaxa_smaller<-head(names(selectRegroup$goodIndpRefTaxWithCount),
                                    n=nRef_smaller)


    fin_ref_1<-selectRegroup$finalIndpRefTax
    ref_taxa_1<-selectRegroup$refTaxa
    selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                             testCovInOrder=testCovInOrder,
                             testCovInNewNam=testCovInNewNam,nRef=nRef_smaller,
                             paraJobs=paraJobs,
                             refTaxa=refTaxa_smaller,
                             sequentialRun=sequentialRun,
                             allFunc=allFunc,
                             refReadsThresh=refReadsThresh,
                             SDThresh=SDThresh,
                             SDquantilThresh=SDquantilThresh,
                             balanceCut=balanceCut,
                             Mprefix=Mprefix,
                             covsPrefix=covsPrefix,
                             binPredInd=binaryInd,
                             adjust_method=adjust_method,
                             seed=seed)


    fin_ref_2<-selectRegroup$finalIndpRefTax
    ref_taxa_2<-selectRegroup$refTaxa
    while_loop_ind<-identical(fin_ref_1,fin_ref_2) || identical(ref_taxa_1,ref_taxa_2)

    if(while_loop_ind==FALSE){
      message(round(100* (loop_num+1)/3,0), " percent of phase 1 analysis has been done")
    }
    if(while_loop_ind==TRUE){
      message("100 percent of phase 1 analysis has been done")
    }
  }
  results$selecCountOverall=selectRegroup$selecCountOverall
  colnames(results$selecCountOverall)<-microbName
  results$selecCountMatIndv=selectRegroup$selecCountMatIndv
  finalIndpRefTax=microbName[taxaNames%in%(selectRegroup$finalIndpRefTax)]
  results$finalRefTaxonQualified=selectRegroup$refTaxonQualified
  results$goodIndpRefTaxLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxLeastCount)]
  results$goodIndpRefTaxWithCount=selectRegroup$goodIndpRefTaxWithCount
  names(results$goodIndpRefTaxWithCount)=microbName[sapply(names(selectRegroup$goodIndpRefTaxWithCount),function(x) which(taxaNames%in%x))]

  results$goodIndpRefTaxWithEst=selectRegroup$goodIndpRefTaxWithEst
  names(results$goodIndpRefTaxWithEst)=microbName[sapply(names(selectRegroup$goodIndpRefTaxWithEst),function(x) which(taxaNames%in%x))]


  results$goodRefTaxaCandi=microbName[taxaNames%in%(selectRegroup$goodRefTaxaCandi)]
  results$randomRefTaxa=microbName[taxaNames%in%(selectRegroup$refTaxa)]
  goodIndpRefTax.ascend=sort(results$goodIndpRefTaxWithCount)
  goodIndpRefTaxNam=names(goodIndpRefTax.ascend)
  rm(selectRegroup)





  MCPExecuTime = (proc.time()[3] - regul.start.time)/60
  results$MCPExecuTime=MCPExecuTime
  message("Phase 1 analysis used ", round(MCPExecuTime,2)," minutes")

  results$finalizedBootRefTaxon=finalIndpRefTax

  startT=proc.time()[3]
  message("Start Phase 2 parameter estimation")

  qualifyData=data


  if(length(binaryInd_test)>0){
    qualifyData=data[rowSums(data[,taxaNames]>0)>=2,,drop=FALSE]
    allBinPred=paste0(covsPrefix,binaryInd_test)
    nBinPred=length(allBinPred)

    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax=taxaNames[!microbName%in%(results$finalizedBootRefTaxon)]
    unbalanceTaxa=c()
    unbalancePred<-c()
    for(i in AllTaxaNamesNoRefTax){
      for(j in allBinPred){
        twoColumns.ij=qualifyData[,c(i,j)]
        nNonZero=sum(twoColumns.ij[,1]>0)
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        if(sumOfBin%in%c(0,1,(nNonZero-1),nNonZero)){
          unbalanceTaxa=c(unbalanceTaxa,i)
          unbalancePred<-c(unbalancePred,j)
        }
      }
    }
    if (length(unbalanceTaxa)>0) {
      unbalanceTaxa_ori_name<-microbName[sapply(unbalanceTaxa,function(x) which(taxaNames%in%x))]
      unbalancePred_ori_name<-testCovInOrder[sapply(unbalancePred,function(x) which(testCovInNewNam%in%x))]
    } else {
      unbalanceTaxa_ori_name<-NULL
      unbalancePred_ori_name<-NULL
    }
    # unestimableTaxa=unique(unbalanceTaxa)
    rm(allBinPred)
  } else {
    unbalanceTaxa_ori_name<-NULL
    unbalancePred_ori_name<-NULL
  }

  # check zero taxa and subjects with zero taxa reads




  allRefTaxNam=unique(c(results$finalizedBootRefTaxon,goodIndpRefTaxNam))
  nGoodIndpRef=length(allRefTaxNam)
  results$allRefTaxNam=allRefTaxNam

  results$nRefUsedForEsti=min(nGoodIndpRef,nRefMaxForEsti)

  results$estiList=list()
  for(iii in 1:(results$nRefUsedForEsti)){
    message("Start estimation for the ", iii,"th final reference taxon")
    time11=proc.time()[3]
    originTaxNam=allRefTaxNam[iii]
    newRefTaxNam=taxaNames[microbName%in%originTaxNam]
    results$estiList[[originTaxNam]]=bootResuHDCI(data=data,
                                                  refTaxa=newRefTaxNam,
                                                  originRefTaxNam=originTaxNam,
                                                  bootB=bootB,
                                                  binPredInd=binaryInd,covsPrefix=covsPrefix,
                                                  Mprefix=Mprefix,
                                                  unbalanceTaxa_ori_name=unbalanceTaxa_ori_name,
                                                  unbalancePred_ori_name=unbalancePred_ori_name,
                                                  testCovInOrder=testCovInOrder,
                                                  adjust_method=adjust_method,
                                                  microbName=microbName,
                                                  fwerRate=fwerRate,
                                                  paraJobs=paraJobs,
                                                  seed=seed)
    time12=proc.time()[3]
    message("Estimation done for the ", iii,"th final reference taxon and it took ",round((time12-time11)/60,3)," minutes")
  }

  endT=proc.time()[3]

  message("Phase 2 parameter estimation done and took ",round((endT-startT)/60,3)," minutes.")

  ### calculate mean ###
  fin_ref_taxon_name<-names(results$estiList)
  all_cov_sig_list<-list()
  all_cov_list<-list()

  for (i in 1:length(fin_ref_taxon_name)) {
    all_cov_list[[fin_ref_taxon_name[i]]]<-results$estiList[[i]]$all_cov_list
    all_cov_sig_list[[fin_ref_taxon_name[i]]]<-results$estiList[[i]]$sig_list_each
  }
  results$all_cov_list<-all_cov_list
  results$all_cov_sig_list<-all_cov_sig_list


  ref_taxon_name<-names(all_cov_list)

  exclu_1<-!colnames(all_cov_list[[1]]$est_save_mat)%in%ref_taxon_name[2]
  exclu_2<-!colnames(all_cov_list[[2]]$est_save_mat)%in%ref_taxon_name[1]


  est_save_mat_mean<-(all_cov_list[[1]]$est_save_mat[,exclu_1,drop=FALSE]+all_cov_list[[2]]$est_save_mat[,exclu_2,drop=FALSE])/2
  se_mat_mean<-(all_cov_list[[1]]$se_mat[,exclu_1,drop=FALSE]+all_cov_list[[2]]$se_mat[,exclu_2,drop=FALSE])/2
  CI_low_mat_mean<-(all_cov_list[[1]]$CI_low_mat[,exclu_1,drop=FALSE]+all_cov_list[[2]]$CI_low_mat[,exclu_2,drop=FALSE])/2
  CI_up_mat_mean<-(all_cov_list[[1]]$CI_up_mat[,exclu_1,drop=FALSE]+all_cov_list[[2]]$CI_up_mat[,exclu_2,drop=FALSE])/2


  p_value_unadj_mean<-t(apply(est_save_mat_mean/se_mat_mean,1,function(x) (1-pnorm(abs(x)))*2))
  p_value_adj_mean<-t(apply(p_value_unadj_mean,1,function(x) p.adjust(x,method = adjust_method)))
  colname_use<-colnames(est_save_mat_mean)

  sig_ind<-which(p_value_adj_mean<=fwerRate,arr.ind = TRUE,useNames = FALSE)
  est_sig<-est_save_mat_mean[sig_ind]
  CI_low_sig<-CI_low_mat_mean[sig_ind]
  CI_up_sig<-CI_up_mat_mean[sig_ind]
  p_adj_sig<-p_value_adj_mean[sig_ind]
  se_sig<-se_mat_mean[sig_ind]

  cov_sig_index<-sort(unique(sig_ind[,1]))
  sig_list_each_mean<-list()
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
      sig_list_each_mean[[testCovInOrder[cov_sig_index[iii]]]]<-cov_sig_mat
    }
  }

  results$sig_results<-sig_list_each_mean
  full_results<-list()
  for (j in testCovInOrder) {
    est_res_save_all<-cbind(est_save_mat_mean[j,],se_mat_mean[j,],CI_low_mat_mean[j,],
                            CI_up_mat_mean[j,],p_value_adj_mean[j,])
    colnames(est_res_save_all)<-c("estimate","SE est","CI low","CI up","adj p-value")
    full_results[[j]]<-est_res_save_all
  }
  results$full_results<-full_results

  results$nTaxa=nTaxa
  results$nPredics=nPredics

  rm(data)

  # return results

  results$nRef=nRef
  return(results)
}

#
