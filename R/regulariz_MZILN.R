Regulariz_MZILN=function(
  data,
  nRef,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  microbName,
  targetTaxa,
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
  allFunc=allFunc,
  seed
) {
  results=list()
  regul.start.time = proc.time()[3]

  nTestCov=length(testCovInd)

  dataSparsCheck(data=data,Mprefix=Mprefix)

  # load abundance data info
  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)

  refTaxa_reOrder=microbName[taxaNames%in%refTaxa]

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


  betaMatList=list()
  CILowMatList=list()
  CIUpMatList=list()
  results$estiList=list()
  for(iii in 1:nRef){
    time11=proc.time()[3]
    originTaxNam=refTaxa_reOrder[iii]
    newRefTaxNam<-refTaxa[iii]
    bootLassoAlpha_bon<-fdrRate
    results$estiList[[originTaxNam]]=bootResuHDCI(data=data,
                                                  refTaxa=newRefTaxNam,
                                                  originRefTaxNam=originTaxNam,
                                                  bootB=bootB,bootLassoAlpha=bootLassoAlpha_bon,
                                                  binPredInd=binaryInd,covsPrefix=covsPrefix,
                                                  Mprefix=Mprefix,
                                                  unbalanceTaxa_ori_name=unbalanceTaxa_ori_name,
                                                  unbalancePred_ori_name=unbalancePred_ori_name,
                                                  testCovInOrder=testCovInOrder,
                                                  adjust_method=adjust_method,
                                                  microbName=microbName,
                                                  fwerRate=fdrRate,
                                                  paraJobs=paraJobs,
                                                  seed=seed)

    time12=proc.time()[3]
    message("Estimation done for the ", iii,"th denominator taxon: ",refTaxa_reOrder[iii],
            " and it took ",round((time12-time11)/60,2)," minutes")


  }

  fin_ref_taxon_name<-names(results$estiList)
  all_cov_sig_list<-list()
  all_cov_list<-list()

  for (i in 1:length(fin_ref_taxon_name)) {
    all_cov_list[[fin_ref_taxon_name[i]]]<-results$estiList[[i]]$all_cov_list
    all_cov_sig_list[[fin_ref_taxon_name[i]]]<-results$estiList[[i]]$sig_list_each
  }

  tgtaxa_save_list<-list()

  if (length(targetTaxa)!=0) {
    for (i in 1:length(fin_ref_taxon_name)) {
      col_to_keep<-colnames(all_cov_list[[fin_ref_taxon_name[i]]]$est_save_mat) %in% targetTaxa
      est_to_keep<-all_cov_list[[fin_ref_taxon_name[i]]]$est_save_mat[,col_to_keep,drop=FALSE]
      se_to_keep<-all_cov_list[[fin_ref_taxon_name[i]]]$se_mat[,col_to_keep,drop=FALSE]
      CI_low_to_keep<-all_cov_list[[fin_ref_taxon_name[i]]]$CI_low_mat[,col_to_keep,drop=FALSE]
      CI_up_to_keep<-all_cov_list[[fin_ref_taxon_name[i]]]$CI_up_mat[,col_to_keep,drop=FALSE]
      p_value_to_keep<-all_cov_list[[fin_ref_taxon_name[i]]]$p_value_save_mat[,col_to_keep,drop=FALSE]
      ref_each_save_list<-list()
      for (j in 1:nrow(est_to_keep)) {
        res_show_mat<-cbind(est_to_keep[j,],se_to_keep[j,],CI_low_to_keep[j,],CI_up_to_keep[j,],p_value_to_keep[j,])
        colnames(res_show_mat)<-c("estimate","SE est","CI low","CI up","adj p-value")
        rownames(res_show_mat)<-targetTaxa
        ref_each_save_list[[rownames(est_to_keep)[j]]]<-res_show_mat
      }
      tgtaxa_save_list[[fin_ref_taxon_name[i]]]<-ref_each_save_list

    }
  }



  results$full_results<-all_cov_list
  results$sig_results<-all_cov_sig_list
  results$targettaxa_result_list=tgtaxa_save_list
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  rm(data)
  # return results
  return(results)
}
