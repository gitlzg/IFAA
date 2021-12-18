
Regulariz=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  microbName,
  nRef,
  nRefMaxForEsti,
  refTaxa,
  # refTaxa_P2,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  fwerRate,
  bootB,
  standardize,
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

  # load abundance data info
  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)

  # results$x1permut=x1permut

  regul.start.time = proc.time()[3]
  message("Start Phase 1 association identification")

  selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                           testCovInOrder=testCovInOrder,
                           testCovInNewNam=testCovInNewNam,nRef=nRef,
                           paraJobs=paraJobs,
                           refTaxa=refTaxa,fwerRate=fwerRate,
                           standardize=standardize,
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
  while_loop_ind<-F
  loop_num<-0
  while (while_loop_ind==F) {
    loop_num<-loop_num+1
    if (loop_num>=10) {
      break
    }
    if (length(refTaxa)<nRef_smaller) {
      refTaxa_smaller<-c(refTaxa,head(colnames(selectRegroup$selecCountOverall)[order(selectRegroup$selecCountOverall)],
                                         n=nRef_smaller-length(refTaxa)))
    } else {
      refTaxa_smaller<-refTaxa
    }

    fin_ref_1<-selectRegroup$finalIndpRefTax
    ref_taxa_1<-selectRegroup$refTaxa
    selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                             testCovInOrder=testCovInOrder,
                             testCovInNewNam=testCovInNewNam,nRef=nRef_smaller,
                             paraJobs=paraJobs,
                             refTaxa=refTaxa_smaller,fwerRate=fwerRate,
                             standardize=standardize,
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
    # cat(fin_ref_2,"\n",ref_taxa_2,"\n")
    while_loop_ind<-identical(fin_ref_1,fin_ref_2) || identical(ref_taxa_1,ref_taxa_2)
  }
  results$selecCountOverall=selectRegroup$selecCountOverall
  colnames(results$selecCountOverall)<-microbName
  results$selecCountMatIndv=selectRegroup$selecCountMatIndv
  finalIndpRefTax=microbName[taxaNames%in%(selectRegroup$finalIndpRefTax)]
  results$finalRefTaxonQualified=selectRegroup$refTaxonQualified
  results$goodIndpRefTaxLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxLeastCount)]
  results$goodIndpRefTaxWithCount=selectRegroup$goodIndpRefTaxWithCount
  names(results$goodIndpRefTaxWithCount)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxWithCount)]

  results$goodIndpRefTaxWithEst=selectRegroup$goodIndpRefTaxWithEst
  names(results$goodIndpRefTaxWithEst)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxWithEst)]

  results$goodRefTaxaCandi=microbName[taxaNames%in%(selectRegroup$goodRefTaxaCandi)]
  results$randomRefTaxa=microbName[taxaNames%in%(selectRegroup$refTaxa)]
  goodIndpRefTax.ascend=sort(results$goodIndpRefTaxWithCount)
  goodIndpRefTaxNam=names(goodIndpRefTax.ascend)
  rm(selectRegroup)





  MCPExecuTime = (proc.time()[3] - regul.start.time)/60
  results$MCPExecuTime=MCPExecuTime
  message("Phase 1 Associaiton identification is done and used ", round(MCPExecuTime,2)," minutes")

  results$finalizedBootRefTaxon=finalIndpRefTax

  startT=proc.time()[3]
  message("Start Phase 2 parameter estimation")

  unestimableTaxa=c()
  qualifyData=data
  if(length(binaryInd)>0){
    qualifyData=data[rowSums(data[,results$taxaNames]>0)>=2,,drop=FALSE]
    firstBinPredNam=paste0(covsPrefix,binaryInd)
    allBinPred=paste0(covsPrefix,binaryInd:nPredics)
    nBinPred=length(allBinPred)

    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax=taxaNames[!taxaNames%in%(results$finalizedBootRefTaxon)]
    unbalanceTaxa=c()
    for(i in AllTaxaNamesNoRefTax){
      for(j in allBinPred){
        twoColumns.ij=qualifyData[,c(i,j)]
        nNonZero=sum(twoColumns.ij[,1]>0)
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        if(sumOfBin%in%c(0,1,(nNonZero-1),nNonZero)){
          unbalanceTaxa=c(unbalanceTaxa,i)
        }
      }
    }
    unestimableTaxa=unique(unbalanceTaxa)
    rm(allBinPred,unbalanceTaxa)
  }

  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(qualifyData[,taxaNames])==0)
  rm(qualifyData)
  unestimableTaxa=unique(c(unestimableTaxa,taxaNames[TaxaNoReads]))
  results$unEstTaxaPos=which(taxaNames%in%unestimableTaxa)
  rm(TaxaNoReads,unestimableTaxa)



  allRefTaxNam=unique(c(results$finalizedBootRefTaxon,goodIndpRefTaxNam))
  nGoodIndpRef=length(allRefTaxNam)
  results$allRefTaxNam=allRefTaxNam

  results$nRefUsedForEsti=min(nGoodIndpRef,nRefMaxForEsti)

  message("Final Reference Taxa are: ",allRefTaxNam[seq(results$nRefUsedForEsti)])

  results$estiList=list()
  for(iii in 1:(results$nRefUsedForEsti)){
    message("Start estimation for the ", iii,"th final reference taxon: ",allRefTaxNam[iii])
    time11=proc.time()[3]
    originTaxNam=allRefTaxNam[iii]
    newRefTaxNam=taxaNames[microbName%in%originTaxNam]
    # bootLassoAlpha_bon<- 1-fwerRate^(1/(nTaxa-1))
    bootLassoAlpha_bon<-0.05
    results$estiList[[originTaxNam]]=bootResuHDCI(data=data,
                                                       refTaxa=newRefTaxNam,
                                                       originRefTaxNam=originTaxNam,
                                                       bootB=bootB,bootLassoAlpha=bootLassoAlpha_bon,
                                                       binPredInd=binaryInd,covsPrefix=covsPrefix,
                                                       Mprefix=Mprefix,
                                                  testCovInOrder=testCovInOrder,
                                                  adjust_method=adjust_method,
                                                  microbName=microbName,
                                                  fwerRate=fwerRate,
                                                      paraJobs=paraJobs,
                                                       standardize=standardize,
                                                       seed=seed)
    time12=proc.time()[3]
    message("Estimation done for the ", iii,"th final reference taxon: ",allRefTaxNam[iii],
        " and it took ",round((time12-time11)/60,3)," minutes")
  }
  # estiResults=results$estiList[[results$finalizedBootRefTaxon]]

  endT=proc.time()[3]

  message("Phase 2 parameter estimation done and took ",round((endT-startT)/60,3)," minutes.")

  #### Analyze code ####
  # fin_ref_taxon_name<-names(results$estiList)
  # all_cov_sig_list<-list()
  # all_cov_list<-list()
  # nTestcov<-length(testCovInOrder)
  # for (i in 1:results$nRefUsedForEsti) {
  #   boot_est<-results$estiList[[i]]$bootResu$Beta.LPR
  #   boot_CI<-results$estiList[[i]]$bootResu$interval.LPR
  #   ref_taxon_name<-names(results$estiList)[i]
  #   p_value_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
  #   est_save_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
  #   CI_up_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
  #   CI_low_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
  #   se_mat<-matrix(nrow = nTestcov,ncol = (nTaxa-1))
  #
  #   for (ii in 1:nTestcov) {
  #     se_est<-apply(boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)],2,calculate_se)
  #     p_value_unadj<-numeric(length(se_est))
  #     boot_est_par<-boot_est[seq(ii+1,ncol(boot_CI),nPredics+1)]
  #     boot_est_CI_low<-boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)][1,]
  #     boot_est_CI_up<-boot_CI[,seq(ii+1,ncol(boot_CI),nPredics+1)][2,]
  #
  #     for (j in 1:length(boot_est_par)) {
  #       if (se_est[j]==0) {
  #         p_value_unadj[j]<-1
  #       } else {
  #         p_value_unadj[j]<-(1-pnorm(abs(boot_est_par[j]/se_est[j])))*2
  #       }
  #     }
  #     p_value_adj<-p.adjust(p_value_unadj,adjust_method)
  #     p_value_save_mat[ii,]<-p_value_adj
  #     est_save_mat[ii,]<-boot_est_par
  #     CI_low_mat[ii,]<-boot_est_CI_low
  #     CI_up_mat[ii,]<-boot_est_CI_up
  #     se_mat[ii,]<-se_est
  #
  #   }
  #   rownames(p_value_save_mat)<-testCovInOrder
  #   rownames(est_save_mat)<-testCovInOrder
  #   rownames(CI_low_mat)<-testCovInOrder
  #   rownames(CI_up_mat)<-testCovInOrder
  #   rownames(se_mat)<-testCovInOrder
  #
  #   colname_use<-microbName[microbName!=ref_taxon_name]
  #   colnames(p_value_save_mat)<-colname_use
  #   colnames(est_save_mat)<-colname_use
  #   colnames(CI_low_mat)<-colname_use
  #   colnames(CI_up_mat)<-colname_use
  #   colnames(se_mat)<-colname_use
  #
  #   sig_ind<-which(p_value_save_mat<fwerRate,arr.ind = T,useNames = F)
  #   est_sig<-est_save_mat[sig_ind]
  #   CI_low_sig<-CI_low_mat[sig_ind]
  #   CI_up_sig<-CI_up_mat[sig_ind]
  #   p_adj_sig<-p_value_save_mat[sig_ind]
  #   se_sig<-se_mat[sig_ind]
  #
  #   cov_sig_index<-sort(unique(sig_ind[,1]))
  #   sig_list_each<-list()
  #   if (length(cov_sig_index)>0) {
  #     for (iii in 1:length(cov_sig_index)) {
  #       sig_loc<-which(sig_ind[,1]==cov_sig_index[iii])
  #       est_spe_cov<-est_sig[sig_loc]
  #       CI_low_spe_cov<-CI_low_sig[sig_loc]
  #       CI_up_spe_cov<-CI_up_sig[sig_loc]
  #       p_adj_spe_cov<-p_adj_sig[sig_loc]
  #       se_spe_cov<-se_sig[sig_loc]
  #       cov_sig_mat<-matrix(nrow=length(sig_loc),ncol = 5)
  #       colnames(cov_sig_mat)<-c("estimate","SE est","CI low","CI up","adj p-value")
  #       cov_sig_mat[,1]<-est_spe_cov
  #       cov_sig_mat[,2]<-se_spe_cov
  #       cov_sig_mat[,3]<-CI_low_spe_cov
  #       cov_sig_mat[,4]<-CI_up_spe_cov
  #       cov_sig_mat[,5]<-p_adj_spe_cov
  #       rownames(cov_sig_mat)<-colname_use[sig_ind[sig_loc,2]]
  #       sig_list_each[[testCovInOrder[cov_sig_index[iii]]]]<-cov_sig_mat
  #     }
  #
  #
  #     all_cov_sig_list[[fin_ref_taxon_name[i]]]<-sig_list_each
  #   }
  #   all_cov_list[[fin_ref_taxon_name[i]]]$est_save_mat<-est_save_mat
  #   all_cov_list[[fin_ref_taxon_name[i]]]$p_value_save_mat<-p_value_save_mat
  #   all_cov_list[[fin_ref_taxon_name[i]]]$CI_low_mat<-CI_low_mat
  #   all_cov_list[[fin_ref_taxon_name[i]]]$CI_up_mat<-CI_up_mat
  #   all_cov_list[[fin_ref_taxon_name[i]]]$se_mat<-se_mat
  # }
  # results$all_cov_list<-all_cov_list
  # results$all_cov_sig_list<-all_cov_sig_list

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
  # testCovInOrder<-IFAAresul$testCov

  # subset(all_cov_list[[1]]$est_save_mat,select = -exclu_1)

  est_save_mat_mean<-(all_cov_list[[1]]$est_save_mat[,exclu_1,drop=F]+all_cov_list[[2]]$est_save_mat[,exclu_2,drop=F])/2
  se_mat_mean<-(all_cov_list[[1]]$se_mat[,exclu_1,drop=F]+all_cov_list[[2]]$se_mat[,exclu_2,drop=F])/2
  CI_low_mat_mean<-(all_cov_list[[1]]$CI_low_mat[,exclu_1,drop=F]+all_cov_list[[2]]$CI_low_mat[,exclu_2,drop=F])/2
  CI_up_mat_mean<-(all_cov_list[[1]]$CI_up_mat[,exclu_1,drop=F]+all_cov_list[[2]]$CI_up_mat[,exclu_2,drop=F])/2
  p_value_unadj_mean<-t(apply(est_save_mat_mean/se_mat_mean,1,function(x) (1-pnorm(abs(x)))*2))
  p_value_adj_mean<-t(apply(p_value_unadj_mean,1,function(x) p.adjust(x,method = adjust_method)))
  colname_use<-colnames(est_save_mat_mean)

  sig_ind<-which(p_value_adj_mean<fwerRate,arr.ind = T,useNames = F)
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

  results$sig_list_each_mean<-sig_list_each_mean
  results$all_cov_save_mean<-list(est_save_mat_mean,se_mat_mean,CI_low_mat_mean,CI_up_mat_mean,p_value_adj_mean)


  # results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
  # results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
  # results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")
  #
  # results$betaMat.LPR=as(matrix(estiResults$finalBetaEst.LPR,nrow=nPredics),"sparseMatrix")
  # results$CILowMat.LPR=as(matrix(estiResults$CIvecLow.LPR,nrow=nPredics),"sparseMatrix")
  # results$CIUpMat.LPR=as(matrix(estiResults$CIvecUp.LPR,nrow=nPredics),"sparseMatrix")
  # rm(estiResults)
  #
  # estByCovList=list()
  #
  # for(i in 1:nTestCov){
  #   sigTaxaPosition=which(results$selecIndvInOverall[i,]!=0)
  #   nrow=length(sigTaxaPosition)
  #   if(nrow==0)next
  #
  #   ncol=3
  #   estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)
  #
  #   for(j in 1:nrow){
  #     for(k in 1:ncol){
  #       if(k==1)estByCovMat[j,k]=results$betaMat.LPR[i,sigTaxaPosition[j]]
  #       if(k==2)estByCovMat[j,k]=results$CILowMat.LPR[i,sigTaxaPosition[j]]
  #       if(k==3)estByCovMat[j,k]=results$CIUpMat.LPR[i,sigTaxaPosition[j]]
  #     }
  #   }
  #
  #   rownames(estByCovMat)=microbName[results$selecIndvInOverall[i,]!=0]
  #   colnames(estByCovMat)=c("Beta.LPR","LowB95%CI.LPR","UpB95%CI.LPR")
  #
  #   rowsToKeep=which(estByCovMat[,1]!=0)
  #   estByCovMat=estByCovMat[rowsToKeep, , drop = FALSE]
  #
  #   estByCovList[[testCovInOrder[i]]]=estByCovMat
  #   rm(estByCovMat)
  # }
  #
  # if(length(estByCovList)==0){
  #   results$estByCovList="No significant assoication is identified. You may consider a higher FDR level."
  # }else{
  #   results$estByCovList=estByCovList
  # }
  # rm(estByCovList)
  #
  # SigCovByTaxaList=list()
  # for(i in 1:nTaxa){
  #   sigCov=which(results$selecIndvInOverall[,i]!=0)
  #   if(length(sigCov)==0)next
  #   SigCovByTaxaList[[microbName[i]]]=testCovInOrder[sigCov]
  # }
  # if(length(SigCovByTaxaList)==0){
  #   results$SigCovByTaxaList="No significant assoication is identified. You may consider a higher FDR level."
  # }else{
  #   results$SigCovByTaxaList=SigCovByTaxaList
  # }
  # rm(SigCovByTaxaList,microbName)
  #
  # results$reguMethod=reguMethod
  # results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  # return(results)

  rm(data)

  # return results
  # results$fwerRate=fwerRate
  # results$nPermu=nPermu
  results$nRef=nRef
  return(results)
}

#
