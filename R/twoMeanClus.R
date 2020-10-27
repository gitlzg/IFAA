
twoMeanClus=function(
  matrix,
  nRef
){
  results=list()
  if(sum(matrix<0)>0) stop("Count matrix contains negative values")

  matrix=as.matrix(t(as.matrix(matrix)))

  n.col=ncol(matrix)
  n.row=nrow(matrix)
  finalValVec=rep(NA,n.col)
  twoMeanSucVec=rep(NA,n.col)
  finalSelec=matrix(NA,nrow=n.row,ncol=n.col)
  for (i in 1:n.col){
    vecI=matrix[,i]
    cluster2=try(kmeans(vecI,centers=c(0,nRef)),silent=TRUE)
    twoMeanSucVec[i]=1-inherits(cluster2,"try-error")

    if(inherits(cluster2,"try-error")){
      finalSelec[1,i]=1
      finalSelec[2:n.row,i]=0
    }
    else{
      hiClus=which(cluster2$centers==max(cluster2$centers))
      finalSelec[,i]=(cluster2$cluster==hiClus)+0
    }
    finalValVec[i]=min(vecI[(finalSelec[,i]==1)])
    rm(vecI)
  }

  results$finalSelec=t(as.matrix(finalSelec))

  rm(finalSelec)
  results$twoMeanCut=finalValVec
  rm(finalValVec)
  results$twoMeanSucVec=twoMeanSucVec
  rm(twoMeanSucVec)
  return(results)
}
