
##' @export



dataRecovTrans=function(
  data,
  ref,
  Mprefix,
  covsPrefix,
  xOnly=F,
  yOnly=F
){
  results=list()
  
  # load A and log-ratio transformed RA
  data.and.init=AIcalcu(data=data,ref=ref,Mprefix=Mprefix,
                        covsPrefix=covsPrefix)
  rm(data)
  
  taxaNames=data.and.init$taxaNames
  A=data.and.init$A
  logRatiow=data.and.init$logRatiow
  nSub=data.and.init$nSub
  nTaxa=data.and.init$nTaxa
  xData=data.and.init$xData
  nPredics=data.and.init$nPredics
  twoList=data.and.init$twoList
  lLast=data.and.init$lLast
  L=data.and.init$l
  lengthTwoList=data.and.init$lengthTwoList
  rm(data.and.init) 
  
  nNorm=nTaxa-1
  xDimension=nPredics+1 # predictors+intercept
  
  # create omegaRoot
  omegaRoot=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if(lLast[i]==nTaxa){
      #omega[[i]]=Diagonal(L[i]-1)
      omegaRoot[[i]]=Diagonal(L[i]-1)
    } else{ if(L[i]==2){
      #omega[[i]]=1/2
      omegaRoot[[i]]=sqrt(0.5)
    }else{
      dim=L[i]-1
      a=(1+(dim-2)/2)/(1/2*(1+(dim-1)/2))
      b=-1/(1+(dim-1)/2)
      #omega[[i]]=1/2*((a-b)*Diagonal(dim)+b*matrix(1,nrow=dim,ncol=dim))
      
      # calculate the square root of omega assuming it is exchangeable
      aStar=dim^2/((dim-1)^2)
      bStar=b*(dim-2)/(dim-1)-a*(dim^2-2*dim+2)/((dim-1)^2)
      cStar=(0.5*b-0.5*a*(dim-2)/(dim-1))^2
      cSquare=(-bStar+sqrt(bStar^2-4*aStar*cStar))/(2*aStar)
      if(cSquare<0) {stop("no solution for square root of omega")}
      d=sqrt((0.5*a-cSquare)/(dim-1))
      if (is.na(d)){stop("no solution for off-diagnal elements for square root of omega")}
      c=(0.5*b-(dim-2)*(d^2))/(2*d)
      omegaRoot[[i]]=-((c-d)*Diagonal(dim)+d*matrix(1,nrow=dim,ncol=dim))
    }
    }
  }
  rm(L,lLast)
  
  if(xOnly){
    # create X_i in the regression equation using Kronecker product
    xDataWithInter=as.matrix(cbind(rep(1,nSub),xData))
    rm(xData)
    
    xInRegres=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      xInRegres[[i]]=as(t(as.matrix(kronecker(Diagonal(nNorm),xDataWithInter[i,]))),"sparseMatrix")
    } 
    rm(xDataWithInter)
    
    xDataTilda=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      xDataTilda[[i]]=omegaRoot[[i]]%*%A[[i]]%*%xInRegres[[i]]
    }
    rm(omegaRoot,logRatiow,xInRegres)
    
    # stack xTilda
    for (j in 1:lengthTwoList){
      i=twoList[j]
      if (j==1) {xTildalong=xDataTilda[[i]]
      } else {
        xTildalong=rbind(xTildalong,xDataTilda[[i]])
      }
    }
    
    results$xTildalong=xTildalong
    rm(xTildalong)
    return(results)
  }
  
  if(yOnly){
    Utilda=list()
    for (j in 1:lengthTwoList){
      i=twoList[j]
      Utilda[[i]]=omegaRoot[[i]]%*%logRatiow[[i]]
    }
    rm(omegaRoot,logRatiow)
    
    # stack Utilda
    for (j in 1:lengthTwoList){
      i=twoList[j]
      if (j==1) {UtildaLong=Utilda[[i]]
      } else {
        UtildaLong=rbind(UtildaLong,Utilda[[i]])
      }
    }
    
    rm(Utilda)
    results$UtildaLong=UtildaLong
    rm(UtildaLong)
    return(results)
  }
  
  # create X_i in the regression equation using Kronecker product
  xDataWithInter=data.matrix(cbind(rep(1,nSub),xData))
  colnames(xDataWithInter)[1]="Inter"
  rm(xData)
  
  xInRegres=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    xInRegres[[i]]=as(t(as.matrix(kronecker(Diagonal(nNorm),xDataWithInter[i,]))),"sparseMatrix")
  } 
  
  rm(xDataWithInter)
  
  # transfer U and xInRegres using omega root
  xDataTilda=list()
  Utilda=list()
  
  for (j in 1:lengthTwoList){
    i=twoList[j]
    Utilda[[i]]=omegaRoot[[i]]%*%logRatiow[[i]]
    xDataTilda[[i]]=omegaRoot[[i]]%*%A[[i]]%*%xInRegres[[i]]
  }
  rm(omegaRoot,logRatiow,xInRegres)
  
  # stack Utilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {UtildaLong=Utilda[[i]]
    } else {
      UtildaLong=rbind(UtildaLong,Utilda[[i]])
    }
  }
  rm(Utilda)
  
  # stack xTilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {xTildalong=xDataTilda[[i]]
    } else {
      xTildalong=rbind(xTildalong,xDataTilda[[i]])
    }
  }
  rm(xDataTilda,lengthTwoList)
  
  # return objects 
  results$UtildaLong=UtildaLong
  rm(UtildaLong)
  results$xTildalong=xTildalong
  rm(xTildalong)
  
  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}
