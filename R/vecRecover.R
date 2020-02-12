##' @export


groupBetaToFullBeta=function(
  nTaxa,
  nPredics,
  unSelectList,
  newBetaNoInt
){

  results=list()
  unSelectList=unique(sort(unSelectList))
  nUnSelec=length(unSelectList)
  nAlphaSelec=nTaxa*nPredics
  nNewBetaNoInt=length(newBetaNoInt)

  if(nNewBetaNoInt!=((nTaxa-nUnSelec)*nPredics))
  {
    stop("Beta dimension from grouped analyis does not match the expected number")
  }

  finalBeta=rep(0,nAlphaSelec)

  for (i in 1:nUnSelec)
  {
    if(i==1){
      unSelec.i=1
      unSelec.i1=unSelectList[i]
    }else{
      unSelec.i=unSelectList[i-1]
      unSelec.i1=unSelectList[i]
    }
    if (nUnSelec==1) {
      if(unSelec.i1==1) {
        finalBeta[(nPredics+1):nAlphaSelec]=newBetaNoInt
      }
      if(unSelec.i1==nTaxa) {
        finalBeta[1:(nAlphaSelec-nPredics)]=newBetaNoInt
      }
      if((unSelec.i1>1) & (unSelec.i1<nTaxa)) {
        finalBeta[1:((unSelec.i1-1)*nPredics)]=newBetaNoInt[1:((unSelec.i1-1)*nPredics)]
        finalBeta[(unSelec.i1*nPredics+1):nAlphaSelec]=newBetaNoInt[((unSelec.i1-1)*nPredics+1):nNewBetaNoInt]
      }
    }
    else {
      distance.i=unSelec.i1-unSelec.i
      if (distance.i==0) next
      if (distance.i>1){
        if(i==1){
          finalBeta[1:((unSelec.i1-1)*nPredics)]=newBetaNoInt[1:((unSelec.i1-1)*nPredics)]
        } else {
          finalBeta[(unSelec.i*nPredics+1):((unSelec.i1-1)*nPredics)]=newBetaNoInt[((unSelec.i-i+1)*nPredics+1):((unSelec.i1-i)*nPredics)]
          if ((i==nUnSelec) & (unSelec.i1<nTaxa)){
            finalBeta[(unSelec.i1*nPredics+1):nAlphaSelec]=newBetaNoInt[((unSelec.i1-i)*nPredics+1):nNewBetaNoInt]
          }
        }
      }

      if((distance.i==1) & (i==1)){
        finalBeta[1:((unSelec.i1-1)*nPredics)]=newBetaNoInt[1:((unSelec.i1-1)*nPredics)]
      }
      if ((distance.i==1) & (i==nUnSelec) & (unSelec.i1<nTaxa)){
        finalBeta[(unSelec.i1*nPredics+1):nAlphaSelec]=newBetaNoInt[((unSelec.i1-i)*nPredics+1):nNewBetaNoInt]
      }
      if ((distance.i==1) & (i>1) & (i<nUnSelec)) next

      if ((distance.i==1) & (i=nUnSelec) & (unSelec.i1=nTaxa)) next
    }
  }

  results$finalBeta=finalBeta
  return(results)
}
