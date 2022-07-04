groupBetaToFullBeta = function(nTaxa,
                               nPredics,
                               unSelectList,
                               newBetaNoInt) {
  results = list()
  unSelectList = unique(sort(unSelectList))
  nUnSelec = length(unSelectList)
  nAlphaSelec = nTaxa * nPredics
  nNewBetaNoInt = length(newBetaNoInt)
  
  if (nNewBetaNoInt != ((nTaxa - nUnSelec) * nPredics)) {
    stop("Error: Beta dimension from grouped analyis does not match the expected number")
  }
  
  if (nTaxa < max(unSelectList) | 1 > min(unSelectList)) {
    stop("Error: unSelectList out of range")
  }
  finalBeta = newBetaNoInt
  for (i in unSelectList) {
    finalBetaTemp = rep(NA, (length(finalBeta) + nPredics))
    lengthTemp = length(finalBetaTemp)
    
    if (i == 1) {
      finalBetaTemp[seq_len(nPredics)] = 0
      finalBetaTemp[(nPredics + 1):lengthTemp] = finalBeta
      finalBeta = finalBetaTemp
    }
    
    if (i > 1) {
      if ((i * nPredics) <= (length(finalBeta))) {
        finalBetaTemp[seq_len(((i - 1) * nPredics))] = finalBeta[seq_len(((i - 1) *
                                                                            nPredics))]
        finalBetaTemp[((i - 1) * nPredics + 1):(i * nPredics)] = 0
        finalBetaTemp[(i * nPredics + 1):lengthTemp] = finalBeta[((i - 1) *
                                                                    nPredics + 1):(length(finalBeta))]
      } else{
        finalBetaTemp[seq_len(((i - 1) * nPredics))] = finalBeta
        finalBetaTemp[((i - 1) * nPredics + 1):lengthTemp] = 0
      }
      finalBeta = finalBetaTemp
    }
  }
  rm(finalBetaTemp)
  results$finalBeta = finalBeta
  rm(finalBeta)
  return(results)
}
