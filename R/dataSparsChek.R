
dataSparsCheck=function(
  data,
  Mprefix
){
  results=list()

  # get the original sample size
  nSub=nrow(data)
  MVarNamLength=nchar(Mprefix)

  # get taxa variable names
  micros = sapply(substr(colnames(data),1,MVarNamLength), function(x) {grep(Mprefix, x)})
  microPositions=which(micros == 1)
  rm(micros)

  taxaNames=colnames(data)[microPositions]
  rm(microPositions)

  w=data[,taxaNames,drop=FALSE]
  rm(data,taxaNames)
  overallSparsity=round(100*sum(w==0)/(nrow(w)*ncol(w)),2)
  cat(overallSparsity,"percent of microbiome sequencing reads are zero","\n")

  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads=length(which(Matrix::colSums(w)==0))
  if(numTaxaNoReads>0){
    cat("There are",numTaxaNoReads,"taxa without any sequencing reads and
        excluded from the analysis","\n")
  }
  rm(numTaxaNoReads)

  numSubNoReads=length(which(Matrix::rowSums(w)==0))
  if(numSubNoReads>0){
    cat("There are",numSubNoReads,"subjects without any sequencing reads and
        excluded from the analysis","\n")
  }
  rm(numSubNoReads,w)
}
