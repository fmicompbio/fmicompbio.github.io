cosie <- function(x, paramFile=NULL, exprExt=2) {
  # calculate corrected splicing indices (v 1.1)
  # arguments:
  #   x         : log2 probeset (exon) expression levels, one of:
  #                 - ExpressionSet (with feature and sample names)
  #                 - matrix (with row and column names)
  #   paramFile : filename with corrector parameters
  #   exprExt   : extend expression ranges observed in training set by
  #               this many log2-units on both sides to determine
  #               correctable expression levels
  # return value:
  #   list of three components:
  #      correctedSI : p-times-s matrix with corrected splicing indices
  #                    (p: correctable probesets, s: samples)
  #      transcriptClusterLevels : t-times-s matrix with transcript cluster levels
  #                                used in calculating splicing indices
  #                                (t: transcript clusters, s: samples)
  #      tcIds : p-element vector with transcript cluster ids corresponding to probesets

  # get sample and probeset ids
  if(is.matrix(x)) {
    sampleIds <- if(is.null(colnames(x))) as.character(1:ncol(x)) else colnames(x)
    psIds <- rownames(x)
  } else if(class(x) == 'ExpressionSet') {
    sampleIds <- sampleNames(x)
    psIds <- featureNames(x)
    x <- exprs(x)
  } else {
    stop("x needs to be a matrix or object of class 'ExpressionSet'", call.=FALSE)
  }
  nSamples <- length(sampleIds)

  # check input parameters
  if(is.null(paramFile)) {
    stop("no parameter file 'paramFile' specified. Aborting.", call.=FALSE)
  }
  if(any(!is.finite(x))) {
    stop("x must not have any missing or non-finite values", call.=FALSE)
  }
  if(max(x) > 20) {
    warning("x contains values larger than 20, and therefore is likely not to be in log2", call.=FALSE, immediate.=TRUE)
  }

  # calculate mean probeset expression over samples
  psExprMin <- 3.5 # only used to report expressed but uncorrectable probesets
  meanPsExpr <- rowMeans(x)
  
  # load parameter file
  param <- scan(paramFile, sep='\t', dec='.', quiet=TRUE, multi.line=FALSE, 
                what=list(psId="", tcId="",
                  status=integer(0), geneTrainFrom=numeric(0), geneTrainTo=numeric(0),
                  a=numeric(0), b=numeric(0), c=numeric(0), d=numeric(0)))
  x2param <- match(param$psId, psIds)
  if(any(is.na(x2param))) {
    #missing <- param$psId[is.na(x2param)]
    #stop(length(missing), " probesets missing in expression data x: ", missing[1], ", ...")
    toRm.tcId <- unique(param$tcId[is.na(x2param)])
    toRm <- which(param$tcId %in% toRm.tcId)
    toRm.psId <- param$psId[toRm]
    warning(paste('Removing',length(toRm.tcId), 'transcript clusters (', length(which(is.na(x2param))), 'probesets) with probeset composition deviating from training data'), call.=FALSE, immediate.=TRUE)
    for(nm in names(param))
      param[[nm]] <- param[[nm]][-toRm]
    x2param <- match(param$psId, psIds)
    rm(toRm.tcId, toRm, toRm.psId)
  }
  psExprInTestButNotInTraining <- length(which(param$status==1 & meanPsExpr[match(param$psId,psIds)]>psExprMin))

  # select needed probesets
  toKeep <- which(param$status==0 | param$status==3)
  x2param <- x2param[toKeep]
  probesetsThatHaveToBeCorrected <- length(toKeep) + psExprInTestButNotInTraining
  for(nm in names(param))
    param[[nm]] <- param[[nm]][toKeep]
  param$tcId <- factor(param$tcId)
  
  # order x/psIds/meanPsExpr according to param
  x <- x[x2param,]
  psIds <- psIds[x2param]
  meanPsExpr <- meanPsExpr[x2param]

  # calculate transcript cluster (gene) expression levels
  ps.by.tc <- split(1:nrow(x), param$tcId)
  calcGnLevel <- function(j) { colMeans(x[j,]) }
  geneLevels <- matrix(sapply(ps.by.tc, calcGnLevel, USE.NAMES=FALSE),
                       ncol=nSamples, byrow=TRUE, dimnames=list(names(ps.by.tc), sampleIds))

  # check for correctable probesets
  toCorrect <- which(apply(geneLevels,1,min)[param$tcId]>param$geneTrainFrom - exprExt &
                     apply(geneLevels,1,max)[param$tcId]<param$geneTrainTo + exprExt &
                     param$status!=3)
  tcIds <- unique(param$tcId[toCorrect])
  probesetsThatCouldBeCorrected <- length(toCorrect)

  # correct probesets
  x <- x - geneLevels[param$tcId,] # subtract gene level from probeset level in all experiments
  x <- x[toCorrect,]
  meanPsExpr <- meanPsExpr[toCorrect]
  for(nm in names(param))
    param[[nm]] <- param[[nm]][toCorrect]
  
  csi <- matrix(NA_real_, ncol=nSamples, nrow=probesetsThatCouldBeCorrected,
                dimnames=list(psIds[toCorrect], sampleIds))
  for(i in 1:nSamples) {
    csi[,i] <- x[,i] + (param$a/param$d*log(exp(-param$d*geneLevels[param$tcId,i]) + exp(-param$d*param$c)) - param$b)
  }
  csi <- csi -rowMeans(csi) +meanPsExpr # set the mean to the probeset mean
  
  # return results
  correctedPercent <- round(100*probesetsThatCouldBeCorrected/probesetsThatHaveToBeCorrected, 2)
  cat(paste("Number of probesets that need correction:", probesetsThatHaveToBeCorrected, "\n"))
  cat(paste("Number of probesets that were corrected :", probesetsThatCouldBeCorrected, "(", correctedPercent, "%)\n\n"))
  list(correctedSI=csi, transcriptClusterLevels=geneLevels[tcIds,], tcIds=param$tcId)                
}


write.cosie <- function(x, fname) {
  dat <- data.frame(psId=rownames(x[[1]]), tcId=as.character(x[[3]]))
  dat <- cbind(dat, x[[1]], x[[2]][as.character(x[[3]]),])
  nSamples <- ncol(x[[1]])
  colnames(dat)[3:(nSamples+2)] <- sprintf("SI_%s", colnames(x[[2]]))
  colnames(dat)[(nSamples+3):(2*nSamples+2)] <- sprintf("tcLevel_%s", colnames(x[[2]]))
  write.table(dat, file=fname, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

