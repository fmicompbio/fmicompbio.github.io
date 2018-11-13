# load data (GSE36916)
dat  <- as.matrix(read.delim("Fig2bcd_GSE36916_normalized_forPhaseFitting.txt.gz", row.names=1))
ex1  <- dat[,grep("^exon.Nascent",colnames(dat))]  # Nascent RNA
ex2  <- dat[,grep("^exon.RNASeq",colnames(dat))]   # mRNA (exons)
intr <- dat[,grep("^intron.RNASeq",colnames(dat))] # mRNA (introns)

# define function to fit amplitudes/phases
fitAmplPhase <- function(TLE, lambda=6) {
    # TLE: expression matrix (rows: genes, columns: timepoints)
    # lambda: number of timepoints per cycle
    RM <- cbind(X1=cos(2*pi/lambda*(0:(ncol(TLE)-1))), X2=-sin(2*pi/lambda*(0:(ncol(TLE)-1))))
    CAB <- matrix(0, nrow(TLE), ncol=(ncol(RM)+1), dimnames=list(rownames(TLE), c("Const","A","B")))
    for(i in 1:nrow(TLE))
        CAB[i,] <- lm(TLE[i,] ~ RM)$coefficient
    data.frame(amplitude=sqrt(rowSums(CAB[,2:3]^2)), phase=-atan2(CAB[,3],CAB[,2]), row.names=rownames(CAB))
}

# fit phases and amplitudes
AmplPhase <- list()
AmplPhase[['exNascent']] <- fitAmplPhase(ex1, 6)
AmplPhase[['exTotal']]   <- fitAmplPhase(ex2, 6)
AmplPhase[['intrTotal']] <- fitAmplPhase(intr, 6)
AmplPhaseM <- do.call(cbind, AmplPhase) # combine

