# load edgeR library
library(edgeR)

# input files and parameters
exsFile <- "Fig3abc_GSE33252_rawcounts_exonic.txt.gz"
insFile <- "Fig3abc_GSE33252_rawcounts_intronic.txt.gz"
gnsFile <- "nonOverlappingGenes_mm9_stranded.txt.gz"
conditions <- c("ES","ES","TN","TN") # correspond to the columns in exsFile and insFile; the first condition will be the reference

# read in count tables and remove the first column (feature width)
cntEx <- read.delim(exsFile)[,-1]
cntIn <- read.delim(insFile)[,-1]
genes_no_overlap <- scan(gnsFile, what="")

# find non-overlapping genes with sufficient exonic and intronic counts (genes.sel)
cntEx.norm <- t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx)) # normalize samples to avearge sequencing depth for exons
cntIn.norm <- t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn)) # normalize samples to avearge sequencing depth for introns
genes.sel <- rownames(cntEx) %in% genes_no_overlap & rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5

# mean-variance relationship for exons and introns (just take the two ES replicates)
yEx <- DGEList(counts=cntEx[genes.sel,1:2])
yIn <- DGEList(counts=cntIn[genes.sel,1:2])
#design.simple <- model.matrix(~ conditions)
#rownames(design.simple) <- colnames(cntEx)

#yEx <- estimateDisp(yEx, design.simple)
#yIn <- estimateDisp(yIn, design.simple)
yEx <- estimateDisp(yEx)
yIn <- estimateDisp(yIn)

yEx <- calcNormFactors(yEx)
yIn <- calcNormFactors(yIn)

png("EISA_edgeR_R-code_plotMeanVars.png", height=500, width=900, pointsize=16)
par(mfrow=c(1,2))
xlims <- c(16,60e3)
ylims <- c(40,5.8e8)
plotMeanVar(yEx,main="Exons",
            show.raw.vars=FALSE, # gray circles
            show.tagwise.vars=TRUE, # blue circles
            show.binned.common.disp.vars=FALSE, # red crosses
            show.ave.raw.vars=TRUE, # brown crosses
            NBline=TRUE, # dark blue line
            xlim=xlims, ylim=ylims)
legend(x="topleft", bty="n", lwd=2, lty=c(1,1,NA,NA), pch=c(NA,NA,1,4), col=c("black","dodgerblue3","lightskyblue","darkred"), legend=c("Poisson","Neg.Bin.","tagwise vars","average raw vars"))
plotMeanVar(yIn, main="Introns",
            show.raw.vars=FALSE, # gray circles
            show.tagwise.vars=TRUE, # blue circles
            show.binned.common.disp.vars=FALSE, # red crosses
            show.ave.raw.vars=TRUE, # brown crosses
            NBline=TRUE, # dark blue line
            xlim=xlims, ylim=ylims)
legend(x="topleft", bty="n", lwd=2, lty=c(1,1,NA,NA), pch=c(NA,NA,1,4), col=c("black","dodgerblue3","lightskyblue","darkred"), legend=c("Poisson","Neg.Bin.","tagwise vars","average raw vars"))
dev.off()
