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

# combine exon and intron raw counts for edgeR containing genes with sufficient counts in both exonic and intronic levels
cnt <- cbind(Ex=as.data.frame(cntEx[genes.sel,]), In=cntIn[genes.sel,])

# edgeR workflow
factorRegion    <- factor(rep(c("ex","in"),each=ncol(cntEx)), levels=c("in", "ex")) # define experimental factor exon/intron 
factorCondition <- factor(rep(conditions,2), levels=unique(conditions)) # define experimental factor 'conditions'

y <- DGEList(counts=cnt, genes=data.frame(ENTREZID=rownames(cnt))) # define DGEList object
y <- calcNormFactors(y) # determine normalization factors
design <- model.matrix(~ factorRegion * factorCondition) # design matrix with interaction term
rownames(design) <- colnames(cnt)

y <- estimateDisp(y, design) # estimate dispersion
fit <- glmFit(y, design) # fit generalized linear model
lrt <- glmLRT(fit) # calculate likelihood-ratio between full and reduced models
#final table with significance level for each gene 
tt <- topTags(lrt, n=nrow(y))
head(tt$table)
