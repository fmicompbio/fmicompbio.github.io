# install necessary bioconductor packages
# (only first time, see www.bioconductor.org)
source("http://www.bioconductor.org/biocLite.R")
biocLite("oligo")
biocLite("pd.huex.1.0.st.v2")

# normalize CEL files using oligo package (version 1.10.0 or newer)
library(oligo)
celFiles <- list.celfiles('/path/to/my/celfiles', full.name=TRUE)
raw <- read.celfiles(celFiles, pkgname="pd.huex.1.0.st.v2")
eset <- rma(raw, target="probeset")

# alternatively, load normalized log2-expression data from text file
# to obtain an ExpressionSet object with probeset expression levels:
#library(affy)
#eset <- readExpressionSet('/path/to/my/expression/file.txt')

# load cosie functions
source('cosie.R')

# calculate corrected splicing indices
cosieOut <- cosie(eset, '/path/to/cosie/parameter/file.txt')

# write cosie results to tab-delimited text file with columns:
#   psId      : probeset identifier
#   tcId      : transcript cluster identifier
#   SI_*      : corrected pre-splicing index
#               (one column per sample, sample name in '*',
#                subtract two of these columns in order to
#                obtain the differential splicing index)
#   tcLevel_* : trancript cluster expression level
#               (one column per sample, sample name in '*')
write.cosie(cosieOut, 'cosieOut.txt')

