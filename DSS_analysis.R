
# Install package
BiocManager::install("DSS", ask = FALSE, update = FALSE)

# Load modules
library(DSS)
require(bsseq)

# Step 1. Read in text files and create an object of BSseq class
path = file.path(system.file(package="DSS"), "extdata")
dat1.1 = read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 = read.table(file.path(path, "cond1_2.txt"), header=TRUE)
dat2.1 = read.table(file.path(path, "cond2_1.txt"), header=TRUE)
dat2.2 = read.table(file.path(path, "cond2_2.txt"), header=TRUE)
BSobj = makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2), c("C1","C2", "N1", "N2") )[1:1000,]
BSobj

# Step2. Perform statistical test for DML calling DMLtest function
#        1. Estimates mean methylation levels for all CpG sites
#        2. Estimates dispersions at each CpG sites
#        3. Conduct Wald test