# Install DSS package from Bioconductor
# (ask = FALSE and update = FALSE to avoid full system updates)
BiocManager::install("DSS", ask = FALSE, update = FALSE)

# Load required libraries
library(DSS)
require(bsseq)

# Step 1: Read in example methylation count files
# 4 example BS-seq samples: 2 from each of 2 groups (e.g. control vs treatment)
path = file.path(system.file(package="DSS"), "extdata")
dat1.1 = read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 = read.table(file.path(path, "cond1_2.txt"), header=TRUE)
dat2.1 = read.table(file.path(path, "cond2_1.txt"), header=TRUE)
dat2.2 = read.table(file.path(path, "cond2_2.txt"), header=TRUE)

# Create a BSseq object
# This merges all CpG sites found in any sample into one common coordinate space
# Missing sites are handled with 0 coverage
# Here we limit to the first 1000 loci for faster testing
BSobj = makeBSseqData(
  list(dat1.1, dat1.2, dat2.1, dat2.2),
  c("C1","C2", "N1", "N2")
  )[1:1000,]

# Check the BSseq object: how many CpGs, how many samples
BSobj

# Step 2: Differential Methylation Loci (DML) test
# This tests whether CpG sites show statistically significant differences
# in methylation between two sample groups: C1/C2 (e.g., Control) vs N1/N2 (e.g., KO)

# DMLtest() does 3 things at each CpG:
#   1. Estimates mean methylation level in each group (mu1 and mu2)
#   2. Estimates variability between replicates (dispersion: phi1 and phi2)
#   3. Performs a Wald test for group difference (returns test statistic and p-value)

# -- Run test WITHOUT smoothing (based on raw values per CpG)
dmlTest = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))

# -- Run test WITH smoothing (uses methylation info from nearby CpGs)
#    This helps get more stable results when methylation is noisy
dmlTest.sm = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"), 
                     smoothing=TRUE)

# Step 3: Identify DMLs based on statistical test results

# Null hypothesis: the diff. in methylation levels is 0 between groups (mu1 = mu2)
# No minimum effect size required
dmls <- callDML(dmlTest, p.threshold = 0.001)
head(dmls)

# Null hypothesis: the difference in methylation is less than or equal to delta (|mu1 - mu2| ≤ 0.1)
# This filters for biologically meaningful changes (at least 10%)
# 1 - posterior probability is used like a p-value threshold
dmls2 <- callDML(dmlTest, delta = 0.1, p.threshold = 0.001)
head(dmls2)

# Step 4: Identify DMRs by grouping nearby significant CpGs (DMLs)

# Basic DMR calling:
# - Uses p-value threshold to define significant CpGs
# - Groups significant CpGs that are close together into regions
dmrs <- callDMR(dmlTest, p.threshold = 0.01)
head(dmrs)

# Alternative version using delta:
# - Only consider CpGs with at least 10% methylation difference
# - More stringent: makes sure the DMRs have biologically meaningful size
dmrs2 <- callDMR(dmlTest, delta = 0.1, p.threshold = 0.05)
head(dmrs2)

# Plot the first DMR: shows methylation levels and coverage for each CpG
showOneDMR(dmrs[1, ], BSobj)
