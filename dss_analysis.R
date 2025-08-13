## DSS Analysis - Detect DMRs and DMLs between two sample groups
# R script version of the Rmd file

# Input arguments
# Must be tab delimited 4 col input DSS expects
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript dss_analysis.R <input1:Reference> <input2:TEST>")
}

input1 <- args[1] # Reference group
input2 <- args[2] # Test group


# Makes sure the necessary packages are installed before running this document
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("DSS", quietly = TRUE)) {
    BiocManager::install("DSS")
}

# Load required libraries
library(DSS)
require(bsseq)

cat("Comparing methylation between:\n")
cat("  - Reference group:", input1, "\n")
cat("  - Test group:", input2, "\n")

# Read Files
# Load methylation count data: chr, pos, total reads (N), methylated reads (X)
ref_file <- input1
test_file <- input2

ref <- read.table(ref_file, header = FALSE)
test <- read.table(test_file, header = FALSE)

colnames(ref) <- c("chr", "pos", "N", "X")
colnames(test) <- c("chr", "pos", "N", "X")

# Create BSseq object
# This object holds methylation information in the structure DSS requires
BSobj <- makeBSseqData(list(ref, test), c("REF", "TEST"))

# View summary of the BSseq object
cat("BSobj Summary:\n")
BSobj
cat("\n")

# Run DML test with smoothing (recommended for WGBS data and when there are no replicates)
dmlTest.sm <- DMLtest(BSobj, group1 = "REF", group2 = "TEST", smoothing = TRUE, ncores = 4)

cat(
    "Total mehtylated loci that passed coverage threshold and are tested:",
    nrow(dmlTest.sm), "\n\n"
)

dmlTest.sm$end <- dmlTest.sm$pos + 2

canon_chrs <- paste0("chr", 1:22)

dmlTest.sm$chr <- factor(dmlTest.sm$chr, levels = canon_chrs)
# Order by chromosomes (with factor level order) and then by start position
dmlTest.sm <- dmlTest.sm[order(dmlTest.sm$chr, dmlTest.sm$pos), ]

cpg_loci_id_prefix <- "CpG"
dmlTest.sm$CPG_ID <- seq_len(nrow(dmlTest.sm))
dmlTest.sm$name <- paste0(cpg_loci_id_prefix, "_", dmlTest.sm$CPG_ID, "_", toupper(dmlTest.sm$chr))

options(scipen = 0)

# Null hypothesis: the difference in methylation is less than or equal to delta (|mu1 - mu2| ≤ 0.1)
# Code below filters for biologically meaningful changes (at least 10%)
dmls <- callDML(dmlTest.sm, delta = 0.1, p.threshold = 0.01)
cat("Output lines from DML detection\n")
head(dmls)

# For DMLs add +2 to cover the whole CpG site as the locus
dmls$end <- dmls$pos + 2

canon_chrs <- paste0("chr", 1:22)

dmls$chr <- factor(dmls$chr, levels = canon_chrs)
# Order by chromosomes (with factor level order) and then by start position
dmls <- dmls[order(dmls$chr, dmls$pos), ]

options(scipen = 0)
# Detect DMRs based on smoothed DML test results
# Call DMRs using p-value threshold only (no effect size filter)
dmrs <- callDMR(dmlTest.sm, delta = 0.1, p.threshold = 0.01)
cat("Output lines from DMR detection:\n")
head(dmrs)

showOneDMR(dmrs[1, ], BSobj)

# Extend end coordinate to cover the full CpG site at the end (2bp) for IGV visualisation
dmrs$end <- dmrs$end + 2

# Sort DMRs by chromosome and start position
# Set chromosome column as a factor with canonical order
canon_chrs <- paste0("chr", 1:22)
dmrs$chr <- factor(dmrs$chr, levels = canon_chrs)
# Order by chromosomes (with factor level order) and then by start position
dmrs <- dmrs[order(dmrs$chr, dmrs$start), ]

ref_name <- tools::file_path_sans_ext(basename(input1), compression = TRUE)
test_name <- tools::file_path_sans_ext(basename(input2), compression = TRUE)

# Define group label for use in DMR names
test_label <- paste0(test_name)

# Create DMR names
dmr_id_prefix <- paste0(test_label, "_DSS_DMR")
dmrs$DMR_ID <- seq_len(nrow(dmrs))
dmrs$name <- paste0(dmr_id_prefix, "_", dmrs$DMR_ID, "_", toupper(dmrs$chr))


dml_id_prefix <- paste0(test_label, "_DSS_DML")
dmls$DML_ID <- seq_len(nrow(dmls))
dmls$name <- paste0(dml_id_prefix, "_", dmls$DML_ID, "_", toupper(dmls$chr))

dmls_hyper <- subset(dmls, diff < 0)
dmls_hypo <- subset(dmls, diff > 0)

dmls_hyper_id_prefix <- paste0(test_label, "_DSS_DML_HYPER")
dmls_hypo_id_prefix <- paste0(test_label, "_DSS_DML_HYPO")

dmls_hyper$name <- sub("DSS_", "_DML_HYPER_", dmls_hyper$name)
dmls_hypo$name <- sub("DSS_", "_DML_HYPO_", dmls_hypo$name)

# based on the direction of methylation change (diff.Methy):
#   - Negative = hypermethylation in tested group
#   - Positive = hypomethylation in tested group

dmrs_hyper <- subset(dmrs, diff.Methy < 0)
dmrs_hypo <- subset(dmrs, diff.Methy > 0)

dmrs_hyper_id_prefix <- paste0(test_label, "_DSS_DMR_HYPER")
dmrs_hypo_id_prefix <- paste0(test_label, "_DSS_DMR_HYPO")

dmrs_hyper$name <- sub("DMR_", "DMR_HYPER_", dmrs_hyper$name)
dmrs_hypo$name <- sub("DMR_", "DMR_HYPO_", dmrs_hypo$name)

cat("Output Summary:\n")
cat("  Total methylated Loci:", nrow(dmlTest.sm), "\n")
cat("  Total DMLs:", nrow(dmls), "\n")
cat("    Hypermethylated:", nrow(dmls_hyper), "\n")
cat("    Hypomethylated:", nrow(dmls_hypo), "\n")
cat("  Total DMRs:", nrow(dmrs), "\n")
cat("    Hypermethylated:", nrow(dmrs_hyper), "\n")
cat("    Hypomethylated:", nrow(dmrs_hypo), "\n")

options(scipen = 999) # to prevent scientific notation on tables

today <- format(Sys.Date(), "%d_%m_%Y")

# Create file prefix
prefix <- paste0("HG38_", today, "_", toupper(test_name), "_VS_", toupper(ref_name))

all_loci_file <- paste0(prefix, "_CpG_Loci.bed")
dml_file <- paste0(prefix, "_DMLs.bed")
# Additional hypo/hyper output files
dml_file_hyper <- paste0(prefix, "_DMLs_HYPER.bed")
dml_file_hypo <- paste0(prefix, "_DMLs_HYPO.bed")


dmr_file <- paste0(prefix, "_DMRs.bed")
# Additional hypo/hyper output files
dmr_file_hyper <- paste0(prefix, "_DMRs_HYPER.bed")
dmr_file_hypo <- paste0(prefix, "_DMRs_HYPO.bed")

# Write results to files
write.table(dmlTest.sm[, c("chr", "pos", "end", "name")],
    file = all_loci_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmls[, c("chr", "pos", "end", "name")],
    file = dml_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmls_hyper[, c("chr", "pos", "end", "name")],
    file = dml_file_hyper,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmls_hypo[, c("chr", "pos", "end", "name")],
    file = dml_file_hypo,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmrs[, c("chr", "start", "end", "name")],
    file = dmr_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmrs_hyper[, c("chr", "start", "end", "name")],
    file = dmr_file_hyper,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(dmrs_hypo[, c("chr", "start", "end", "name")],
    file = dmr_file_hypo,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Outputs
cat("Saving all methylated loci to:", all_loci_file, "\n")
cat("Saving DMLs to:", dml_file, "\n")
cat("Saving hyper DMLs to:", dml_file_hyper, "\n")
cat("Saving hypo DMLs to:", dml_file_hypo, "\n")
cat("Saving DMRs to:", dmr_file, "\n")
cat("Saving hyper DMRs to:", dmr_file_hyper, "\n")
cat("Saving hypo DMRs to:", dmr_file_hypo, "\n")
cat("Done.")
