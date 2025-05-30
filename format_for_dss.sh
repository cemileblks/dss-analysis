#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: format_for_dss.sh
# Purpose: Convert Oxford Nanopore bedMethyl files to DSS input format
# Usage: ./format_for_dss.sh input_file.bed
# Output: A new file named dss_data_<input_file>.bed with 4 columns:
#         - chr (column 1)
#         - pos (column 2, 0-based start position)
#         - N (column 10 = Nvalid_cov, total reads at site)
#         - X (column 12 = Nmod, number of methylated reads)
#
# Note:
#   - DSS expects chr, pos, N, X
# -----------------------------------------------------------------------------

# Check if user provided a file
if [ -z "$1" ]; then
  echo "Usage: $0 input_file.bed"
  echo "Converts a bedMethyl file into DSS input format (chr, pos, N, X)"
  exit 1
fi

INPUT_FILE=$1
# Extract base name (remove path) amd change file extenstion to dssin to signify file input for dss
BASENAME=$(basename "$INPUT_FILE")
OUTFILE="dss_data_${BASENAME%.*}.dssin"

echo "Processing file: $INPUT_FILE"
echo "Saving output to: $OUTFILE"

# Extract needed columns: chr, pos, Nvalid_cov (col 10), Nmod (col 12)
awk 'BEGIN{OFS="\t"} {print $1, $2, $10, $12}' "$INPUT_FILE" > "$OUTFILE"

# Confirm success
echo "Conversion complete. Output file: $(realpath "$OUTFILE")"
ls -lh "$OUTFILE"

# UPDATES TO MAKE:
# CHECK IF IT IS OXFORD NANAOPORE BEDMETHYL FILE (CHECK NUMBER OF COLOUMNS ARE 18 AND AND IF POS, N AND X ARE NUMBERS) 
# SORT BY CHROMOSOME AND POSITION BEFORE SAVING TO OUTPUT FILE 
# CHANGE THE OUTPUT FILE TO .DSS INSTEAD OF BED FILES BECAUSE THEY ARE NOT BED FILES
# GZIP THE FILE OUTPUT BECAUSE R CAN READ GZIPPED FILES AND TO SAVE SPACE