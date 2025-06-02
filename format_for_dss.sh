#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: format_for_dss.sh
# Purpose: Convert Oxford Nanopore bedMethyl files to DSS input format
# Usage: ./format_for_dss.sh input_file.bed
# Output: A new file named <input_file>.dssin.gz with 4 columns DSS expects:
#         - chr (column 1)
#         - pos (column 2, 0-based start position)
#         - N (column 10 = Nvalid_cov, total reads at site)
#         - X (column 12 = Nmod, number of methylated reads)
#
# Notes:
#   - DSS expects chr, pos, N, X
#   - Script sorts by chromosome and position
#   - Script compresses output with gzip
# -----------------------------------------------------------------------------

# Check if user provided a file
if [ -z "$1" ]; then
  echo "Usage: $0 input_file.bed"
  echo "Converts a bedMethyl file into DSS input format (chr, pos, N, X)"
  exit 1
fi

INPUT_FILE=$1
# Extract base name (remove path)
BASENAME=$(basename "$INPUT_FILE")
BASENAME_NOEXT="${BASENAME%.*}" # reomve extension
OUTFILE="${BASENAME_NOEXT}.dssin.gz"

echo "Processing file: $INPUT_FILE"
echo "Saving output to: $OUTFILE"

# Extract needed columns: chr, pos, Nvalid_cov (col 10), Nmod (col 12). Sort and compress.
awk 'BEGIN{OFS="\t"} {print $1, $2, $10, $12}' "$INPUT_FILE" | sort -k1,1V -k2,2n | gzip > "$OUTFILE"

# Confirm success
echo "Conversion complete. Output file: $(realpath "$OUTFILE")"
ls -lh "$OUTFILE"
