#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: bedmethyl_to_bed.sh
# Purpose: Convert bedMethyl files to a standard BED format with full CpG site
# Usage:   ./bedmethyl_to_bed.sh input_file.bed
# Output:  A new file named <input_file>_CpG.bed with 4 columns:
#           - chr (column 1)
#           - start (column 2)
#           - end+1 (column 3 incremented by 1 to capture full CpG site)
#           - name or other identifier (column 4)
#
# Notes:
#   - This script increments the 3rd column to include the full CpG dinucleotide.
#   - Assumes tab-delimited bedMethyl input with at least 4 columns.
# -----------------------------------------------------------------------------

if [ -z "$1" ]; then
  echo "Usage: $0 input_file.bed"
  echo "Converts a bedMethyl file to bed file with each CpG site"
  exit 1
fi

INPUT_FILE=$1
# Extract base name (remove path)
BASENAME=$(basename "$INPUT_FILE")
BASENAME_NOEXT="${BASENAME%.*}" # reomve extension
OUTFILE="${BASENAME_NOEXT}_CpG.bed"

echo "Processing file: $INPUT_FILE"
echo "Saving output to: $OUTFILE"

# Check and process input
awk 'BEGIN{OFS="\t"} NF>=3 && $3 ~ /^[0-9]+$/ {print $1, $2, $3+1, $4}' "$INPUT_FILE" > "$OUTFILE"

if [ $? -eq 0 ]; then
  echo "Conversion complete. Output file: $(realpath "$OUTFILE")"
  ls -lh "$OUTFILE"
else
  echo "Error during processing. Output file not created."
fi