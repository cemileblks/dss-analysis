#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: tidy_dssin.sh
# Purpose: Filter DSS input files to retain only autosomal canonical chromosomes (chr1–chr22)
# Usage: ./tidy_dssin.sh <input_file>.dssin.gz
# Input: A .dssin.gz file (4-column gzipped DSS input: chr, pos, N, X)
# Output: A gzipped file <input_file>_tidy.dssin.gz containing only chr1 to chr22
#
# Notes:
#   - This removes random contigs (e.g., chr1_KI...) and non autosomal chromosomes (chrX, chrY, chrM)
#   - Keeps lines where the chromosome name starts with chr1 to chr22 followed by a whitespace (tab or space)
# -----------------------------------------------------------------------------

# Check if user provided a file
if [ -z "$1" ]; then
  echo "Usage: $0 input_file.dssin.gz"
  echo "Tidies DSS input file to retain only canonical autosomal chromosomes (chr1-chr22)"
  exit 1
fi

INPUT_FILE=$1

# Validate file extension
if [[ "$INPUT_FILE" != *.dssin.gz ]]; then
  echo "Error: Input file must have a .dssin.gz extension"
  exit 1
fi

OUTFILE="${INPUT_FILE%.dssin.gz}_tidy.dssin.gz"

echo "Processing: $INPUT_FILE"
echo "Output will be: $OUTFILE"

# Process:
# - Decompress input (zcat)
# - Filter for chr1–chr22 only (canonical autosomes)
# - Recompress output

zcat "$INPUT_FILE" | grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' | gzip > "$OUTFILE"

# Confirm success
echo "Done. Tidy DSS input file saved to:"
realpath "$OUTFILE"
ls -lh "$OUTFILE"