#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input_file1.bedGraph input_file2.bedGraph"
  echo "Subtracts methylation values (col 4) from two bedGraph files"
  exit 1
fi

INPUT_FILE1="$1"
INPUT_FILE2="$2"

# Output filename
BASENAME1=$(basename "$INPUT_FILE1" .bedGraph)
BASENAME2=$(basename "$INPUT_FILE2" .bedGraph)
OUTFILE="${BASENAME1}_vs_${BASENAME2}_diff.bedGraph"

# Temporary filtered files
FILTERED1="${INPUT_FILE1%.bedGraph}_filtered.bedGraph"
FILTERED2="${INPUT_FILE2%.bedGraph}_filtered.bedGraph"

# Step 1: Filter to autosomal canonical chromosomes only
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "$INPUT_FILE1" > "$FILTERED1"
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "$INPUT_FILE2" > "$FILTERED2"

# Step 2: Subtract methylation values
bedtools unionbedg -i "$FILTERED1" "$FILTERED2" -filler 0 | \
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 - $5}' | sort -k1,1V -k2,2n > "$OUTFILE"

# Clean up
rm "$FILTERED1" "$FILTERED2"

echo "Done. Output written to: $OUTFILE"