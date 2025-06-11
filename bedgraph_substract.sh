#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: bedgraph_subtract.sh
# Purpose: Subtract methylation values between two bedGraph files to produce
#          a differential methylation track.
#
# Usage: ./bedgraph_subtract.sh <sample1.bedGraph> <sample2.bedGraph>
#
# Inputs:
#   - input_file1.bedGraph: First bedGraph file (e.g. KO sample)
#   - input_file2.bedGraph: Second bedGraph file (e.g. WT sample)
#     Each bedGraph must contain 4 columns: chr, start, end, methylation_value
#
# Outputs:
#   - A file named <sample1>_vs_<sample2>_diff.bedGraph containing:
#     chr, start, end, methylation_difference (sample1 - sample2)
#
# Notes:
#   - Filters to canonical autosomal chromosomes only (chr1–chr22)
#   - Input files must not contain duplicated intervals (start/end) per chrom
#   - Sorted multipe times for bedtools process and to avoid duplication errors
# -----------------------------------------------------------------------------

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

# Filter to autosomal canonical chromosomes (chr1–chr22) and sort by chromosome and start position (imp. for bedtools)
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "$INPUT_FILE1" | sort -k1,1V -k2,2n > "$FILTERED1"
grep -E '^chr([1-9]|1[0-9]|2[0-2])[[:space:]]' "$INPUT_FILE2" | sort -k1,1V -k2,2n > "$FILTERED2"

# Subtract methylation values using bedtools
# - 'unionbedg' aligns by coordinates and fills in 0s if cpg site does not exist in the other file
# - 'awk' subtracts sample2 from sample1
# - sort to clean output for downstream use (e.g., bigWig conversion for igv visualisation)
bedtools unionbedg -i "$FILTERED1" "$FILTERED2" -filler 0 | \
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 - $5}' | sort -k1,1V -k2,2n > "$OUTFILE"

# Clean up temp files
rm "$FILTERED1" "$FILTERED2"

echo "Done. Output written to: $OUTFILE"