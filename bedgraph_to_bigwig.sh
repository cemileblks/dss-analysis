#!/usr/bin/bash

# -----------------------------------------------------------------------------
# Script: bedgraph_to_bigwig.sh
# Purpose: Convert a bedGraph file to a BigWig file after adjusting intervals
#          to span the whole CpG site
#
# Usage: ./bedgraph_to_bigwig.sh <input.bedGraph> <chrom.sizes>
#
# Inputs:
#   - input.bedGraph: A 4 column bedGraph file (chr, start, end, value)
#   - chrom.sizes: A file containing chromosome sizes (2 columns, e.g. from UCSC hg38)
#
# Output: A BigWig file: <input>.bw
#
# Notes:
#   - This script assumes each methylation site covers a single base. It adds +1 to the 'end' field.
#   - bedGraph must be sorted by chr/start before conversion.
#   - If you do not have the 'bedGraphToBigWig' binary, dowload instructions can be found at
#     https://genome.ucsc.edu/goldenpath/help/bigWig.html
#   - Adapted from submit_percent_bedgraph_to_bigWig_CpG.sh and https://genome.ucsc.edu/goldenpath/help/bigWig.html#Ex3
# -----------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input.bedGraph chrom.sizes"
  echo "Converts bedGraph to BigWig, adjusting intervals to full CpG coverage"
  exit 1
fi

# Input files
BEDGRAPH="$1"
CHROMSIZES="$2"

# Output filenames
BASENAME=$(basename "$BEDGRAPH" .bedGraph)
BIGWIG="${BASENAME}.bw"
SORTED_BEDGRAPH="${BASENAME}_modified_n_sorted.bedGraph"

# Adjust end coordinate to include the whole CpG site and sort
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3+1, $4}' "$BEDGRAPH" | sort -k1,1V -k2,2n > "$SORTED_BEDGRAPH"

# Convert to BigWig
bedGraphToBigWig "$SORTED_BEDGRAPH" "$CHROMSIZES" "$BIGWIG"

echo "Done. BigWig file created: $BIGWIG"