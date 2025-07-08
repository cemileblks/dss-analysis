#!/usr/bin/bash
# Exit if any part of code fails
set -euo pipefail

# Input arguments
CPG_FILE="$1"
HYPER_FILE="$2"
HYPO_FILE="$3"
OUT_DIR="${4:-dml_feature_overlap_output}" # Oprional output dir
ANNOT_DIR="./data/annotation_data"

# Features (in order of precedence for exclusive assignment) 
FEATURES=("tss" "downstream_1kb" "5utrs" "3utrs" "exons" "introns" "intergenic")

# Create output directory
mkdir -p "$OUT_DIR"

# Output CSV files
CSV_RAW="$OUT_DIR/overlap_counts.csv"
CSV_EXCLUSIVE="$OUT_DIR/overlap_counts_exclusive.csv"

# Header of CSV file
echo "Feature,All_CpG_Count,HyperDMLs,HypoDMLs" > "$CSV_RAW"
echo "Feature,All_CpG_Count,HyperDMLs,HypoDMLs" > "$CSV_EXCLUSIVE"

# Function to count overlaps (non-exclusive)
count_overlaps() {
    local query=$1
    local target=$2
    # How many unique regions in query overlap the features in target
    bedtools intersect -u -a "$query" -b "$target" | wc -l
}

# Non-exclusive overlap counts loop
for feature in "${FEATURES[@]}"; do
    FEATURE_FILE="$ANNOT_DIR/hg38_${feature}.bed"

    CPG_COUNT=$(count_overlaps "$CPG_FILE" "$FEATURE_FILE")
    HYPER_COUNT=$(count_overlaps "$HYPER_FILE" "$FEATURE_FILE")
    HYPO_COUNT=$(count_overlaps "$HYPO_FILE" "$FEATURE_FILE")

    echo "$feature,$CPG_COUNT,$HYPER_COUNT,$HYPO_COUNT" >> "$CSV_RAW"
done

# Function to assign loci exclusively
assign_exclusive() {
    local input=$1
    local prefix=$2
    local remaining="$OUT_DIR/${prefix}_remaining.bed"

    cp "$input" "$remaining"

    for feature in "${FEATURES[@]}"; do
        FEATURE_FILE="$ANNOT_DIR/hg38_${feature}.bed"
        OUTFILE="$OUT_DIR/${prefix}_${feature}.bed"

        # find all the sites that have unique overlaps
        bedtools intersect -u -a "$remaining" -b "$FEATURE_FILE" > "$OUTFILE"
        # remove all the sites (-v) from the remaining cpg sites
        bedtools intersect -v -a "$remaining" -b "$FEATURE_FILE" > "${remaining}.tmp"
        mv "${remaining}.tmp" "$remaining"
    done

    rm -f "$remaining"
}

# Assign CpGs and DMLs to features (exclusively)
assign_exclusive "$CPG_FILE" "exclusive_cpg"
assign_exclusive "$HYPER_FILE" "exclusive_hyper"
assign_exclusive "$HYPO_FILE" "exclusive_hypo"

# count exclusive overlaps in files and wrie to csv file
for feature in "${FEATURES[@]}"; do
    CPG_EX="$OUT_DIR/exclusive_cpg_${feature}.bed"
    HYPER_EX="$OUT_DIR/exclusive_hyper_${feature}.bed"
    HYPO_EX="$OUT_DIR/exclusive_hypo_${feature}.bed"

    CPG_COUNT=$(wc -l < "$CPG_EX" 2>/dev/null || echo 0) # the case handles if there are no overlaps and returns 0
    HYPER_COUNT=$(wc -l < "$HYPER_EX" 2>/dev/null || echo 0)
    HYPO_COUNT=$(wc -l < "$HYPO_EX" 2>/dev/null || echo 0)

    echo "$feature,$CPG_COUNT,$HYPER_COUNT,$HYPO_COUNT" >> "$CSV_EXCLUSIVE"
done

echo "Overlap count complete."
echo "Non-exclusive:     $CSV_RAW"
echo "Exclusive (unique): $CSV_EXCLUSIVE"