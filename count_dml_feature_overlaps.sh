#!/usr/bin/bash

# Usage: ./count_dml_feature_overlaps.sh cpg.bed hyper.bed hypo.bed output_dir/

# Input files
CPG_FILE=$1
HYPER_DML_FILE=$2
HYPO_DML_FILE=$3
OUT_DIR=${4:-"dml_feature_overlap_output"}

# Annotation directory
ANNOT_DIR="./data/annotation_data"

# Create output directory if not exists
mkdir -p "$OUT_DIR"

# Define features to test
FEATURES=("tss" "exons" "introns" "5utrs" "3utrs" "downstream_1kb" "intergenic")

# Write header to output CSV
OUT_CSV="$OUT_DIR/overlap_counts.csv"
# echo "Feature,All_CpG_Count,HyperDMLs,HypoDMLs" >"$OUT_CSV"

# Function to count unique overlaps
count_overlaps() {
    local loci_file=$1
    local feature_file=$2

    bedtools intersect -u -a "$loci_file" -b "$feature_file" | wc -l
}

# Process each feature
for feature in "${FEATURES[@]}"; do
    ANNOT_FILE="$ANNOT_DIR/hg38_${feature}.bed"

    if [[ ! -f "$ANNOT_FILE" ]]; then
        echo "Warning: Missing annotation file: $ANNOT_FILE" >&2
        continue
    fi

    # Compute counts
    CPG_COUNT=$(count_overlaps "$CPG_FILE" "$ANNOT_FILE")
    HYPER_COUNT=$(count_overlaps "$HYPER_DML_FILE" "$ANNOT_FILE")
    HYPO_COUNT=$(count_overlaps "$HYPO_DML_FILE" "$ANNOT_FILE")

    # Append to CSV
    echo "$feature,$CPG_COUNT,$HYPER_COUNT,$HYPO_COUNT" >>"$OUT_CSV"
done

# CPG_COUNT=$(wc -l < "$CPG_FILE")
# HYPER_COUNT=$(wc -l < "$HYPER_DML_FILE")
# HYPO_COUNT=$(wc -l < "$HYPO_DML_FILE")
# echo "total,$CPG_COUNT,$HYPER_COUNT,$HYPO_COUNT" >>"$OUT_CSV"


echo "Done. Output written to: $OUT_CSV"
