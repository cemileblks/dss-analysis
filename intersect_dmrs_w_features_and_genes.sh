#!/usr/bin/bash

set -euo pipefail

# Defaults
HYPER_FILE=""
HYPO_FILE=""
DMR_FILE=""
OUT_DIR="dmr_feature_gene_mapping_output"
ANNOT_DIR="./data/annotation_data"
# List of genomic features to test against DMRs
FEATURES=("tss" "downstream_1kb" "5utrs" "3utrs" "exons" "introns" "intergenic")
# Path to 4 column BED file of protein-coding genes (Ensembl ID in 4th column)
GENE_BED="$ANNOT_DIR/hg38_protein_coding_genes.bed"

# Help message function
usage() {
    echo "Usage: $0 [--hyper FILE] [--hypo FILE] [--dmr FILE] [--out DIR]"
    echo "  --hyper   Hyper DMR BED file (optional)"
    echo "  --hypo    Hypo DMR BED file (optional)"
    echo "  --dmr     All DMR BED file (optional, for combined DMRs)"
    echo "  --wp      Assign DMRs to features with precedence (TSS > ... > intergenic)"
    echo "  --out     Output directory (default: dmr_feature_gene_mapping_output)"
    exit 1
}

PRECEDENCE=false

# Argument parsing loop
while [[ $# -gt 0 ]]; do
    case "$1" in
    --hyper)
        HYPER_FILE="$2"
        shift 2
        ;;
    --hypo)
        HYPO_FILE="$2"
        shift 2
        ;;
    --dmr)
        DMR_FILE="$2"
        shift 2
        ;;
    --out)
        OUT_DIR="$2"
        shift 2
        ;;
    --wp)
        PRECEDENCE=true
        shift
        ;;
    *) # if no args, call usage and exit
        usage
        ;;
    esac
done

# Input validation
if [[ -z "$HYPER_FILE" && -z "$HYPO_FILE" && -z "$DMR_FILE" ]]; then
    echo "Error: You must specify at least one of --hyper, --hypo, or --dmr"
    usage
fi

# Create output directory
mkdir -p "$OUT_DIR"

echo "Starting DMR-feature-gene intersection..."

# Define function for one DMR type
# Intersects a given DMR file with all genomic features
# and maps the overlapping regions to genes.
process_dmr_type() {
    # $1: type label ("hyper", "hypo", or "all")
    local TYPE="$1"
    # $2: path to the DMR BED file
    local FILE="$2"

    for FEATURE in "${FEATURES[@]}"; do
        FEATURE_BED="$ANNOT_DIR/hg38_${FEATURE}.bed"

        INTERSECT_OUT="$OUT_DIR/${TYPE}_dmrs_in_${FEATURE}.bed"
        GENE_HITS_OUT="$OUT_DIR/${TYPE}_${FEATURE}_gene_hits.tsv"
        GENE_IDS_OUT="$OUT_DIR/${TYPE}_${FEATURE}_gene_ids.txt"

        echo " - ${TYPE^^} DMRs in ${FEATURE}..."
        # Intersect DMRs with the current genomic feature
        bedtools intersect -wa -wb -a "$FILE" -b "$FEATURE_BED" >"$INTERSECT_OUT"

        # Intersect the above result with gene coordinates
        bedtools intersect -wa -wb -a "$INTERSECT_OUT" -b "$GENE_BED" >"$GENE_HITS_OUT"

        # Extract unique Ensembl gene IDs from the last column of the result
        awk '{print $NF}' "$GENE_HITS_OUT" | sort | uniq >"$GENE_IDS_OUT"
    done
}

assign_dmr_exclusive() {
    local TYPE="$1"
    local FILE="$2"
    local REMAINING="$OUT_DIR/${TYPE}_remaining.bed"

    cp "$FILE" "$REMAINING"

    for FEATURE in "${FEATURES[@]}"; do
        FEATURE_BED="$ANNOT_DIR/hg38_${FEATURE}.bed"

        INTERSECT_OUT="$OUT_DIR/${TYPE}_dmrs_in_${FEATURE}.bed"
        GENE_HITS_OUT="$OUT_DIR/${TYPE}_${FEATURE}_gene_hits.tsv"
        GENE_IDS_OUT="$OUT_DIR/${TYPE}_${FEATURE}_gene_ids.txt"

        echo " - ${TYPE^^} DMRs in ${FEATURE} (exclusive)..."
        # Intersect only the remaining DMRs
        bedtools intersect -a "$REMAINING" -b "$FEATURE_BED" -u > "$INTERSECT_OUT"

        # Intersect with genes
        bedtools intersect -wa -wb -a "$INTERSECT_OUT" -b "$GENE_BED" > "$GENE_HITS_OUT"

        # Extract gene IDs
        awk '{print $NF}' "$GENE_HITS_OUT" | sort | uniq > "$GENE_IDS_OUT"

        # Remove assigned DMRs from remaining
        bedtools intersect -v -a "$REMAINING" -b "$FEATURE_BED" > "${REMAINING}.tmp"
        mv "${REMAINING}.tmp" "$REMAINING"
    done

    rm -f "$REMAINING"
}

if [[ -n "$HYPER_FILE" ]]; then
    if $PRECEDENCE; then
        assign_dmr_exclusive "hyper" "$HYPER_FILE"
    else
        process_dmr_type "hyper" "$HYPER_FILE"
    fi
fi

if [[ -n "$HYPO_FILE" ]]; then
    if $PRECEDENCE; then
        assign_dmr_exclusive "hypo" "$HYPO_FILE"
    else
        process_dmr_type "hypo" "$HYPO_FILE"
    fi
fi

if [[ -n "$DMR_FILE" ]]; then
    if $PRECEDENCE; then
        assign_dmr_exclusive "all" "$DMR_FILE"
    else
        process_dmr_type "all" "$DMR_FILE"
    fi
fi

echo "Done. Results saved in $OUT_DIR"
