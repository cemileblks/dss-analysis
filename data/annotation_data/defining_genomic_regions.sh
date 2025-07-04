#!/usr/bin/bash

# This script defines exons, introns, TSS, UTRs, downstream, and intergenic regions from Ensembl GTF

# Parameters
RELEASE="114"
GENOME_NAME="hg38"
GTF_FILE="Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
FA_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FAI_FILE="${FA_FILE}.fai"
GENOME_FILE="genome_file.txt"
GENOME_FILE_CANONICAL="genome_file_canonical.txt"
GENOME_FILE_AUTOSOMAL="genome_file_autosomal.txt"
GTF_PROTEIN_CODING="${GENOME_NAME}_protein_coding_only.gtf.gz"

EXONS="${GENOME_NAME}_exons.bed"
INTRONS="${GENOME_NAME}_introns.bed"
TSS="${GENOME_NAME}_tss.bed"
UTR5="${GENOME_NAME}_5utrs.bed"
UTR3="${GENOME_NAME}_3utrs.bed"
DOWNSTREAM="${GENOME_NAME}_downstream_1kb.bed"
# everything that does not include the rest of the files above
# (note: the exons or other features of IncRNA etc are also considered intergenic in this case)
INTERGENIC="${GENOME_NAME}_intergenic.bed"
GENIC_MERGED="genic_regions_merged.bed"

# Download data
echo "Downloading Ensembl GTF (release $RELEASE)..."
wget -nc ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/homo_sapiens/${GTF_FILE} -P ./

echo "Downloading reference genome..."
wget -nc "ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/homo_sapiens/dna/${FA_FILE}.gz"
gunzip -c "${FA_FILE}.gz" > "$FA_FILE"

# Index FASTA and extract genome sizes
samtools faidx "$FA_FILE"
cut -f1,2 "$FAI_FILE" > "$GENOME_FILE"
grep -E '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)[[:space:]]' "$GENOME_FILE" | sort -k1,1 > "$GENOME_FILE_CANONICAL"
grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' "$GENOME_FILE" | sort -k1,1 > "$GENOME_FILE_AUTOSOMAL"
echo "----------------------------------------"

# Count the different types of features in the GTF file
echo "Counting the different types of features in the input file:"
zcat "$GTF_FILE" | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn
echo "----------------------------------------"

# Define genomic features
echo "Filtering for protein-coding features..."
zcat "$GTF_FILE" |
    grep "protein_coding" | gzip >"$GTF_PROTEIN_CODING"

echo "Extracting all protein-coding gene coordinates (gene_id only)..."
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t"}
     $3 == "gene" {
         match($0, /gene_id "([^"]+)"/, gid);
         if (gid[1] != "") {
             print $1, $4-1, $5, gid[1];
         }
     }' |
    sortBed >"${GENOME_NAME}_protein_coding_genes.bed"

echo "Extracting and merging exons:"
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t";} $3=="exon" {print $1, $4-1, $5}' |
    sortBed |
    mergeBed -i - >"$EXONS"

# Extract gene regions and subtract exon regions to get introns
echo "Extracting introns... (Subtracting exons from gene regions)"
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t";} $3=="gene" {print $1, $4-1, $5}' |
    sortBed |
    subtractBed -a stdin -b "$EXONS" >"$INTRONS"

# Let's intersect exons and introns (this should produce no output as exons and introns shouldn't overlap)
echo "Intersecting exons and introns to check (should produce no output):"
intersectBed -a "$EXONS" -b "$INTRONS"

echo "Extracting promoters (±200bp TSS)..."
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t"}
     $3=="gene" {
         chrom=$1;
         start=($7=="+" ? $4-201 : $5-200);
         end=($7=="+" ? $4+199 : $5+200);
         if (start < 0) start = 0;
         print chrom, start, end}' |
    bedtools sort -i - -g "$GENOME_FILE_CANONICAL" |
    mergeBed -i - >"$TSS"

echo "Extracting 5' UTRs..."
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t";} $3=="five_prime_utr" {print $1, $4-1, $5}' |
    sortBed |
    mergeBed -i - >"$UTR5"

echo "Extracting 3' UTRs..."
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t";} $3=="three_prime_utr" {print $1, $4-1, $5}' |
    sortBed |
    mergeBed -i - >"$UTR3"

echo "Extracting immediate downstream regions..."
zcat "$GTF_PROTEIN_CODING" |
    grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' |
    awk 'BEGIN{OFS="\t"}
     $3=="gene" {
         chrom=$1;
         strand=$7;
         if (strand == "+") {
             start = $5 + 1;
             end = $5 + 1000;
         } else {
             start = $4 - 1000 - 1;
             end = $4 - 1;
             if (start < 0) start = 0;
         }
         print chrom, start, end;}' |
    bedtools sort -i - -g "$GENOME_FILE_CANONICAL" |
    mergeBed -i - >"$DOWNSTREAM"

# Merge all genic regions
echo "Merging all genic regions..."
cat "$EXONS" "$INTRONS" "$TSS" "$UTR5" "$UTR3" "$DOWNSTREAM" |
    sort -k1,1 -k2,2n |
    mergeBed -i - >"$GENIC_MERGED"

# Define intergenic regions
echo "Extracting intergenic regions..."
bedtools complement -i "$GENIC_MERGED" -g "$GENOME_FILE_AUTOSOMAL" >"$INTERGENIC"

# Convert chromosome names to UCSC style
echo "Converting chromosome names to UCSC style (only if not already prefixed)..."
for f in hg38_*.bed; do
    awk 'BEGIN{OFS="\t"}
       {
         if ($1 ~ /^chr/) {
           print
         } else {
           $1 = ($1 == "MT" ? "chrM" : "chr"$1);
           print
         }
       }' "$f" >tmp && mv tmp "$f"
done

# Cleanup large files
echo "Cleaning up temporary files..."
rm -f "$FA_FILE" "$FAI_FILE" "$GTF_PROTEIN_CODING" "$GENIC_MERGED"

echo "----------------------------------------"
echo "Done. Output files:"
ls -1 hg38_*.bed
