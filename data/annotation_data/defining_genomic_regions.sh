#!/usr/bin/bash

# $INPUT_FILE=$1

# echo "How many different types of features are there in the input file?"

# zcat Homo_sapiens.GRCh38.114.gtf.gz | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn

# zcat Homo_sapiens.GRCh38.114.gtf.gz | \
# awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | \
# sortBed | mergeBed -i - > Homo_sapiens.GRCh38.114_exon_merged.bed


# zcat Homo_sapiens.GRCh38.114.gtf.gz | \
# awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | \
# sortBed | subtractBed -a stdin -b Homo_sapiens.GRCh38.114_exon_merged.bed | \
# gzip > Homo_sapiens.GRCh38.114_intron.bed
 
# #let's intersect the two files
# #this shouldn't produce any output
# intersectBed -a Homo_sapiens.GRCh38.114_exon_merged.bed -b Homo_sapiens.GRCh38.114_intron.bed

# The input GTF file is compressed in .gz format
GTF_FILE="Homo_sapiens.GRCh38.114.gtf.gz"
EXON_MERGED_FILE="Homo_sapiens.GRCh38.114_exons_merged2.bed"
INTRON_FILE="Homo_sapiens.GRCh38.114_introns2.bed"
INTERGENIC_FILE="Homo_sapiens.GRCh38.114_intergenic_regions.bed"

# Count the different types of features in the GTF file
# echo "Counting the different types of features in the input file:"
# zcat "$GTF_FILE" | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn

# echo "Filtering input file to genes and protein coding genes"
# # zcat "$GTF_FILE" |
# # grep "protein_coding" |


# Extract exons and merge overlapping/adjacent regions
# echo "Extracting and merging exons:"
# zcat "$GTF_FILE" | \
# grep "protein_coding" | \
# awk 'BEGIN{OFS="\t";} $3=="exon" {print $1, $4-1, $5}' | \
# sortBed | \
# mergeBed -i - > "$EXON_MERGED_FILE"

# # Extract gene regions and subtract exon regions to get introns
# echo "Extracting introns by subtracting exons from gene regions:"
# zcat "$GTF_FILE" | \
# grep "protein_coding" | \
# awk 'BEGIN{OFS="\t";} $3=="gene" {print $1, $4-1, $5}' | \
# sortBed | \
# subtractBed -a stdin -b "$EXON_MERGED_FILE" > "$INTRON_FILE"

# # Let's intersect exons and introns (this should produce no output as exons and introns don't overlap)
# echo "Intersecting exons and introns to check (should produce no output):"
# intersectBed -a "$EXON_MERGED_FILE" -b "$INTRON_FILE"

# echo "getting chromosome lengths file"
# wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
# cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > genome_file.txt
# grep -E '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)[[:space:]]' genome_file.txt | sort -k1,1 > genome_file_canonical.txt
GENOME_FILE="genome_file_canonical.txt"


# Extract intergenic regions (complement of protein coding genes)
echo "Extracting intergenic regions (complement of genes):"
zcat "$GTF_FILE" | grep -E '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)[[:space:]]' | \
grep "protein_coding" | \
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1, $4-1, $5}' | \
sortBed | \
complementBed -i stdin -g "$GENOME_FILE" > "$INTERGENIC_FILE"
