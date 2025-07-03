#!/usr/bin/bash

wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz -P ./

# The input GTF file is compressed in .gz format
GTF_FILE="Homo_sapiens.GRCh38.114.gtf.gz"
EXON_MERGED_FILE="Homo_sapiens.GRCh38.114_exons_merged.bed"
INTRON_FILE="Homo_sapiens.GRCh38.114_introns.bed"
TSS_FILE="Homo_sapiens.GRCh38.114_tss_merged.bed"
THREE_PRIME_UTR_MERGED_FILE="Homo_sapiens.GRCh38.114_3utrs.bed"
FIVE_PRIME_UTR_MERGED_FILE="Homo_sapiens.GRCh38.114_5utrs.bed"
DOWNSTREAM_FILE="Homo_sapiens.GRCh38.114_downstream_1kb_merged.bed"
# everything that does not include the rest of the files above 
# (note: the exons or other features of IncRNA etc are also considered intergenic in this case)
INTERGENIC_FILE="Homo_sapiens.GRCh38.114_intergenic_regions.bed" 

# Count the different types of features in the GTF file
echo "Counting the different types of features in the input file:"
zcat "$GTF_FILE" | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn

# Extract exons and merge overlapping/adjacent regions
echo "Extracting and merging exons:"
zcat "$GTF_FILE" | \
grep "protein_coding" | \
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1, $4-1, $5}' | \
sortBed | \
mergeBed -i - > "$EXON_MERGED_FILE"

# # Extract gene regions and subtract exon regions to get introns
echo "Extracting introns by subtracting exons from gene regions:"
zcat "$GTF_FILE" | \
grep "protein_coding" | \
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1, $4-1, $5}' | \
sortBed | \
subtractBed -a stdin -b "$EXON_MERGED_FILE" > "$INTRON_FILE"

# # Let's intersect exons and introns (this should produce no output as exons and introns shouldn't overlap)
echo "Intersecting exons and introns to check (should produce no output):"
intersectBed -a "$EXON_MERGED_FILE" -b "$INTRON_FILE"

# echo "getting chromosome lengths file"
wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > genome_file.txt
grep -E '^([1-9]|1[0-9]|2[0-2]|X|Y|MT)[[:space:]]' genome_file.txt | sort -k1,1 > genome_file_canonical.txt
GENOME_FILE="genome_file_canonical.txt"

echo "Extracting TSS regions around from gene start pos (±200bp):"
zcat "$GTF_FILE" | \
grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' | \
grep "protein_coding" | \
awk 'BEGIN{OFS="\t"}
     $3=="gene" {
         chrom=$1;
         start=($7=="+" ? $4-201 : $5-200);
         end=($7=="+" ? $4+199 : $5+200);
         if (start < 0) start = 0;
         print chrom, start, end}' | \
bedtools sort -i - -g genome_file_canonical.txt | \
mergeBed -i - > "$TSS_FILE"

# Extract 5utrs and 3utrs
echo "Extracting and merging 5 prime utrs:"
zcat "$GTF_FILE" | \
grep "protein_coding" | \
grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' | \
awk 'BEGIN{OFS="\t";} $3=="five_prime_utr" {print $1, $4-1, $5}' | \
sortBed | \
mergeBed -i - > "$FIVE_PRIME_UTR_MERGED_FILE"

echo "Extracting and merging 3 prime utrs:"
zcat "$GTF_FILE" | \
grep "protein_coding" | \
grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' | \
awk 'BEGIN{OFS="\t";} $3=="three_prime_utr" {print $1, $4-1, $5}' | \
sortBed | \
mergeBed -i - > "$THREE_PRIME_UTR_MERGED_FILE"

# Extract Immediate downstream
echo "Extracting immediate downstream regions of protein coding genes"
zcat "$GTF_FILE" | \
grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' | \
grep "protein_coding" | \
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
         print chrom, start, end;}' | \
bedtools sort -i - -g genome_file_canonical.txt | \
mergeBed -i - > "$DOWNSTREAM_FILE"

echo "Combining genic regions (exons, introns, UTRs, promoters, downstream)"
cat "$EXON_MERGED_FILE" \
    "$INTRON_FILE" \
    "$TSS_FILE" \
    "$THREE_PRIME_UTR_MERGED_FILE" \
    "$FIVE_PRIME_UTR_MERGED_FILE" \
    "$DOWNSTREAM_FILE" | \
sort -k1,1 -k2,2n | \
mergeBed -i - > genic_regions_merged.bed

# Extract intergenic regions (complement of protein coding genes)
echo "Extracting intergenic regions (complement of genic features):"

grep -E '^([1-9]|1[0-9]|2[0-2])[[:space:]]' genic_regions_merged.bed | \
bedtools complement -i - -g "$GENOME_FILE" > "$INTERGENIC_FILE"