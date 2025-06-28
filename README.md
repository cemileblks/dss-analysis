# Oxford Nanopore Methylation Data Analysis – DSS Pipeline

![Status](https://img.shields.io/badge/status-in--progress-yellow)
![License](https://img.shields.io/badge/license-TBD-lightgrey)
![R version](https://img.shields.io/badge/R-4.0+-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%2F%20Unix-lightgrey)
![Bioinformatics](https://img.shields.io/badge/domain-bioinformatics-brightgreen)
![Languages](https://img.shields.io/github/languages/top/cemileblks/dss-analysis)

## Description

This project is part of my MSc Bioinformatics project in undertanding DNA methylation patterns in mammalian cells using OXford Nanopore sequensing data. It explores whether methylation patterns differ significantly between wild-type and mutant mammalian cells lacking key DNA methyltransferases. The pipeline is built around the DSS package in R to:

- Detect differentially methylated loci (DMLs) and regions (DMRs)
- Filter for canonical chromosomes
- Evaluate biological relevance and statistical significance of findings

The motivation behind this work is to better understand how gene body methylation may influence gene expression. Long read sequencing offers a promising way to observe these patterns at high resolution across the entire genome. Through this project so far, I gained experience working with large methylation datasets, developing shell workflows, and performing statistical analyses using R.

Note: This project is still in progress.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [License](#license)

## Installation

Clone the repository:
```
git clone https://github.com/cemileblks/dss-analysis.git
cd dss-analysis
```
Install required software:

- R (v4.0+)
- R packages: `DSS`, `bsseq`, `GenomicRanges`, `rtracklayer`
- Command-line tools: `bedtools`, `awk`, `bash`
- IGV or UCSC Genome Browser (for visualisations)

## Usage

This pipeline takes `.bedMethyl` and `.bedGraph` files from Oxford Nanopore methylation calls and outputs both statistics (DMLs/DMRs) and genome browser compatible `.bed` files.

.bedMethyl -> format_for_dss.sh -> tidy_dssin.sh -> final_dss_analysis.R -> list of tested cpg loci.bed hyper_dmls.bed hypo_dmls.bed hyper_dmrs.bed hypo_dmrs.bed
           
           -> bedmethyl_to_cpg.bed.sh 
           (converst bedmethl to bed that
           covers the whole cpg site in genome browser)

.bedGraph files -> bedgraph_substract.sh         -> bedgraph_to_bigwig.sh (fixes c to CpG, sorts, bigwig file output (input_diff.bw)) -> IGV or genome browser visualisation to see if results are satisfactory, the other bedmehtyl_to_cpg also goes here as well as the detected hypo and hyper dmls and dmrs. The list of cpg loci from previous step is for statistical justification of the locations of dmls and dmrs in thge genome
                    (substracts two bedgraph files to obtain
                    differences in methylation levels and 
                    filters chromosomes to autosomal canonical chrs)


flowchart LR
    A[.bedMethyl] --> B[format_for_dss.sh]
    B --> C[tidy_dssin.sh]
    C --> D[final_dss_analysis.R]
    D --> E1[tested_cpg_loci.bed]
    D --> E2[hyper_dmls.bed / hypo_dmls.bed]
    D --> E3[hyper_dmrs.bed / hypo_dmrs.bed]

    A --> F[bedmethyl_to_cpg.bed.sh]
    F --> G[IGV-ready CpG BED format]

    H[.bedGraph WT + 3BKO] --> I[bedgraph_subtract.sh]
    I --> J[bedgraph_to_bigwig.sh]
    J --> K[input_diff.bw (IGV/Genome browser)]

## Credits

- [DSS R package documentation](https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html)  
  Used for differential methylation analysis in R.

- [CLASS 2023 tutorials](https://github.com/boulderrinnlab/CLASS_2023)  
  Helpful for understanding genomic overlaps and using `GenomicRanges`.

- Various YouTube walkthroughs  
  Especially those covering `GenomicRanges` usage and methylation data visualization in R.

Special thanks to my MSc supervisor and lab group for their guidance through this project.

## License

This project is shared for academic and educational purposes only and does not currently use a formal open-source license.  
If you plan to reuse or modify the code, please cite the project appropriately.