# Bioinfo_scripts
Misc scripts than can be useful for bioinformatics

> Not all scripts were created by myself. Some have been adapted from scripts I have found. Credit can be found within each script

Not all scripts are ready to be used as is. They should be adapted for your needs.

## [snippets.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/snippets.sh)

Miscellaneous one-liners for everyday use

## [collect_HISAT2_mapping_stats.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/collect_HISAT2_mapping_stats.sh)

Collect HISAT2 mapping statistics from multiple summary files located in the current directory and organized them into a table

## [create_busco_plot.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/create_busco_plot.py)

Create busco summary plot for several busco results as seen here <https://busco.ezlab.org/busco_userguide.html#companion-scripts>.

## [create_coverage_subset.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/create_coverage_subset.R)

Create coverage table for mafCoverage results

## [create_qualimap_summary.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/create_qualimap_summary.sh)

Create summary table from qualimap results for several samples

## [ensembl_rapidrelease_homologue_query.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/ensembl_rapidrelease_homologue_query.py)
Query the ensembl rapid release homologue gene page to identify homologues of a given list of genes

## [genestats.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/genestats.sh)

Calculates several gene structure measures based on a GFF annotation file

## [get_assembly_stats.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/get_assembly_stats.py)

Get assembly statistics of a given fasta file

## [get_sniffles_summary.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/get_sniffles_summary.sh)

Create summary table of Sniffles results

## [get_struct_var_vcf_subset.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/get_struct_var_vcf_subset.sh)

Create summary table from svim and svim-asm results

## [makeAGPfromFasta.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/get_struct_var_vcf_subset.sh)

Create AGP file from a fasta file

## [misc_Rfucntions.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/misc_Rfunctions.R)

Collection of reusable R functions

## [mummerplot_first30.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/mummerplot_first30.sh)

Create mummerplot of the first X chromosomes/scaffolds

## [mummerplot_top40.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/mummerplot_top40.sh)

Create mummerplot from list of chromosomes/scaffolds

## [plot_assoc.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/plot_assoc.R)

Generate manhattan plot from GWAS results

## [plot_interactive_pca.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/plot_interactive_pca.py)

Create interactive PCA plot

## [plot_sample_mapping_summary.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/plot_sample_mapping_summary.R)

Plot histogram of mapping stats: mean coverage, mapping rate and mean mapping quality.
Meant to work with the output of create_qualimap_summary.sh

## [plot_vcfstats_sample_summary.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/plot_vcfstats_sample_summary.R)

Plot histogram of num ref homozygous, num non ref hom, num heterozygous, ts/tv, indels from the bcftools stats output

## [qualimap_summary.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/qualimap_summary.sh)

Get qualimap summary table

## [remove_duplicates_fasta.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/remove_duplicates_fasta.py)

Finds duplicated sequences in a fasta file, keeps only one copy and combines the name of the duplicated scaffolds.

## [rename_fasta.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/rename_fasta.py)

Rename fasta sequences according to a two column file. First column is current name, second column is new name

## [reverse_complement_fasta.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/reverse_complement_fasta.py)

reverse complement the sequences in the fasta file whose names are in scaffs_to_convert.txt

## [RMout_to_bed.pl](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/RMout_to_bed.pl)
Read a repeat masker output file and write TEs in a bed file

## [run_gwas_gemma.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_gwas_gemma.sh)

Run and plot GWAS with Gemma: Run GEMMA Association Tests with Univariate Linear Mixed Models

## Looking for signatures of selection

### [run_scan.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_scan.sh) and [run_scan.](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_scan.r)

Convert data from input file to an object of class haploh.  
Compute iHH, iES and inES over a whole chromosome

### [run_ihs.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_ihs.sh) and [run_ihs.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_ihs.R)

Compute iHS (standardized ratio of iHH values of two alleles). Create manhattan plot with results

### [run_xpehh.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_xpehh.sh) and [run_xpehh.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/run_xpehh.R)

Compute XP-EHH (standardized ratio of iES of two populations). Create manhattan plot with results

## [samplot_SV.sh](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/samplot_SV.sh)

Plot SV larger than minsize from a SV vcf summary table

## [snakemake_rule.code-snippets](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/snakemake_rule.code-snippets)

Visual Studio Code snippet to create snakemake rule

## [sniffles_var_summary.R](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/sniffles_var_summary.R)

Get SV summary table from VCF with sniffles results

## [structural_var_summary.R](structural_var_summary.R)
Get summary of the structural variants (from SVIM) found in a given vcf to stdout

## [vcf2eigenstrat.py](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/vcf2eigenstrat.py)

Convert a vcf file to eigenstrat format. Removes multi-alleleic and indel sites
