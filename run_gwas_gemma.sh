# RUN GEMMA Association Tests with Univariate Linear Mixed Models (GWAS)

# module load plink/1.9-180913

## Starter VCF
# VCF=<vcf>

## Create PED files from VCF
# plink --vcf $VCF --make-bed --out <prefix> --double-id --allow-extra-chr

## change phenotypes in .fam file
# open .fam and in the last column change the number to 0 if it's a control or 1 if it's a case

### GEMMA
## Create relationship matrix
gemma -bfile <prefix> -gk 1 -o <prefix>_rel_matrix
gemma -bfile <prefix> -k output/<prefix>_rel_matrix.cXX.txt -lmm 4 -o gwas_gemma

## To plot the GWAS, first open the .assoc.txt file and change the name of the columns:
# chr = CHR
# ps = BP 
# p_score = P
# rs = SNP

### Make sure your CHR column contains only numbers
## in my case, the names look like this 
# scaffold1,26779259,r265z7927k5a0m100_f401Z26771332   
## so I need to remove "scaffold" and everything after the first comma
## remove scaffold:
# sed 's/scaffold//' -i output/gwas_gemma.assoc.txt
## remove everything after the first comma:
# sed 's/,[^[:blank:]]*//' -i output/gwas_gemma.assoc.txt

## use the R script plot_assoc.R to plot like
# Rscript plot_assoc.R output.pdf output/gwas_gemma.assoc.txt "Case vs control"
