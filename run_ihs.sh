#!/bin/bash
#SBATCH --time=3-0:00:0
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --error=ihs_%j.err
#SBATCH --output=ihs_%j.out
#SBATCH --job-name=ihs
#SBATCH --mem=80000

module load R/3.5.3

# Population 1
vcf=<yourvcf.vcf.gz>
prefix=<prefix>

time Rscript run_ihs.R --vcf=$vcf --output=$prefix
