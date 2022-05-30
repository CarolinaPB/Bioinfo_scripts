#!/bin/bash
#SBATCH --time=10-0:00:0
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --error=scan_%j.err
#SBATCH --output=scan_%j.out
#SBATCH --job-name=scan

# First step for the ihs or xpehh analysis

vcf=<vcf>
prefix=<prefix>

echo $vcf
echo $prefix

Rscript run_ihs.R --vcf=$vcf --output=$prefix
