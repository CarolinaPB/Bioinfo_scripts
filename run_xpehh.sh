#!/bin/bash
#SBATCH --time=3-0:00:0
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --error=xpehh_%j.err
#SBATCH --output=xpehh_%j.out
#SBATCH --job-name=xpehh
#SBATCH --mem=40000

module load R/3.5.3

# First run ihs

# compare pop1 with pop2
rds1=<wgscan_pop1.rds>
rds2=<wgscan_pop2.rds>
pop1=<pop1>
pop2=<pop2>

Rscript run_xpehh.R --rds1=$rds1 --rds2=$rds2 --pop1=$pop1 --pop2=$pop2