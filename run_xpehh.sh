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
csv1=<wgscan_pop1.rds>
csv2=<wgscan_pop2.rds>
pop1=<pop1>
pop2=<pop2>

Rscript run_xpehh.R --csv1=$csv1 --csv2=$csv2 --pop1=$pop1 --pop2=$pop2