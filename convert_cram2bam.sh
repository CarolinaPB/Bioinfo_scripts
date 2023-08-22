#!/bin/bash
#SBATCH --time=8:00:0
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --error=bam2cram_%j.err
#SBATCH --output=bam2cram_%j.out
#SBATCH --job-name=bam2cram
#SBATCH --mem=24G

for f in *.cram
do
    name=$(echo $f | sed 's/.cram//g')
if [[ ! -f $name.bam ]]
    then
        samtools view -b -@12 -o $name.bam $f
        samtools index -@12 $name.bam
    fi
done
