########################################################
# Doc to keep random snippets that can be useful again #
########################################################

# to get the fasta sequence of one scaffold
## get the lines between scaffolds_20 and scaffolds_21 (including headers)
sed -n '/Scaffolds_20 /,/Scaffolds_21/p' scaffolds.fasta > scaffolds_20.fasta
## remove all lines that match the pattern (>) 
## This leaves you with a file that has only sequence and no headers
grep -v ">" scaffolds_87.fasta > temp.txt; mv temp.txt scaffolds_87.fasta


# compare files to see if they are different
cmp --silent scaffolds_20.fasta scaffolds_87.fasta || echo "files are different"
Or 
cmp --silent  scaffolds_20.fasta scaffolds_87.fasta && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'

# delete the last line of a file
head -n -1 scaffolds_20.fasta > temp.txt ; mv temp.txt scaffolds_20.fasta

# delete first line of file
sed '1,1d' <file>

# delete last line of file
sed '$d' <file>

# To get length of fasta sequences
## index fasta
samtools faidx <file>
## get name of sequences followed by size
cut -f1-2 <file>.fai

# get slurm job info
sacct --format='JobID,JobName,Elapsed,MaxVMSize,State'

# show next X lines in grep search
grep -A10 <pattern> <file>

# show previous X lines in grep search
grep -B10 <pattern> <file>

# get top 40 scaffolds
sort -k2,2 -nr lengths_scaffolds.txt | head -n40

# rename fasta according to a two column file. first column is current name, second column is new name
python rename_fasta.py <fasta> --ids <chr file> > <out>

# sort fasta by size
conda activate random # Has bbmap sortbyname.sh
sortbyname.sh in=file.fa out=sorted.fa length descending

# Create bed file out of .fai
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' fasta.fai > fasta.bed

# hold job hpc
scontrol hold name=JOBNAME (or comma separated list of job ids)
scontrol release name=JOBNAME

# check predicted job start time
squeue --user=USERNAME --start

# check if file exists [ ! -f $f.bai ]
for f in a2*.RG.bam 
    do
        if [ ! -f $f.bai ]; then
            samtools index -@ 16 $f
        fi
    done

# split string at ., get first value
echo $f | awk -F'[_.]' '{print $1}'

# count number of reads mapped by bwa
grep "Processed" <BWA_LOG> | awk '{print $3}' | awk '{s+=$1} END {print s}'

# count number of reads in fasta.gz file
zgrep -c "^>" <reads.fasta.gz>

# Assembly stats (all give slightly different results)
conda activate py2.7
N50 -x <fasta>

or

conda activate random
statswrapper.sh <fasta1>,<fasta2> (for several files)
stats.sh in=<fasta>

or 

conda activate quast
quast --threads 16 --eukaryote --split-scaffolds --large <fasta> 
