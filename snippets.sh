########################################################
# Doc to keep random snippets that can be useful again #
########################################################

# to get the fasta sequence of one scaffold
## get the lines between scaffolds_20 and scaffolds_21 (including headers)
sed -n '/Scaffolds_20 /,/Scaffolds_21/p' scaffolds.fasta > scaffolds_20.fasta
## remove all lines that match the pattern (>) 
## This leaves you with a file that has only sequence and no headers
grep -v ">" scaffolds_87.fasta > temp.txt; mv temp.txt scaffolds_87.fasta

# Subset fasta file by chromosome name
samtools faidx <file>.fa <chr id>

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
sacct --format='JobID,JobName,Elapsed,MaxVMSize,State,Start,End'

--starttime YYYY-MM-DD[THH:MM[:SS]]
# show next X lines in grep search
grep -A10 <pattern> <file>

# show previous X lines in grep search
grep -B10 <pattern> <file>

# When grepping several files: show in which file the pattern was found in
grep -H <pattern> *.err

# find if variable is in file using grep x=exact matches=-F search in file. use -Fxq if grep output not needed
grep -Fx "$name" normal_BM_unique.txt

# grep to match exactly three underscores
grep -xiE '([^_]*_){3}[^_]*' 

# get top 40 scaffolds
sort -k2,2 -nr lengths_scaffolds.txt | head -n40

# rename fasta according to a two column file. first column is current name, second column is new name
python rename_fasta.py <fasta> --ids <chr file> > <out>

# sort fasta by size
conda activate random # Has bbmap sortbyname.sh
sortbyname.sh in=file.fa out=sorted.fa length descending

# Create bed file out of .fai
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' fasta.fai > fasta.bed

##### SLURM/HPC #####
# hold job hpc
scontrol hold name=JOBNAME (or comma separated list of job ids)
scontrol release name=JOBNAME

# check predicted job start time
squeue --user=USERNAME --start

# check job time limit
squeue -l -j <jobID>

# Show allocated resources and details of a job
scontrol show jobid <jobID>

# Get job efficiency data 
/cm/shared/apps/accounting/job_efficiency <JOBID>
seff <jobID>

# for interactive slurm session
sinteractive -c <num_cpus> --mem <amount_mem> --time <minutes>

# add to run programs compiled with an AVX512 instruction in it on SLURM
#SBATCH --constraint=avx512 

# add dependencies - run this script after JOBID has finished successfuly
#SBATCH  --dependency=afterok:<JOBID>:<JOBID>

# start job only jobs with JOBNAME have finished
sbatch --dependency=singleton --jobname=<JOBNAME> script.sh

# submit job that's not in a script
sbatch --job-name="<jobname>" --wrap="<command>"

# get file basename from file path (no extension)
$(basename $f) | sed "s/\..*//"

# get file basename from file path (with extension)
$(basename $f)

# get directory name above file
$(basename $(dirname $some_path))

# give permissions to a specific user
setfacl -m u:username:rwx myfolder

# Add permissions
chmod +rwx <filename>

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
# same for bwamem2
grep "Processed" <BWA_LOG> | awk '{print $4}' | awk '{s+=$1} END {print s}'


# mapping rate from bam file
module load bamtools
bamtools stats -in <bam> > <bam.stats>

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

# remove directories in the current directory that were created today
# print the directory names to confirm you want to delete them
find . -mtime -1 -type d -exec echo {} \;
# delete directories
find . -mtime -1 -type d -exec rm -r {} \;

# remove files in the current directory that were created today
# start by printing the files to confirm you want to delete them
find . -maxdepth 1 -mtime -1 -type f -name t* -print
# delete
find . -maxdepth 1 -mtime -1 -type f -name t* -delete

# remove everything after first space
sed 's/\s.*$//' <file>
## After first dot
sed 's/\..*$//'
## everything before first pipe
sed 's/.*|//'
## everything after the first two patterns (after the second underscore)
sed -E 's/((_*[^_]*){2}).*/\1/'
## everything after the first underscore
sed 's/_.*//'
## Remove comma followed by anything which is not a white-spaces and all matches.
sed 's/,[^[:blank:]]*//'

# Get number of each element in column 3
cut -f3 <file> | sort | uniq -c

# List number of directories in the current directory
ls -l | grep -c ^d

# To start jupyter lab from the HPC
conda activate py37  # has jupyter lab installed with conda
jupyter-lab --port 8888 --no-browser

#In local conda powershell:
ssh -t <user>@login.anunna.wur.nl -L 8888:localhost:8888

# time a command in jupyter notebook
%timeit <other code>

# Repeatmasker output to .bed
scripts/RMout_to_bed.pl <infile> <prefix>

# bed to gff3
conda activate genestructure # has genometools package
gt bed_to_gff3 <file> > <output file>

# unzip tar.gz
tar -xvf <file>

# untar tar.gz into specific directory
tar -zvxf file.tar --directory /path/to/directory

# create symbolic link
ln -s <directory/file to link> <name of link>

# list chromosomes in vcf file
zcat <file.vcf.gz> | grep "^[^#]" | cut -f 1 | uniq | sort -n


##### Snakemake ###################################################
# snakemake: create rule workflow pdf
snakemake --forceall --rulegraph | dot -Tpdf > workflow.pdf
snakemake --forceall --rulegraph | dot -Tpng > workflow.png
snakemake --forceall --dag | dot -Tpng > dag.png

# rerun rules with updated params
snakemake -n -R `snakemake --list-params-changes`

# rereun rules with updated code
snakemake -n -R `snakemake --list-code-changes`

# remove all snakemake output files
snakemake some_target --delete-all-output

# snakemake print only table with jobs to run
snakemake -np --quiet

# export conda env to environment.yaml
conda env export > environment.yml

# remove package from conda env
conda remove <package>

# create requirements.txt from conda env
conda list -e > requirements.txt

# install snakemake
conda install -c conda-forge mamba

mamba install -c conda-forge -c bioconda snakemake

##############################################################


# get number of homozygous SNPs in a vcf file 
zcat vcf1 | grep -v "#" | grep "1/1:" | wc -l

# get number of each SVTYPE in VCF (takes a while)
module load vcftools 
vcf-query -f  '%INFO/SVTYPE\n' <vcf> | sort |uniq -c

##Or use bcftools
for f in .vcf.gz
    do 
        echo $f
        bcftools query -f '[%INFO/SVTYPE\n]' $f | sort | uniq -c 
    done

# list directories
ls -d */

# get path to conda env
echo $CONDA_PREFIX

# get chosen fields out of VCF file based on file with CHR POS columns
## uniq is used since a line is printed for each sample, and those lines have the same info
bcftools query -R <chr_pos_file.txt> -f'[%CHROM\t%POS\t%INFO/CSQ\n]' <vcf> | uniq


# subset vcf to contain only one sample, don't include contigs or sex chromosomes
# Get number of SVs
bcftools view -S <file with one sample name> <vcf> | awk '$1 != "Contig"' | awk '$1 != "X"' | awk '$1 != "Y"' | awk '$1 != "MT"' | bcftools view -i 'AC>0' | perl -ne 'print "$1\n" if /[;\t]SVTYPE=([^;\t]+)/' | sort | uniq -c

# get sample names from VCF
bcftools query -l <vcf>

# subsample VCF based on text files with sample names
# use --force-samples if any of the samples is not in the original vcf
bcftools view --threads 16 -S <samples.txt> <vcf> | \
bcftools view -i 'AC>0' | \
bcftools +fill-tags -Oz > <output.vcf.gz>

# automatically install list of packages in R
list.of.packages <- c("optparse", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Git
## discard local changes to a certain file
git checkout <file>

# extract first sequence from fasta file
awk '/^>/{if(N)exit;++N;} {print;}' <file>

# format snakemake file (https://github.com/snakemake/snakefmt)
snakefmt Snakefile

# check what to do to follow best snakemake practices
snakemake --lint

# one line fasta
seqtk seq -l0 {input} > {output}

# Get specific sequence interval from fasta file and bed file (3 columns: chr   start   end)
bedtools getfasta -fi <file.fasta> -bed <file.bed> -fo <output.fa>

# concatenate several files with the same header. Keep the header
awk 'FNR>1 || NR==1' results/*_files.txt

# rsync from hpc to local machine
/usr/bin/rsync -av -e "ssh -v -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" --progress gw2hpct02:/hpc/my_group/my_username/data /tmp/data_from_HPC/

# rsync from local machine to HPC
/usr/bin/rsync -av -e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" --progress /tmp/data_to_HPC gw2hpct02:/hpc/my_group/my_username/data/ 

# rsync from storage to hpc
rsync --archive \
 --copy-links --progress \
 --log-file=rsync.log \
 --rsh "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" \
 "<name>":"/path/files/to/copy" \
 "/local/path"
 
 # rsync - copy files in file from directory
 rsync --archive \
 --copy-links --progress \
 --log-file=rsync.log \
 --rsh "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" \
 --files-from=files_to_sync.txt \
 "<name>":"/" \ # source
 "/local/path" #destination

# Automatically create repo changelog from git commits
git log --pretty="- %s" > CHANGELOG.md

# Check which step failed in cromwell pipeline
## returns all "rc" files with non-zero exit status
find -name rc -exec grep -vH 0 {} \;

# replace string in filename with another string (accepts filenames like filename*)
rename <string to replace> <new string> <filename>

# list files of tar.gz
tar -ztvf <file.tar.gz>

# check disk usage
du -sh .

# count number of lines in each file in directory
wc -l `find .  -type f`

# remove duplicate lines from several files
for f in *.tsv
    do
        name=$(echo $f | sed 's/.tsv//g')
        awk '!seen[$0]++' $f > $name.nodups.tsv
    done

# Count number of lines in (gzipped) files that don't contain "#". Include file name
zgrep  -Hcv -- '#' *.vcf.gz > <output_file>

# submit batch script that runs using container
## - create submit.sh that contains the SBATCH fields and the following line to 
## execute the script you want to run using the container
apptainer exec -B "<path/to/mount>" docker://<containerID> bash <yourscript.sh>

## - submit job with sbatch submit.sh

# Access existing job
srun --jobid=<JOBID> --pty /usr/bin/bash

# find files where the name contains certain pattern (can be in nested directories) and copy to another directory
# here, grep is used to make sure that we only take the "example" files from the <subdirectory> directory and not from
# another subdirectory where there could also be files containing "example"
find -type f -name "example*" | grep <subdirectory> | xargs -I {} cp {} /path/to/destination

# print out line 3 (or any other number) of file
sed -n "3p" <file>

# extract VEP column names from VCF header
zgrep "##INFO=<ID=CSQ" <file.vcf.gz> | grep -o 'Format:.*$' | sed 's/Format: //' | sed 's/|/,/g' | sed 's/">//'
