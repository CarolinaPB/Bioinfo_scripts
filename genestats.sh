#!/usr/bin/env bash
# from https://gist.github.com/darencard/fcb32168c243b92734e85c5f8b59a1c3
usage()
{
cat << EOF
genestats
Version 1.0 (10 January, 2018)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script calculates several gene structure measures based on a GFF annotation file.
The bgzip and tabix programs (http://www.htslib.org/doc/tabix.html) and bedtools
(https://bedtools.readthedocs.io/en/latest/) must be installed and in the \$PATH.
The user supplies the GFF3 file as a single argument to the command and the output is a
12 column tab-delimited text written to STDOUT with the following columns:
1. transcript ID
2. transcript sequence length
3. number of exons
4. total exon sequence length
5. number of introns
6. total intron sequence length
7. number of CDS chunks
8. total CDS sequence length
9. number of 5' UTR sequences
10. total 5' UTR sequence length
11. number of 3' UTR sequences
12. total 3' UTR sequence length

To work optimally and produce reliable results, the following tags in column 3 must be
included within the GFF file:
mRNA = transcript
exon = exons
CDS = coding sequence
five_prime_UTR = 5' UTR sequence
three_prime_UTR = 3' UTR sequence

The first three tags should be present in any GFF file but the latter two may not be
depending on the source of the file. In those cases, estimates of UTR features will be
erroneous, but the other features should be reported correctly. Note that temporary files
are created in the working directory as the program runs.

USAGE:
genestats <file.gff>

EOF
}

if [[ -z $1 ]] || [[ $1 == "-h" ]] || [[ $1 == "-help" ]]
then
	usage
	exit 1
fi

# this script basically pastes together 6 queries of the GFF file for each mRNA annotation
# for the mRNA, exons, introns, CDS, and UTRs (x2)
# for each, tabix is used to rapidly pull out the feature from the GFF file and
# the appropriate feature lines are then used to produce counts and sequence lengths
# for the introns, bedtools is used to essentially find the features not matching exons
# (i.e. inverse match)

cat $1 | sort -k1,1 -k4,4n -k5,5n | bgzip -c > $1.gz && \
tabix -p gff $1.gz && \
cat $1 | \
awk '{ if ($3 == "mRNA") print $0 }' | \
awk -F "ID=|;" '{ print $1, $2 }' | \
while read line; do \
parent=$("echo $line | awk '{ print $9 }'"); \
chrom=$("echo $line | awk '{ print $1 }'"); \
start=$("echo $line | awk '{ print $4 }'"); \
end=$("echo $line | awk '{ print $5 }'"); \
paste <(echo $parent) \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "mRNA") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "exon") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "mRNA") print $0 }' > mRNA.txt && \
	tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "exon") print $0 }' > exons.txt && \
	bedtools subtract -a mRNA.txt -b exons.txt | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print count, sum }' && \
	rm mRNA.txt exons.txt) \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "CDS") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "five_prime_UTR") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "three_prime_UTR") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += ($5 - $4 + 1) } END { print count, sum }'); \
done