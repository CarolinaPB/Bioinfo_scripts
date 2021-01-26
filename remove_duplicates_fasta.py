#######################################################################
# find duplicated sequences in a fasta file,                          #
# keep only one copy and combine the name of the duplicated scaffolds #
#######################################################################

## Usage:
# python remove_duplicates_fasta.py -f <file>

from Bio import SeqIO
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--File", help = "Fasta file to have duplicated sequences removed")
parser.add_argument("-o", "--Output", help = "Output file", default="deduplicated.fasta")

args = parser.parse_args()


dedup_records = defaultdict(list)
for record in SeqIO.parse(args.File, "fasta"):
    # Use the sequence as the key and then have a list of id's as the value
    dedup_records[str(record.seq)].append(record.id)
with open(args.Output, 'w') as output:
    for seq, ids in dedup_records.items():
        # Join the ids and write them out as the fasta
        output.write(">{}\n".format('|'.join(ids)))
        output.write(seq + "\n")
        