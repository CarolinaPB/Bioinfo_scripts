from Bio.Seq import Seq
from Bio import SeqIO

#################################################################################################
# reverse complement the sequences in the fasta file whose names are in scaffs_to_convert.txt   #
#################################################################################################

# reverse complements sequences that are specified in a file
# for each sequence in the fasta file, if the name of the sequence is in the names file, the program will get it's reverse complement and write it to the new file. if it's not, it writes the original sequence to the fasta file


to_write = []

for record in SeqIO.parse("scaffolds.reduced.fa", "fasta"):
    with open('scaffs_to_convert.txt') as scaffs_to_convert:  
        if record.id in scaffs_to_convert.read():
           seq_to_convert = record.seq
           converted_seq =  seq_to_convert.reverse_complement()
           record.seq = seq_to_convert.reverse_complement()
           
           to_write.append(record)
           print(record.id)

        else:
            to_write.append(record)
            print()
            
            
SeqIO.write(to_write, "converted.fasta", "fasta")
            