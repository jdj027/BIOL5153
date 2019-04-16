#! /usr/bin/env python3

# inputs
#gff_file   = 'watermelon.gff'
#fasta_file = 'watermelon.fsa'

import argparse
#from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "Parses a .fasta file using the features saved in a .gff file")

# add positional arguments
parser.add_argument("fasta_file_load", help="Input the name of the desired .fasta file as a sting.")
parser.add_argument("gff_file_load", help='Input the name of the desired .gff file as a sting.')

# parse the command line
args = parser.parse_args()

print('Successfully loaded ' + str(args.gff_file_load) + '.')
print('Successfully loaded ' + str(args.fasta_file_load) + '.')

# open the FASTA file
fasta = open(args.fasta_file_load, 'r')
next(fasta)
genome = fasta.read().rstrip("\n")

# open the GFF file
gff = open(args.gff_file_load, 'r')

# Read GFF file
for line in gff:
    # Remove newlines
    line = line.rstrip('\n')
    # Extract gff fields
    # [sequence, source, feature, begin, end, length, strand, phase, attributes]
    fields = line.split('\t')
    start = int(fields[3])-1
    stop = int(fields[4])-1

    # extract the DNA sequence from the genome for this feature
    substring = genome[start:stop]

    # print the DNA sequence for this feature
    print('DNA Sequence of feature:')
    print(substring)

    # print('Reverse Complement of feature:')
    # print(substring.reverse_complement.seq)

    # calculate the GC content for this feature, and print it to the screen
    g_content = substring.count('G')
    c_content = substring.count('C')
    sub_len = len(substring)
    gc_content = ((g_content + c_content)*100) / sub_len
    print('GC content of this feature: ' + str(gc_content) + '%')

fasta.close()
gff.close()
