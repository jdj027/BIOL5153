#! /usr/bin/env python

# For testing, use inputs watermelon.fsa and watermelon.gff

import argparse
from Bio.Seq import Seq


def get_args():
    # create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
    parser = argparse.ArgumentParser(description = "Parses a .fasta file using the features saved in a .gff file")

    # add positional arguments
    parser.add_argument("fasta_file_load", help="Input the name of the desired .fasta file as a sting.")
    parser.add_argument("gff_file_load", help='Input the name of the desired .gff file as a sting.')

    # Return parsing variable
    return parser.parse_args()

def open_files(args):
# open the FASTA file
    fasta = open(args.fasta_file_load, 'r')
    genome = fasta.read().rstrip("\n")
    print('Successfully loaded ' + str(args.fasta_file_load) + '.')
    print('Total length of the genome: ' + str(len(genome)))

    # open the GFF file
    gff = open(args.gff_file_load, 'r')
    print('Successfully loaded ' + str(args.gff_file_load) + '.')

    return(gff,genome)

def parse_gff(line,genome):
    # Remove newlines
    line = line.rstrip('\n')
    # Extract gff fields
    # [sequence, source, feature, begin, end, length, strand, phase, attributes]
    fields = line.split('\t')
    start = int(fields[3])-1
    stop = int(fields[4])-1

    dna_sign = fields[6]
    # extract the DNA sequence from the genome for this feature
    substring = genome[start:stop]

    # print the DNA sequence for this feature
    print('DNA Sequence of feature:')
    print(substring)

    print('Directionality: ' + str(fields[6]))

    return(substring,dna_sign)

def rev_comp(substring,dna_sign):
    if dna_sign == '-':
        backward_dna = Seq(substring)
        print('Reverse Compliment: ')
        print(backward_dna.reverse_complement())

def gc_calc(substring):
    # calculate the GC content for this feature, and print it to the screen
    g_content = substring.count('G')
    c_content = substring.count('C')
    sub_len = len(substring)
    gc_content = ((g_content + c_content)*100) / sub_len
    print('GC content of this feature: ' + str(gc_content) + '%')

def main():
    # get the argument before calling main
    args = get_args()
    [gff,genome] = open_files(args)
    print(gff)
    for line in gff:
        [substring,dna_sign] = parse_gff(line,genome)
        rev_comp(substring,dna_sign)
        gc_calc(substring)

# execute the program by calling main
if __name__=="__main__":
    main()
