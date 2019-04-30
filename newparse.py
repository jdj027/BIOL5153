#! /usr/bin/env python

# For testing, use inputs watermelon.fsa and watermelon.gff
# In class way to solve parsegff problem - April 24, 2019
import argparse
from Bio import SeqIO
import csv
from collections import defaultdict
import re

def get_args():
    parser = argparse.ArgumentParser(description = 'Parses a GFF file')

    parser.add_argument("gff_file", help="GFF3-formatted file")
    parser.add_argument("fasta_file", help="FASTA file corresponding to the GFF3 file")

    return parser.parse_args()

def parse_fasta():
    genome_object = SeqIO.read(args.fasta_file, 'fasta')
    return genome_object.seq

def parse_gff(genome):
    # make a dictionary to hold cds sequence
    # key = exon # | value = sequence
    coding_seqs = defaultdict(dict)

    # open the tags file
    with open(args.gff_file, 'r') as gff:

        #create a csv reader object
        reader = csv.reader(gff, delimiter='\t')

        # read in the file, line by line
        for line in reader:

            # skip blank lines
            if not line:
                continue

            # skip commented lines
            # elif re.match("^#", line):
            #    continue

            # this line returns start/end coods for the feature
            else:
                begin = int(line[3])-1
                end = int(line[4])
                strand = line[6]
                feature_type = line[2]
                attributes = line[8]
                species = line[0]

                # Extract sequence for this feature
                feature_sequence = genome[begin:end]
                # Use to check for extracting the correct genome size
                # print(len(feature_sequence),line[5])

                # Reverse compliment the seq if necessary
                if(strand == '-'):
                    feature_sequence = rev_comp(feature_sequence)

                # calculate gc content
                gc_content = gc(feature_sequence)

                # Want to write out CDS sequences with exons
                if( feature_type == 'CDS'):
                    # split the attributes field into seperate parts to get gene info
                    exon_info = attributes.split(' ; ');

                    # extract the gene name
                    gene_name = exon_info[0].split()[1]

                    # test if there is an entry in index 2 - holds the value 'exon' for genes with introns
                    # if there is no value in index 2, the gene has no intron & we can just print it
                    if len(exon_info[0].split()) > 2:
                        # take the exon number
                        exon_number = exon_info[0].split()[-1]

                        if gene_name in coding_seqs:
                            # store the coding sequence
                            coding_seqs[gene_name][exon_number] = feature_sequence

                        else:
                            # if first time finding this gene, declare dictionary entry
                            coding_seqs[gene_name] = {}

                            # then store the coding seq like before
                            coding_seqs[gene_name][exon_number] = feature_sequence
                    else:
                        # print sequence in FASTA format
                        print('>' + line[0].replace('','_') + '_' + gene_name)
                        print(feature_sequence)

        # GFF file read, loop over codng_seqs to print CDS sequence
        for gene, exons in coding_seqs.items():
            # loop over all exons for this gene
            # SORT EXONS FIRST

            cds_for_this_gene = ''
            #print fasta header
		    print('>' + line[0].replace(' ', '_') + '_' + gene)


		    for exon_num, exon_seq in sorted(exons.items()):
			    # print(gene, exon_num, exon_seq)
			    cds_for_this_gene += exon_seq

		# print the CDS sequence for this gene
		print(cds_for_this_gene)


def rev_comp(seq):
    return seq.reverse_complement()
    #print('1 - ' + seq)
    #print('2 - ' + seq_revcomp)

def gc(seq):
    seq = seq.upper()
    count_of_G = seq.count('G')
    count_of_C = seq.count('C')

    return (count_of_G + count_of_C)/len(seq)

def main():
    genome = parse_fasta()
    parse_gff(genome)

# get the args before calling Main
args = get_args()

# Main code
if __name__=="__main__":
    main()
