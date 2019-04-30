#! /usr/bin/env python

import argparse
import csv
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description = "This script removes false frequency codes")

    parser.add_argument("tags_file", help = 'the list of real telemetry frequencies and codes')
    parser.add_argument("data_file", help = 'telemetry data')

    return parser.parse_args() # Change to return instead of storing in a

def parse_tags():
    # codes dictionary: key = frequency, value = list of real codes
    codes = defaultdict(dict)

    # open the tags file
    with open(args.tags_file) as tags:

        #create a csv reader object
        reader = csv.reader(tags, delimiter='\t')

        # skip the header line
        header = next(reader)

        # read in the file, line by line
        for line in reader:

            # skip blank lines
            if not line:
                continue
            # otherwise, return the data we want by appending to a list
            else:
                # add a conditional to make sure the key exists
                if line[0] in codes:
                    codes[line[0]].append(line[1])
                else:
                    # define it as a list using []
                    codes[line[0]] = []
                    codes[line[0]].append(line[1])

        # check our work
        for freq,code in codes.items():
            print(freq, code)
    # Make sure to save the dictionary
    return codes


def parse_data(code_dict):
    # open, read, and parse telemetry data file
    with open(args.data_file, 'r') as data:
        for line in data:
            row = line.split()

            # skip the header line
            if row[0] == 'Date':
                # Keep the header for output data
                print(line)
                continue
            else:
                if row[5] in code_dict[row[4]]:
                    print(line)
                else:
                    continue


def main():
    # Save the output dictionary as code_dict
    code_dict = parse_tags()
    # Give the dictionary to parse_data
    parse_data(code_dict)

# get the argument before calling main
args = get_args()

# execute the program by calling main
if __name__=="__main__":
    main()
