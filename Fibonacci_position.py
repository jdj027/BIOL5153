#! /usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description = "This script returns the Fibonacci sequence at a specified position")

# one mandatory argument - num input
parser.add_argument("num", type = int, help = "Input the desired position of the Fibonacci Sequence")

# 2 optionals for how to print output as simple/verbose
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action = "store_true", help = 'Print verbose output') # Store_true makes default value false
group.add_argument("-s", "--simple", action = "store_true", help = 'Print simple output (default)')

args = parser.parse_args()

# calculate Fibonacci
a,b = 0,1

for i in range(args.num):
    a,b = b, a+b;

if args.verbose:
    # Verbose output
    print("For Position " + str(args.num) + "," + " The Fibonacci number is " + str(a))
else:
    # Simple output
    print(args.num,a)
