#! /usr/bin/env python3

import argparse

def get_args():
    parser = argparse.ArgumentParser(description = "This script returns the Fibonacci sequence at a specified position")

    # one mandatory argument - num input
    parser.add_argument("num", type = int, help = "Input the desired position of the Fibonacci Sequence")

    # 2 optionals for how to print output as simple/verbose
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action = "store_true", help = 'Print verbose output') # Store_true makes default value false
    group.add_argument("-s", "--simple", action = "store_true", help = 'Print simple output (default)')

    return parser.parse_args() # Change to return instead of storing in a variable

def fib(n):
    # calculate Fibonacci
    a,b = 0,1

    for i in range(n):
        a,b = b, a+b;

    return(a)

def print_output(position, fib_num):
    if args.verbose:
        # Verbose output
        print("For Position " + position + "," + " The Fibonacci number is " + fib_num)
    else:
        # Simple output
        print(str(args.num),fib_num)

def main():
    fib_num = fib(args.num)
    print_output(args.num, fib_num)
    print(position)

# get the argument before calling main
args = get_args()

# execute the program by calling main
if __name__=="__main__":
    main()
