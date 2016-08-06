#!/usr/bin/env python

import sys, os
import argparse


def main():
    parser = argparse.ArgumentParser(description="It's gonna rain...",
                                usage="python rainmaker.py -f data_table.txt")

    parser.add_argument("-f", "--filename", dest="filename", required=True,
                    help="input file with data", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
    args = parser.parse_args()



def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("Cannot find that data table: %s" % arg)
    else:
        print("Successfully found input file: %s" % arg)
        return open(arg, 'r')  # return an open file handle



def parse_input_table(input_table):
    if os.path.isfile(input_table):
        print("Found the input file: " + input_table)
    else: 
        sys.exit("I c ")


def alive():
    response = "I'm alive!"
    return response


if __name__ == '__main__':
    main()
    



