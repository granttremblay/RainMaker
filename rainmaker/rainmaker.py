#!/usr/bin/env python

import sys, os
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="It's gonna rain...",
                                usage="python rainmaker.py -f data_table.txt")

    parser.add_argument("-f", "--filename", dest="filename", required=True,
                    help="input file with data", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
    args = parser.parse_args()
    filename = args.filename.name

    dirtytable = parse_data_table(filename)
    print(dirtytable)
    


def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("Cannot find that data table: %s" % arg)
    else:
        print("Successfully found input file: %s" % arg)
        return open(arg, 'r')  # return an open file handle

def parse_data_table(intable):
    a = np.genfromtxt(intable, names=True)
    return a



def alive():
    response = "I'm alive!"
    return response


if __name__ == '__main__':
    main()
    



