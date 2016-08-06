#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
from astropy.io import ascii


def main():
    parser = argparse.ArgumentParser(description="It's gonna rain...",
                                     usage="python rainmaker.py -f data_table.txt")

    parser.add_argument("-f", "--filename", dest="filename", required=True,
                        help="input file with data", metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-n", "--name_of_cluster", dest="name_of_cluster", required=False,
                        default='ABELL_2597', help="Name of the cluster, as written in the table")

    args = parser.parse_args()
    filename = args.filename.name

    # Be flexible with the cluster name. If they entered a space,
    # replace it with a '_', then convert it to UPPERCASE (as required by the ACCEPT table)
    clustername = args.name_of_cluster.replace(" ", "_").upper()

    data = parse_data_table(filename, clustername)
    print(data)
    


def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("Cannot find that data table: %s" % arg)
    else:
        print("Successfully found input file: %s" % arg)
        return open(arg, 'r')      # return an open file handle




def parse_data_table(intable, cluster_name):
    data = ascii.read(intable)     # This creates an Astropy TABLE object
                                   # http://docs.astropy.org/en/v1.2.1/table/index.html#astropy-table

    mask = data['Name'] == cluster_name
    
    # Check if that mask is an array of only FALSE, which means the name wasn't found
    if not any(mask):
        sys.exit("... but the cluster name was not found in the data table.")
    else: 
        masked_data = data[mask]
        return masked_data


def alive():
    response = "I'm alive!"
    return response


if __name__ == '__main__':
    main()
    



