#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
from astropy.io import ascii

def main():
    filename, cluster_name = parse_arguments()

    data = parse_data_table(filename, cluster_name)

    print(data)

def parse_arguments():
    '''
    Use argparse voodoo to parse command line arguments. 
    It requires that the user input -f data_table.txt and -n name_of_cluster, 
    where name_of_cluster can flexibly be entered as, e.g., "abell 2597" or 
    "ABELL_2597". 

    If the given filename does not exist, or if the cluster name is not found 
    in the table, it exits with an error. 

    Calls is_valid_file()
    '''
    parser = argparse.ArgumentParser(description="It's gonna rain...",
                                     usage="python rainmaker.py -f data_table.txt")

    parser.add_argument("-f", "--filename", dest="filename", required=True,
                        help="input file with data", metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-n", "--name_of_cluster", dest="name_of_cluster", required=True,
                        help="Name of the cluster, e.g. 'Abell 2597', 'ZWICKY_2701', etc.")

    args = parser.parse_args()
    filename = args.filename.name

    # Be flexible with the cluster name. If they entered a space,
    # replace it with a '_', then convert it to UPPERCASE (as required by the ACCEPT table)
    cluster_name = args.name_of_cluster.replace(" ", "_").upper()

    return filename, cluster_name 


    

def is_valid_file(parser, arg):
    '''
    A simple check to ensure that the file given in the command line argument
    actually exists. If not, exit with an error (should probably improve this by 
    querying for a new filename instead of failing). 
    '''
    if not os.path.isfile(arg):
        parser.error("Cannot find that data table: %s" % arg)
    else:
        print("\nSuccessfully found input file: %s" % arg)
        return open(arg, 'r')      # return an open file handle


def parse_data_table(filename, cluster_name):
    data = ascii.read(filename)     # This creates an Astropy TABLE object
                                   # http://docs.astropy.org/en/v1.2.1/table/index.html#astropy-table

    mask = data['Name'] == cluster_name
    
    # Check if that mask is an array of only FALSE, which means the name wasn't found
    if not any(mask):
        print("\n...but the cluster name was not found in the data table :( \n")
        sys.exit("Look here for appropriate names: http://www.pa.msu.edu/astro/MC2/accept/ \n")
    else: 
        masked_data = data[mask]
        return masked_data





def alive():
    response = "I'm alive!"
    return response


if __name__ == '__main__':
    main()
    



