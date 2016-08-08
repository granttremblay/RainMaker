#!/usr/bin/env python

'''
Rainmaker fits log density and temperature profiles to the ACCEPT
tables from Cavagnolo et al.

-Grant Tremblay (Yale University)
'''


import sys
import os
import argparse
import numpy as np
from astropy.io import ascii


def main():
    filename, cluster_name = parse_arguments()

    data = parse_data_table(filename, cluster_name)

    print(data)


def parse_arguments():
    '''Set up and parse command line arguments.'''

    parser = argparse.ArgumentParser(description="It's gonna rain...",
                                     usage="rainmaker.py -f data_table.txt")

    parser.add_argument("-f", "--filename", dest="filename", required=True,
                        help="input file with data", metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-n", "--name_of_cluster",
                        dest="name_of_cluster",
                        required=True,
                        help="Name of the cluster.")

    args = parser.parse_args()
    filename = args.filename.name

    # Be flexible with the cluster name.
    # If they entered a space, replace it with a '_',
    # then convert to UPPERCASE (to match ACCEPT table)

    cluster_name = args.name_of_cluster.replace(" ", "_").upper()

    return filename, cluster_name


def is_valid_file(parser, arg):
    '''Check to ensure existence of the file.'''

    if not os.path.isfile(arg):
        parser.error("Cannot find that data table: %s" % arg)
    else:
        print("\nSuccessfully found input file: %s" % arg)
        return open(arg, 'r')      # return an open file handle


def parse_data_table(filename, cluster_name):
    '''Match input cluster name to that in table, return that object's data'''

    data = ascii.read(filename)     # This creates a flexible Astropy TABLE

    # 'tcool5/2' is a bad column name. Change it if there.
    if 'tcool5/2' in data.columns:
        data.rename_column('tcool5/2', 'tcool52')

    # 'tcool3/2' is also a bad column name. Change it if there.
    if 'tcool3/2' in data.columns:
        data.rename_column('tcool3/2', 'tcool32')

    data, obs_by_name, clusters_in_table = search_for_cluster_name(data, cluster_name)


def search_for_cluster_name(data, cluster_name):
    '''Takes input astropy TABLE object'''

    obs_by_name = data.group_by('Name')
    clusters_in_table = obs_by_name.groups.keys

    cluster_found = cluster_name in clusters_in_table['Name']

    if cluster_found is True:
        print("Matched cluster name to one in the table: " + cluster_name)
        mask = data['Name'] == cluster_name
        masked_data = data[mask]
        print(masked_data)
        return masked_data, obs_by_name, clusters_in_table

    else:
        print("That cluster name wasn't found: " + cluster_name)
        cluster_name = str(input("Try again: ")).replace(' ','_').upper()
        cluster_found = cluster_name in clusters_in_table['Name']


    # elif any(clusters_in_table['Name'] == cluster_name):
    #     print("Matched given cluster name to one in the table: " + cluster_name)
    #     mask = data['Name'] == cluster_name
    #     masked_data = data[mask]
    #     print(masked_data)
    #     return masked_data, obs_by_name, clusters_in_table 


    # if not any(clusters_in_table['Name'] == cluster_name):
    #     print("That cluster name wasn't found: " + cluster_name)
    #     sys.exit("Not found.")


    # elif any(clusters_in_table['Name'] == cluster_name):
    #     print("Matched given cluster name to one in the table: " + cluster_name)
    #     mask = data['Name'] == cluster_name
    #     masked_data = data[mask]
        
    #     return masked_data, obs_by_name, clusters_in_table        


    # mask = data['Name'] == cluster_name
    
    # # Check if that mask is an array of only FALSE, which means the name wasn't found
    # if not any(mask):
    #     print("\n...but the cluster name was not found in the data table :( \n")
    #     sys.exit("Look here for appropriate names: http://www.pa.msu.edu/astro/MC2/accept/ \n")
    # else: 
    #     masked_data = data[mask]
    #     return masked_data


    

def alive():
    response = "I'm alive!"
    return response


if __name__ == '__main__':
    main()
    



