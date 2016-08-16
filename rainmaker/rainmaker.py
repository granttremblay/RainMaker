#!/usr/bin/env python

'''
Rainmaker fits log density and temperature profiles to the ACCEPT
tables from Cavagnolo et al.

-Grant Tremblay (Yale University)
'''


import sys
import os
import time
import argparse
import numpy as np
from astropy.io import ascii
import astropy.units as u
import matplotlib.pyplot as plt


def main():
    '''The main program runs the whole sequence.'''

    # Parse command line arguments. Iterate with user if cluster not found.
    filename, cluster_name = parse_arguments()

    # DATA is an astropy TABLE object,
    # filtered to show all properties of a given cluster
    # Can be split by e.g. data['Rin'], data['Mgrav'], etc.
    data = parse_data_table(filename, cluster_name)

    logTemp(data)


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

    data = filter_by_cluster(data, cluster_name)

    return data


def filter_by_cluster(data, cluster_name):
    '''Takes input astropy TABLE object'''

    obs_by_name = data.group_by('Name')
    clusters_in_table = obs_by_name.groups.keys

    cluster_found = cluster_name in clusters_in_table['Name']

    while cluster_found is False:
            new_cluster_name = input("Cluster [" + cluster_name +
                                     "] not found, try again: ")
            if new_cluster_name.startswith('"') and new_cluster_name.endswith('"'):
                cluster_name = new_cluster_name[1:-1].replace(' ', '_').upper()
            else:
                cluster_name = new_cluster_name.replace(' ', '_').upper()
            cluster_found = cluster_name in clusters_in_table['Name']

    if cluster_found is True:
        print("Matched cluster name to one in the table: " + cluster_name)
        mask = data['Name'] == cluster_name
        masked_data = data[mask]
        print(masked_data)
        return masked_data


def fit_polynomial(x, y, deg, whatIsFit):
    '''Fits a DEG-order polynomial in x, y space'''

    print("-----------------------------------------------------")
    print("Fitting " + make_number_ordinal(deg) + " order polynomial to " + whatIsFit )
    print("-----------------------------------------------------")

    coeffs = np.polyfit(x, y, deg)

    chi2 = np.sum((np.polyval(coeffs, x) - y)**2)

    return coeffs, chi2


def logTemp(data):
    ''' 
    Fit the logarithmic electron density profile ln(n_e) (in cm^-3)
    to a polynomial in log r (in Mpc) of degree 'deg'.
    The array 'coeffs' returns the coefficients of that polynomial fit.
    '''
    whatIsFit = "log electron density (cm^-3) in log radius (Mpc)"

    nbins = len(data['Rin'])
    r = (data['Rin'] + data['Rout']) * 0.5
    logr = np.log10(r)

    lognelec = np.log10(data['nelec'])
    logneerr = np.log10(data['neerr'] / data['nelec'])
    yerror = logneerr

    deg = 3
    coeffs, chi2 = fit_polynomial(logr, lognelec, deg, whatIsFit)

    return coeffs


def plotter():

    plt.rcParams.update({'font.size': 22,
                         'axes.labelsize': 20,
                         'legend.fontsize': 16,
                         'xtick.labelsize': 18,
                         'ytick.labelsize': 18,
                         'axes.linewidth': 2})


def alive():
    response = "I'm alive!"
    return response


def make_number_ordinal(num):
    '''Take number, turn into ordinal. E.g., "2" --> "2nd" '''

    SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}

    if 10 <= num % 100 <= 20:
        suffix = 'th'
    else:
        # the second parameter is a default.
        suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix


def rainmaker_notebook_init(filename, cluster_name_raw):
    '''Run this in a Jupyter Notebook for exploration'''
    cluster_name = cluster_name_raw.replace(" ","_").upper()
    data = parse_data_table(filename, cluster_name)

    return data


if __name__ == '__main__':
    start_time = time.time()
    main()
    print("------ Finished in %s seconds -------"
          % (round((time.time() - start_time), 3)))
