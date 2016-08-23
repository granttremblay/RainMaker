#!/usr/bin/env python

'''
Rainmaker fits log density and temperature profiles to the ACCEPT
tables from Cavagnolo et al.

It works by fitting third-order polynomials in log space to the
temperature and pressure profiles, which is robust and usually
gives a pretty good fit.  The logarithmic pressure profile is
then  analytically differentiated to determine
rg(r) = - (kT / mu m_p) (d ln P / d ln r),
from which one gets the free-fall time.

In some instances, the best-fit pressure profile gets very flat
or even reverses near the center, which is handled by preventing
the derivative from going all the way to zero.  And to account
for the presence of a BCG, the minimum value of g is set to be
that of a singular isothermal sphere with a velocity dispersion
of 250 km/s.  This correction is important only within ~ 10 kpc,
if at all.

Notice that
t_c ~ kT / n \Lambda(T)  and  t_ff ~ r / (kT | dlnP / dlnr |)^1/2 ,
so their ratio is
t_c / t_ff ~ (1/nr) (kT)^3/2 [\Lambda(T)]^{-1}  | dlnP / dlnr |^1/2


At radii of several tens of kpc in a cool-core cluster,
the product of density times radius is roughly constant,
so the value of this ratio will track the temperature gradient.
However, we find that the density profile at smaller radii flattens
out in our deprojected profiles, which causes the product nr to
grow with radius and the timescale ratio to drop with radius.

So the key thing to look at is the behavior of the product of
radius and density at small radii.



-Grant Tremblay (Yale University)
'''

import sys
import os
import time
import argparse

import numpy as np
import scipy

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

    logTemp_fit(data)


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

    while not cluster_found:
            new_cluster_name = input("Cluster [" + cluster_name +
                                     "] not found, try again: ")
            if new_cluster_name.startswith('"') and new_cluster_name.endswith('"'):
                cluster_name = new_cluster_name[1:-1].replace(' ', '_').upper()
            else:
                cluster_name = new_cluster_name.replace(' ', '_').upper()
            cluster_found = cluster_name in clusters_in_table['Name']

    if cluster_found:
        print("Matched cluster name to one in the table: " + cluster_name)
        mask = data['Name'] == cluster_name
        masked_data = data[mask]
        print(masked_data)
        return masked_data


def fit_polynomial(x, y, deg, yerror, whatIsFit):
    '''
    Fits a DEG-order polynomial in x, y space.
    A 3rd order polynomial is a cubic function
    '''

    print("-----------------------------------------------------")
    print("Fitting " + make_number_ordinal(deg) +
          " order polynomial to " + whatIsFit)
    print("-----------------------------------------------------")

    coeffs, covariance = np.polyfit(x, y, deg, full=False, cov=True)
    #chi2 = np.sum((np.polyval(coeffs, x) - y)**2)

    print(coeffs)
    print(covariance)
    return coeffs

def logTemp_fit(data):
    '''
    Fit the logarithmic electron density profile ln(n_e) (in cm^-3)
    to a polynomial in log r (in Mpc) of degree 'deg'.
    The array 'coeffs' returns the coefficients of that polynomial fit.
    '''
    whatIsFit = "ln kT (keV) in log radius (Mpc) of degree DEG"

    deg = 3

    r = (data['Rin'] + data['Rout']) * 0.5
    logr = np.log(r)
    # this is the NATURAL logarithm, ln

    logt = np.log(data['Tx'])
    logterr = np.log(data['Txerr'] / data['Tx'])

    yerror = logterr

    coeffs = fit_polynomial(logr, logt, deg, yerror, whatIsFit)

    logtfit = 0.0 * logr

    for i in np.arange(deg):
        logtfit = logtfit + coeffs[i]*logr**i

    tfit = np.exp(logtfit)

    plt.scatter(r, data['Tx'], marker='o')
    plt.plot(r, tfit)
    plt.show(block=True)

    return coeffs


def plotter(x, y):

    plt.rcParams.update({'font.size': 22,
                         'axes.labelsize': 20,
                         'legend.fontsize': 16,
                         'xtick.labelsize': 18,
                         'ytick.labelsize': 18,
                         'axes.linewidth': 2})

    fig, ax = plt.subplots()
    ax.plot(x, y, 'b--',marker='s', label=r"$y = \alpha^2$")   


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
