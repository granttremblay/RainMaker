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

import os
import time
import argparse

import numpy as np
import numpy.polynomial.polynomial as poly

from astropy.io import ascii
from astropy.table import QTable
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

    parser = argparse.ArgumentParser(description=
                                     "Rainmaker fits ACCEPT profiles to quantify \
                                     parameters relevant to precipitation",
                                     usage="rainmaker.py -f table.txt -n 'Cluster'")

    parser.add_argument("-f", "--filename",
                        dest="filename", 
                        required=False,
                        default="accept_main_table.txt",
                        help="input file with data (default: accept_main_table.txt)",
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-n", "--name_of_cluster",
                        dest="name_of_cluster",
                        required=False,
                        default="Abell 2597",
                        help="Name of the cluster (default: Abell 2597)")

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
        parser.error("Cannot find that data table: {}".format(arg))
    else:
        print("\nTable found    |  {}".format(arg))
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

    data = assign_units(data)

    return data


def filter_by_cluster(data, cluster_name):
    '''Takes input astropy TABLE object'''

    obs_by_name = data.group_by('Name')
    clusters_in_table = obs_by_name.groups.keys

    cluster_found = cluster_name in clusters_in_table['Name']

    while not cluster_found:
            new_cluster_name = input("Cluster (" + cluster_name +
                                     ") not found, try again: ")

            messyname = (new_cluster_name.startswith('"') and
                         new_cluster_name.endswith('"'))

            # If the user entered quotation marks, strip them
            if messyname:
                cluster_name = new_cluster_name[1:-1].replace(' ', '_').upper()
            else:
                cluster_name = new_cluster_name.replace(' ', '_').upper()
            cluster_found = cluster_name in clusters_in_table['Name']

    if cluster_found:
        print("Cluster found  |  " + cluster_name)
        mask = data['Name'] == cluster_name
        masked_data = data[mask]
        return masked_data


def assign_units(data):

    # I could probably do this in a more intelligent manner,
    # but I want to assign units in a clear way!

    Name = data['Name']
    Rin = data['Rin'] * u.Mpc
    Rout = data['Rout'] * u.Mpc
    nelec = data['nelec'] * u.cm**(-3)
    neerr = data['neerr'] * u.cm**(-3)
    Kitpl = data['Kitpl'] * u.keV * u.cm**2
    Kflat = data['Kflat'] * u.keV * u.cm**2
    Kerr = data['Kerr'] * u.keV * u.cm**2
    Pitpl = data['Pitpl'] * u.dyne * u.cm**(-2)
    Perr = data['Perr'] * u.dyne * u.cm**(-2)
    Mgrav = data['Mgrav'] * u.M_sun
    Merr = data['Merr'] * u.M_sun
    Tx = data['Tx'] * u.keV
    Txerr = data['Txerr'] * u.keV
    Lambda = data['Lambda'] * u.erg * u.cm**3 / u.s
    tcool52 = data['tcool52'] * u.Gyr
    tcool52err = data['t52err'] * u.Gyr
    tcool32 = data['tcool32'] * u.Gyr
    tcool32err = data['t32err'] * u.Gyr

    names = ('Name', 'Rin', 'Rout', 'nelec', 'neerr', 'Kitpl',
             'Kflat', 'Kerr', 'Pitpl', 'Perr', 'Mgrav', 'Merr',
             'Tx', 'Txerr', 'Lambda', 'tcool52', 't52err',
             'tcool32', 't32err'
             )

    # Yes, I know I could do this in a for loop. But I want to
    # enable granular control over what columns are ultimately
    # written into the final "Science-ready" data table.

    # Note, this is an astropy QTable instead of a Table, so
    # that I can preserve units. Read more here:
    # http://docs.astropy.org/en/stable/table/mixin_columns.html#quantity-and-qtable
    data = QTable(
                 [Name, Rin, Rout, nelec, neerr, Kitpl,
                  Kflat, Kerr, Pitpl, Perr, Mgrav, Merr,
                  Tx, Txerr, Lambda, tcool52, tcool52err,
                  tcool32, tcool32err], names=names
                )

    return data


def fit_polynomial(x, y, deg, whatIsFit):
    '''
    Fits a DEG-order polynomial in x, y space.
    A 3rd order polynomial is a cubic function

    poly.polyfit() returns coefficients, from 0th
    order first to N-th order last (note that this is
    *opposite* from how np.polyfit behaves!).
    '''

    print("-----------------------------------------------------")
    print("Fitting " + make_number_ordinal(deg) +
          " order polynomial to " + whatIsFit)
    print("-----------------------------------------------------")

    coeffs = poly.polyfit(x, y, deg)
    logfit = poly.polyval(x, coeffs)

    fit = np.exp(logfit)

    return fit


def logTemp_fit(data):
    '''
    Fit the logarithmic electron density profile ln(n_e) (in cm^-3)
    to a polynomial in log r (in Mpc) of degree 'deg'.
    The array 'coeffs' returns the coefficients of that polynomial fit.
    '''
    whatIsFit = "ln kT (keV) in log radius (Mpc)"
    deg = 3

    r = (data['Rin'] + data['Rout']) * 0.5
    logr = np.log(r.value)
    # this is the NATURAL logarithm, ln

    logt = np.log(data['Tx'].value)
    logterr = np.log(data['Txerr'] / data['Tx'])

    fit = fit_polynomial(logr, logt, deg, whatIsFit)

    plt.figure()
    plt.plot(r.to(u.kpc), data['Tx'], marker='o', markersize=5, linestyle='None')

    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylim(1,10)

    #plt.fill_between(r, data['Tx']+data['Txerr'], data['Tx']-data['Txerr'], facecolor='blue', alpha=0.5)
    plt.plot(r.to(u.kpc), fit)
    plt.show(block=True)
 

def coolingFunction(kT):
    '''
    Implement the Tozzi & Norman (2001) cooling function,
    which is an analytic fit to Sutherland & Dopita (1993).

    This is shown in Equation 16 of Parrish, Quataert,
    & Sharma (2009), as well as Guo & Oh (2014).

    See here: arXiv:0706.1274. The equation is:

    $\Lambda(T) = [C_1 \left( \frac{k_B T}{\mathrm{keV}} \right)^{-1.7}
                  + C_2\left( \frac{k_B T}{\mathrm{keV}} \right)^{0.5}
                  + C_3] \times 10^{-22}$
    '''

    keV = u.eV * 1000.0

    # For a metallicity of Z = 0.3 Z_solar,
    C1 = 8.6e-3 * u.erg / (u.cm**3 * u.s)
    C2 = 5.8e-3 * u.erg / (u.cm**3 * u.s)
    C3 = 6.3e-2 * u.erg / (u.cm**3 * u.s)

    alpha = -1.7
    beta = 0.5

    coolingFunction = (
                       (C1 * (kT / keV)**alpha) +
                       (C2 * (kT / keV)**beta) +
                       (C3)
                       )*1e-22

    return coolingFunction


def plotter(x, y):

    plt.rcParams.update({'font.size': 22,
                         'axes.labelsize': 20,
                         'legend.fontsize': 16,
                         'xtick.labelsize': 18,
                         'ytick.labelsize': 18,
                         'axes.linewidth': 2})

    fig, ax = plt.subplots()
    ax.plot(x, y, 'b--', marker='s', label=r"$y = \alpha^2$")


def make_number_ordinal(number):
    '''Take number, turn into ordinal. E.g., "2" --> "2nd" '''

    suffixes = {1: 'st', 2: 'nd', 3: 'rd'}

    if 10 <= number % 100 <= 20:
        suffix = 'th'
    else:
        # the second parameter is a default.
        suffix = suffixes.get(number % 10, 'th')
    return str(number) + suffix


def rainmaker_notebook_init(filename, cluster_name):
    '''Run this in a Jupyter Notebook for exploration'''

    data = parse_data_table(filename, cluster_name.replace(" ", "_").upper())

    return data


if __name__ == '__main__':
    start_time = time.time()
    main()
    runtime = round((time.time() - start_time), 3)
    print("Finished in    |  {} seconds".format(runtime))
