#!/usr/bin/env python

'''
Rainmaker maps the cooling-to-freefall time ratio as a function
of radius in Chandra X-ray observations of hot galaxy
cluster atmospheres. This first iteration uses the main data table
from the  ACCEPT sample: http://www.pa.msu.edu/astro/MC2/accept/

Projected radial X-ray temperature and density profiles
are fit in log space with 3rd-order polynomials. The logarithmic
pressure profile is then analytically differentiated to determine
the gravitational acceleration, from which rainmaker then derives
the freefall time. The cooling time is also computed from the
temperature profile.

Usage:
    $ python rainmaker.py [-f data_table.txt -n "Name of cluster"]

Example:
    $ python rainmaker.py
    This will run the full sequence using Abell 2597 as an example.

    $ python rainmaker.py -f accept_main_table.txt -n "Centaurus"
    $ python rainmaker.py -n "Abell 2151"
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
import matplotlib.style as style


def main():
    '''The main program runs the whole sequence.'''

    # Parse command line arguments. Iterate with user if cluster not found.
    filename, cluster_name = parse_arguments()

    # DATA is an astropy TABLE object,
    # filtered to show all properties of a given cluster
    # Can be split by e.g. data['Rin'], data['Mgrav'], etc.
    data = parse_data_table(filename, cluster_name)

    logTemp_fit(data)
    logPressure_fit(data)


def parse_arguments():
    '''Set up and parse command line arguments.'''

    parser = argparse.ArgumentParser(description=
                                    "Rainmaker fits ACCEPT profiles to quantify \
                                     parameters relevant to precipitation",
                                     usage="rainmaker.py -f table.txt -n name")

    parser.add_argument("-f", "--filename",
                        dest="filename",
                        required=False,
                        default="accept_main_table.txt",
                        help="input data table",
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


def fit_polynomial(data, ln_xray_property, deg, whatIsFit):
    '''
    Fits a DEG-order polynomial in x, y space.
    A 3rd order polynomial is a cubic function

    poly.polyfit() returns coefficients, from 0th
    order first to N-th order last (note that this is
    *opposite* from how np.polyfit behaves!).
    '''
    r, ln_r, r_fine, log10_r_fine, ln_r_fine = extrapolate_radius(data)

    print("Now fitting    |" + "  " + make_number_ordinal(deg) +
          " order polynomial to " + whatIsFit)

    coeffs = poly.polyfit(ln_r, ln_xray_property, deg)

    # polyval() is used to assemble cubic fit:
    # $p(x) = c_0 + c_1 x + c_2 x^2 + c3 x^3$
    # where c_n are the coeffs returned by polyfit()
    ln_fit = poly.polyval(ln_r, coeffs)
    fit = np.exp(ln_fit)

    # Now use these coefficients to extrapolate fit
    # across larger radius

    ln_fit_fine = poly.polyval(ln_r_fine, coeffs)
    fit_fine = np.exp(ln_fit_fine)

    return fit, r, fit_fine, r_fine


def extrapolate_radius(data):
    '''
    The ACCEPT radii are finite. Fix that.
    '''

    r = (data['Rin'] + data['Rout']) * 0.5
    ln_r = np.log(r.value)
    # this is the NATURAL logarithm, ln

    # Generate the radii you wish to extrapolate
    # across in log10 space
    log10_r_fine = np.arange(300.)/100. - 3.

    # Now un-log10 it, give it a unit
    r_fine = (10**log10_r_fine) * u.Mpc

    # Also give its unitless natural log, used for fitting
    # with polyval() and fit_polynomial()'s coefficients
    ln_r_fine = np.log(r_fine.value)

    return r, ln_r, r_fine, log10_r_fine, ln_r_fine


def logTemp_fit(data):
    '''
    Fit the logarithmic electron density profile ln(n_e) (in cm^-3)
    to a polynomial in log r (in Mpc) of degree 'deg'. Plot it.
    '''
    whatIsFit = "log temperature profile"
    deg = 3

    ln_t = np.log(data['Tx'].value)
    ln_terr = np.log(data['Txerr'] / data['Tx'])

    upperbound = data['Tx'] + data['Txerr']
    lowerbound = data['Tx'] - data['Txerr']

    fit, r, fit_fine, r_fine = fit_polynomial(data, ln_t, deg, whatIsFit)

    prettyplot()

    plt.figure()
    plt.plot(r.to(u.kpc), data['Tx'], marker='o', markersize=10, linestyle='None')

    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1,100)
    ax.set_ylim(1,10)

    plt.xlabel('Cluster-centric Radius (kpc)')
    plt.title('Projected X-ray Temperature')

    #plt.errorbar(r.to(u.kpc).value, data['Tx'].value, data['Txerr'].value)

    plt.fill_between(r.to(u.kpc).value, lowerbound.value, upperbound.value, facecolor='gray', alpha=0.5)
    plt.plot(r.to(u.kpc), fit)
    plt.plot(r_fine.to(u.kpc), fit_fine, linestyle='--')
    plt.show(block=True)


def logPressure_fit(data):
    '''
    Fit the logarithmic electron density profile ln(n_e) (in cm^-3)
    to a polynomial in log r (in Mpc) of degree 'deg'. Plot it.
    '''
    whatIsFit = "log pressure profile"
    deg = 3

    ln_p = np.log(data['Pitpl'].value)
    ln_perr = np.log(data['Perr'] / data['Pitpl'])

    upperbound = data['Pitpl'] + data['Perr']
    lowerbound = data['Pitpl'] - data['Perr']

    fit, r, fit_fine, r_fine = fit_polynomial(data, ln_p, deg, whatIsFit)

    prettyplot()

    plt.figure()
    plt.plot(r.to(u.kpc), data['Pitpl'], marker='o', markersize=10, linestyle='None')

    ax = plt.gca()
    # ax.set_yscale('log')
    ax.set_xscale('log')
    # ax.set_xlim(1,300)
    # ax.set_ylim(1,10)

    plt.xlabel('Cluster-centric Radius (kpc)')
    plt.ylabel('Projected X-ray Pressure')
    plt.title('Projected X-ray Pressure')

    plt.fill_between(r.to(u.kpc).value, lowerbound.value, upperbound.value, facecolor='gray', alpha=0.5)
    plt.plot(r.to(u.kpc), fit)
    plt.plot(r_fine.to(u.kpc), fit_fine, linestyle='--')
    plt.show(block=True)


def prettyplot():
    '''Plots should be pretty'''

#    plt.rcParams['font.family'] = 'sans-serif'
#    plt.rcParams['font.serif'] = 'Ubuntu'
#    plt.rcParams['font.monospace'] = 'Ubuntu Mono'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
#    #plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
#    plt.rcParams['legend.fontsize'] = 10
#    plt.rcParams['figure.titlesize'] = 12
#    plt.rcParams['axes.linewidth'] = 2

    style.use('ggplot')


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
