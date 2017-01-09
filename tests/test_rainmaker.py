#!/usr/bin/env python

'''
===================================
Basic unit testing for rainmaker
===================================


Usage:  `py.test` or `python test_rainmaker.py` in the /tests/ directory.

:Author: Dr. Grant R. Tremblay (Yale University)
:Date: August 2016

The test explicitly sets the PATH so that rainmaker
can be properly imported.

http://docs.python-guide.org/en/latest/writing/structure/
'''

import os
import sys
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const

import numpy as np

import unittest

# Set the path explicitly #
sys.path.insert(0, os.path.abspath(__file__+"/../.."))
from rainmaker import rainmaker



class TestBasics(unittest.TestCase):
    '''Test basic functionality to ensure the code is alive'''

    def test_parse_data_table(self):
        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"

        returned_data = rainmaker.parse_data_table(filename, cluster_name)
        self.assertTrue(cluster_name in returned_data['Name'])

    def test_filter_by_cluster(self):

        '''Filtering by cluster name should work'''
        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"
        data = ascii.read(filename)

        masked_data = rainmaker.filter_by_cluster(data, cluster_name)
        self.assertTrue(cluster_name in masked_data['Name'])
        return masked_data

    def test_assign_units(self):

        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        data = rainmaker.parse_data_table(filename, "ABELL_2597")

        massValue = data['Mgrav'][0]
        radiusValue = data['Rin'][0]

        self.assertTrue('solMass' is str(massValue.unit))

        Mpc_divides_correctly = radiusValue.unit / (u.pc * 1.0e6).to(u.Mpc)
        self.assertTrue(Mpc_divides_correctly == 1)

    def test_fit_polynomial(self):
        
        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"
        data = ascii.read(filename)
        data = rainmaker.parse_data_table(filename, "ABELL_2597")

        ln_t = np.log(data['Tx'].value)
        ln_terr = np.log(data['Txerr'] / data['Tx'])

        fit, r, fit_fine, r_fine, temp_coeffs = rainmaker.fit_polynomial(data, ln_t, deg=3, whatIsFit="Test")
        self.assertTrue(len(temp_coeffs) == 4)
        # THIS IS A DUMB TEST. FIX IT! 

    def test_make_number_ordinal(self):

        number = rainmaker.make_number_ordinal(3)
        self.assertTrue('3rd' == number)

        number = rainmaker.make_number_ordinal(28)
        self.assertTrue('28th' == number)


if __name__ == '__main__':
    unittest.main()