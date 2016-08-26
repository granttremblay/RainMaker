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

import unittest

# Set the path explicitly #
sys.path.insert(0, os.path.abspath(__file__+"/../.."))
from rainmaker import rainmaker



class TestBasics(unittest.TestCase):

    '''Test basic functionality to ensure the code is alive'''
    def test_alive(self):
        response = rainmaker.alive()
        self.assertTrue(response == "I'm alive!")

    def test_parse_data_table(self):
        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"

        returned_data = rainmaker.parse_data_table(filename, cluster_name)
        self.assertTrue(cluster_name in returned_data['Name'])

    def test_filter_by_cluster(self):

        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"
        data = ascii.read(filename)

        masked_data = rainmaker.filter_by_cluster(data, cluster_name)
        self.assertTrue(cluster_name in masked_data['Name'])

    def test_assign_units(self):

        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        data = rainmaker.parse_data_table(filename, "ABELL_2597")

        massValue = data['Mgrav'][0]
        radiusValue = data['Rin'][0]

        self.assertTrue('solMass' is str(massValue.unit))

        Mpc_divides_correctly = radiusValue.unit / (u.pc * 1.0e6).to(u.Mpc)
        self.assertTrue(Mpc_divides_correctly == 1)





if __name__ == '__main__':
    unittest.main()