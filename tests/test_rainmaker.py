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
        data = ascii.read(filename)

        returned_data = rainmaker.parse_data_table(filename, cluster_name)
        self.assertTrue(cluster_name in returned_data['Name'])

    def test_filter_by_cluster(self):

        filename = os.getcwd() + "/testdata/accept_main_table.txt"
        cluster_name = "ABELL_2597"
        data = ascii.read(filename)

        masked_data = rainmaker.filter_by_cluster(data, cluster_name)
        self.assertTrue(cluster_name in masked_data['Name'])        
    



if __name__ == '__main__':
    unittest.main()