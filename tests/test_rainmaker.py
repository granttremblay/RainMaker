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

import unittest
import sys,os

# Set the path explicitly #
sys.path.insert(0,os.path.abspath(__file__+"/../.."))
from rainmaker import rainmaker



class TestBasics(unittest.TestCase):
    '''Test basic functionality to ensure the code is alive'''
    def test_alive(self):
        response = rainmaker.alive()
        self.assertTrue(response == "I'm alive!")

if __name__ == '__main__':
	unittest.main()