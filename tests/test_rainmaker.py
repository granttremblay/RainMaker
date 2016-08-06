#!/usr/bin/env python

'''
Basic unit testing for every rainmaker function. 
Run this either with py.test or `python test_rainmaker.py`
in the /tests/ directory. 

The test explicitly sets the PATH so that rainmaker
can be properly imported. 
'''

import unittest

import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))
from rainmaker import rainmaker

class TestBasics(unittest.TestCase):
    '''Test basic functionality to ensure the code is alive'''
    def test_alive(self):
        response = rainmaker.alive()
        self.assertTrue(response == "I'm alive!")

if __name__ == '__main__':
	unittest.main()