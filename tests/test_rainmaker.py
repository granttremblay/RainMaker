import unittest

import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))
from rainmaker import rainmaker

class TestBasics(unittest.TestCase):

    def test_alive(self):
        response = rainmaker.alive()
        self.assertTrue(response == "I'm alive!")

if __name__ == '__main__':
	unittest.main()