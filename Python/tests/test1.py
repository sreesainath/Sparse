#!/usr/bin/env python

import unittest
import sys
sys.path.append("..")


import sparse.headers.wrapper as w


class TestSequenceFunctions(unittest.TestCase):
    def setUp(self):
        pass

    def test_constructor(self):
        s = w.SparseTensor()
	        
    def test_getitem(self):
        s = w.SparseTensor()
        s[2,2] = 5

        self.assertTrue( s[2,2] == 5 )


if __name__ == '__main__':
    unittest.main()


