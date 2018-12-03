#!/usr/bin/env python3
"""
Unit tests for the functions from embl_header.py.
"""

import unittest
from gffemblconverter.feature_table.embl_header import EMBLHeader

class TestEMBLHeader(unittest.TestCase):
    """
    Unit tests for all EMBLHeader functions.
    """

    def setUp(self):
        self.header = EMBLHeader()

    def test_spacer(self):
        """Testing spacer function
        """
        self.assertEqual(self.header.spacer(), "XX\n")

    def test_accession(self):
        """Testing accession line function
        """
        self.assertEqual(self.header.accession_line(False), "AC   None;\nXX\n")
        self.assertEqual(self.header.accession_line(True), "AC * None;\nXX\n")

if __name__ == '__main__':
    unittest.main()
