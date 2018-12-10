#!/usr/bin/env python3
"""
Unit tests for the functions from embl_utilities.py.
"""

import unittest
from gffemblconverter.feature_table.embl_utilities import embl_line, \
    ensure_row_length, ensure_quoted

class TestEMBLUtilities(unittest.TestCase):
    """
    Unit tests for all EMBL utility functions.
    """

    def test_embl_line(self):
        """Testing embl line function
        """
        self.assertEqual(embl_line("AC", "", True), "AC   \nXX\n")
        self.assertEqual(embl_line("AC", "", False), "AC   \n")
        self.assertEqual(embl_line("ACACAC", "", False), "ACAC \n")

    def test_ensure_row_length(self):
        """Testing ensure_row_length function
        """
        self.assertEqual(ensure_row_length(("AC   this is a long line of words"
                                            " that may or may not need to be "
                                            "truncated. Hint: it does.")),
                         ("AC   this is a long line of words that may or may "
                          "not need to be truncated.\nAC   Hint: it does.\n"))

    def test_ensure_quoted(self):
        """Testing ensure_quoted
        """
        self.assertEqual(ensure_quoted('test string', '"'), '"test string"')
        self.assertEqual(ensure_quoted('"test string', '"'), '"test string"')
        self.assertEqual(ensure_quoted('test string"', '"'), '"test string"')
        self.assertEqual(ensure_quoted('"test string"', '"'), '"test string"')
        self.assertEqual(ensure_quoted('test string', '|'), '|test string|')
        self.assertEqual(ensure_quoted(''), '')

if __name__ == '__main__':
    unittest.main()
