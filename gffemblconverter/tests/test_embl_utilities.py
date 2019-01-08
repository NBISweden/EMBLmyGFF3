#!/usr/bin/env python3
"""
Unit tests for the functions from embl_utilities.py.
"""

import unittest
from gffemblconverter.feature_table.embl_utilities import embl_line, \
    ensure_row_length, quoted

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
        # exact length
        self.assertEqual(ensure_row_length(("AC   here are some words that are "
                                            "exactly at the full length for a "
                                            "single line.")),
                         ("AC   here are some words that are exactly at the "
                          "full length for a single line.\n"))
        # long
        self.assertEqual(ensure_row_length(("AC   this is a long line of words"
                                            " that may or may not need to be "
                                            "truncated. Hint: it does.")),
                         ("AC   this is a long line of words that may or may "
                          "not need to be truncated. \nAC   Hint: it does.\n"))
        # other padding
        self.assertEqual(ensure_row_length(("AC            this is a long line "
                                            "of words that may or may not need "
                                            "to be truncated. Hint: it does."),
                                           pad=14),
                         ("AC            this is a long line of words that may "
                          "or may not need to be \nAC            truncated. "
                          "Hint: it does.\n"))
        # other separator
        self.assertEqual(ensure_row_length(("AC   Good times; Tests passing; "
                                            "Everything turns out great; No "
                                            "problems around here."),
                                           split_on="; "),
                         ("AC   Good times; Tests passing; Everything turns out"
                          " great; \nAC   No problems around here.\n"))

    def test_quoted(self):
        """Testing quoted
        """
        self.assertEqual(quoted('test string', '"'), '"test string"')
        self.assertEqual(quoted('"test string', '"'), '"test string"')
        self.assertEqual(quoted('test string"', '"'), '"test string"')
        self.assertEqual(quoted('"test string"', '"'), '"test string"')
        self.assertEqual(quoted('test string', '|'), '|test string|')
        self.assertEqual(quoted(''), '')

if __name__ == '__main__':
    unittest.main()
