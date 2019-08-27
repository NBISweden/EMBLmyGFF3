#!/usr/bin/env python3
"""
Unit tests for the EmblLocation class from embl_writer.py.
"""

import unittest
from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
from EMBLmyGFF3.feature_table.feature import Feature

class TestEmblLocation(unittest.TestCase):
    """
    Unit tests for EmblLocation formatting.

    Examples are from http://www.insdc.org/files/feature_table.html#3.4.3.
    """

    def test_embl_location(self):
        """Testing embl location
        """
        test1 = FeatureLocation(466, 467, 1)
        test2 = FeatureLocation(339, 565, 1)
        test3 = FeatureLocation(BeforePosition(344), 500, 1)
        test4 = FeatureLocation(BeforePosition(0), 888, 1)
        test5 = FeatureLocation(0, AfterPosition(888), 1)

        # within, and between location is not yet supported as I don't know how
        # to use these types in a FeatureLocation to get the desired result.
        # test6 = FeatureLocation(WithinPosition(101,101,110), ??, strand=1)
        # test7 = BetweenPosition(123, 123, 124)
        test8 = FeatureLocation(11, 78, 1) + FeatureLocation(133, 202, 1)
        test9 = FeatureLocation(33, 126, -1)
        test10 = FeatureLocation(4917, 5163, -1) + \
                 FeatureLocation(2690, 4571, -1)

        self.assertEqual(Feature.embl_location(test1), "467")
        self.assertEqual(Feature.embl_location(test2), "340..565")
        self.assertEqual(Feature.embl_location(test3), "<345..500")
        self.assertEqual(Feature.embl_location(test4), "<1..888")
        self.assertEqual(Feature.embl_location(test5), "1..>888")
        # self.assertEqual(Feature.embl_location(test6), "102.110")
        # self.assertEqual(Feature.embl_location(test7), "123^124")
        self.assertEqual(Feature.embl_location(test8), "join(12..78,134..202)")
        self.assertEqual(Feature.embl_location(test9), "complement(34..126)")
        # Currently, biopython defaults to join(complement(), complement())
        self.assertEqual(Feature.embl_location(test10),
                         "complement(join(2691..4571,4918..5163))")
