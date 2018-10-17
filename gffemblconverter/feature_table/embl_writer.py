#!/usr/bin/env python3
"""
This is a sub-class of the FeatureTable class, intended to produce properly
formatted EMBL.
"""

from .feature_table import FeatureTable

class EmblWriter(FeatureTable):
    """
    EMBL formatting subclass of FeatureTable.

    This subclass will format all output as EMBL.
    """

    def __repr__(self):
        return "[EMBL: {}]".format(self.name)

    @staticmethod
    def as_embl(feature):
        """
        Formats a Feature object into a proper EMBL feature string.
        """
        return str(feature)
