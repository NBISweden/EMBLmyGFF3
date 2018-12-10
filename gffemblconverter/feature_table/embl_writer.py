#!/usr/bin/env python3
"""
This is a sub-class of the FeatureTable class, intended to produce properly
formatted EMBL.
"""

from .feature_table import FeatureTable
from .embl_header import EMBLHeader

class EmblWriter(FeatureTable):
    """
    EMBL formatting subclass of FeatureTable.

    This subclass will format all output as EMBL.
    """

    def __init__(self, record, thread_pool=None, header=None):
        if header:
            header.sequence_length = len(record.seq)
            header.record_id = record.id
        super().__init__(record, thread_pool, EMBLHeader(header))

    def __repr__(self):
        output = str(self.header)
        return output

    @staticmethod
    def as_embl(feature):
        """
        Formats a Feature object into a proper EMBL feature string.
        """
        return str(feature)
