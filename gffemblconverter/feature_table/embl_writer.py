#!/usr/bin/env python3
"""
This is a sub-class of the FeatureTable class, intended to produce properly
formatted EMBL.
"""

from .feature_table import FeatureTable
from .embl_header import EMBLHeader
from .embl_utilities import embl_line

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
        output += self.feature_header()
        for feature in self.features:
            output += self.embl_feature(feature)
            break
        return output

    @staticmethod
    def as_embl(feature):
        """
        Formats a Feature object into a proper EMBL feature string.
        """
        return str(feature) + "\n"

    @staticmethod
    def embl_feature(feature):
        """
        The FT (Feature Table) lines provide a mechanism for the annotation of the
        sequence data. Regions or sites in the sequence which are of interest are
        listed in the table. In general, the features in the feature table represent
        signals or other characteristics reported in the cited references. In some
        cases, ambiguities or features noted in the course of data preparation have
        been included.  The feature table is subject to expansion or change as more
        becomes known about a given sequence.

        Feature Table Definition Document:
        A complete and definitive description of the feature table is given
        in the document "The DDBJ/ENA/GenBank Feature Table:  Definition".
        URL: ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.txt
        """
        line_code = "FT"
        information = str(feature)

        return embl_line(line_code, information, add_spacer=False)

    @staticmethod
    def feature_header():
        """
        The FH (Feature Header) lines are present only to improve readability of
        an entry when it is printed or displayed on a terminal screen. The lines
        contain no data and may be ignored by computer programs. The format of
        these lines is always the same.
        """
        return "FH   Key             Location/Qualifiers\nFH\n"
