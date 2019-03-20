#!/usr/bin/env python3
"""
This is a sub-class of the FeatureTable class, intended to produce properly
formatted EMBL.
"""

import logging

from Bio.SeqFeature import FeatureLocation

from .feature_table import FeatureTable
from .embl_header import EMBLHeader
from .feature import Feature

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
        self.add_source()

    def __repr__(self):
        output = str(self.header)
        output += self.feature_header()
        for feature in self.features:
            output += self.embl_feature(feature)
        output += self.embl_sequence_header()
        output += self.embl_sequence_data()
        output += "//"
        return output

    def add_source(self):
        """
        Adds the source feature to the record.
        """
        seq_length = self.header.settings['sequence_length']
        mol_type = self.header.settings['molecule_type']
        species = self.header.settings['species']
        if not seq_length:
            logging.warning(("Couldn't get sequence length, can't create "
                             "'source' qualifier."))
            return

        source = Feature.from_template("source")
        source.set_location(FeatureLocation(0, seq_length))
        source.set_qualifier('mol_type', mol_type)
        source.set_qualifier('organism', species)
        self.features[0:0] = [source]

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

        # The features do their own formatting, so this function is just a
        # wrapper.

        return str(feature)

    def embl_sequence_data(self):
        """
        3.4.18 The Sequence Data Line
        The sequence data line has a line code consisting of two blanks. The
        sequence is written 60 bases per line, in groups of 10 bases separated
        by a blank character, beginning at position 6 of the line. The direction
        listed is always 5' to 3', and wherever possible the non-coding strand
        (homologous to the message) has been stored. Columns 73-80 of each
        sequence line contain base numbers for easier reading and quick
        location of regions of interest. The numbers are right justified and
        indicate the number of the last base on each line.
        An example of a data line is:
            aaacaaacca aatatggatt [...] ctgtttgtta ttagctcatt        60

        The characters used for the bases correspond to the IUPAC-IUB
        Commission recommendations (see appendices).
        """
        counter = 0
        output = ""
        while counter < len(self.record.seq):
            line = str(self.record.seq[counter:counter+60])
            output += "     "
            output += (f"{line[:10]:<10} {line[10:20]:<10} {line[20:30]:<10} "
                       f"{line[30:40]:<10} {line[40:50]:<10} {line[50:60]:<10}")
            counter += len(line)
            output += f"{counter:>10}\n"
        return output

    def embl_sequence_header(self):
        """
        Prints the sequence in EMBL format.

        3.4.17  The SQ Line
        The SQ (SeQuence header) line marks the beginning of the sequence data
        and Gives a summary of its content. An example is:
            SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
        As shown, the line contains the length of the sequence in base pairs
        followed by its base composition.  Bases other than A, C, G and T are
        grouped together as "other". (Note that "BP" is also used for single
        stranded RNA sequences, which is not strictly accurate, but has been
        used for consistency of format.) This information can be used as a check
        on accuracy or for statistical  purposes. The word "Sequence" is present
        solely as a marker for readability.
        """
        seq_len = len(self.record.seq)
        num_a = self.record.seq.upper().count('A')
        num_c = self.record.seq.upper().count('C')
        num_g = self.record.seq.upper().count('G')
        num_t = self.record.seq.upper().count('T')
        num_other = seq_len - num_a - num_c - num_g - num_t
        header = (f"SQ   Sequence {seq_len} BP; "
                  f"{num_a} A; {num_c} C; {num_g} G; {num_t} T; "
                  f"{num_other} other;\n")

        return header

    @staticmethod
    def feature_header():
        """
        The FH (Feature Header) lines are present only to improve readability of
        an entry when it is printed or displayed on a terminal screen. The lines
        contain no data and may be ignored by computer programs. The format of
        these lines is always the same.
        """
        return "FH   Key             Location/Qualifiers\nFH\n"

    def update_locus_tags(self, counter=0):
        """
        Updates all locus tags to the final value instead of a placeholder. This
        is done in a separate step to ensure that the resulting EMBL file has a
        sequential list of locus tags, which would not be guaranteed by the
        threaded implementation.
        """
        for feature in self.features:
            logging.info("updating locus tags")
            for qualifier in feature.qualifiers:
                if qualifier == 'locus_tag':
                    counter += 1 # increment before update so that we start at 1
                    template = qualifier.value[0]
                    qualifier.value[0] = template.format(number=counter)
                    break

        return counter
