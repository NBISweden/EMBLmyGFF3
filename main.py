#!/usr/bin/env python3
"""
EMBL writer for ENA data submission. Note that this implementation is basically
just the documentation at ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
in python form.

GFF convertion is based on specifications from:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

import gzip
import logging
from Bio import SeqIO
from BCBio import GFF

from gffemblconverter.feature_table.feature_table import FeatureTable

SHAMELESS_PLUG = """
##############################################################################
# NBIS 2018 - Sweden                                                         #
# Authors: Martin Norling, Niclas Jareborg, Jacques Dainat                   #
# Please visit https://github.com/NBISweden/EMBLmyGFF3 for more information. #
##############################################################################

"""

def gff_input(args):
    """Convenience functions that opens the files supplied in the args
    structure in the appropriate way and returns them so that they can be
    supplied directly to a FeatureTable.
    """
    infile = gzip.open(args.gff_file) if args.gff_file.endswith(".gz") else open(args.gff_file)
    infasta = gzip.open(args.fasta) if args.fasta.endswith(".gz") else open(args.fasta)
    seq_dict = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

    return {"gff_files":infile, "base_dict":seq_dict}

if __name__ == '__main__':

    import argparse

    PARSER = argparse.ArgumentParser(description=__doc__)

    # Positional arguments
    PARSER.add_argument("gff_file", help="Input gff-file.")
    PARSER.add_argument("fasta", help="Input fasta sequence.")

    # Logging arguments
    PARSER.add_argument("-v", "--verbose",
                        action="count",
                        default=2,
                        help="Increase logging verbosity.")
    PARSER.add_argument("-q", "--quiet",
                        action="count",
                        default=0,
                        help="Decrease logging verbosity.")

    ARGS = PARSER.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)s %(module)s: %(message)s',
                        level=50-10*(ARGS.verbose+ARGS.quiet),
                        datefmt="%H:%M:%S")

    FEATURES = []

    logging.info("Starting record parsing")
    for i, record in enumerate(GFF.parse(**gff_input(ARGS))):
        FEATURES += [FeatureTable(record)]
