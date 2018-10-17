#!/usr/bin/env python3
"""
EMBL writer for ENA data submission. This implementation is basically just the
documentation at ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt in
python form.

GFF convertion is based on specifications from:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

import time
import gzip
import logging
from Bio import SeqIO
from BCBio import GFF
from concurrent.futures import ThreadPoolExecutor

from gffemblconverter.feature_table.embl_writer import EmblWriter

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
    if args.gff_file.endswith(".gz"):
        infile = gzip.open(args.gff_file)
    else:
        infile = open(args.gff_file)

    if args.fasta.endswith(".gz"):
        infasta = gzip.open(args.fasta)
    else:
        infasta = open(args.fasta)

    seq_dict = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

    return {"gff_files":infile, "base_dict":seq_dict}

def feature_table_args(args):
    """Convenience function to extract command line arguments and format them
    to be properly passed to a FeatureTable.
    """

    return {"header":{"data_class":args.data_class}}

if __name__ == '__main__':

    import argparse

    PARSER = argparse.ArgumentParser(description=__doc__)

    # Positional arguments
    PARSER.add_argument("gff_file", help="Input gff-file.")
    PARSER.add_argument("fasta", help="Input fasta sequence.")

    # Feature table header information
    PARSER.add_argument("-d", "--data_class",
                        default=None,
                        help="Data class of the sample.",
                        choices=["CON", "PAT", "EST", "GSS", "HTC", "HTG",
                                 "MGA", "WGS", "TSA", "STS", "STD"])

    # Logging arguments
    PARSER.add_argument("-v", "--verbose",
                        action="count",
                        default=2,
                        help="Increase logging verbosity.")
    PARSER.add_argument("-q", "--quiet",
                        action="count",
                        default=0,
                        help="Decrease logging verbosity.")

    # Script behaviour arguments
    PARSER.add_argument("--shame",
                        action="store_true",
                        help="Suppress the shameless plug.")
    PARSER.add_argument("-t", "--num_threads",
                        type=int, default=1,
                        help="Number of threads to use for conversion")

    ARGS = PARSER.parse_args()

    logging.basicConfig(format=("%(asctime)s %(levelname)s %(module)s: "
                                "%(message)s"),
                        level=50-10*(ARGS.verbose+ARGS.quiet),
                        datefmt="%H:%M:%S")

    FEATURES = []

    if not ARGS.shame:
        print(SHAMELESS_PLUG)

    THREAD_POOL = None
    if ARGS.num_threads > 1:
        THREAD_POOL = ThreadPoolExecutor(max_workers=ARGS.num_threads)

    logging.info("Starting record parsing")
    for i, record in enumerate(GFF.parse(**gff_input(ARGS))):
        FEATURES += [EmblWriter(record, 
                                thread_pool=THREAD_POOL,
                                **feature_table_args(ARGS))]

    for i, feature in enumerate(FEATURES):
        while feature.get_progress() < 1.0:
            time.sleep(0.1)
        print(f"{i} - {feature}")

