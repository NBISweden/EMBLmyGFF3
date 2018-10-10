#!/usr/bin/env python3
"""This file includes the FeatureTable class and supporting functions.
"""

import os
import glob
import json
import logging

from .qualifier import Qualifier
from .feature import Feature

FILEPATH = os.path.dirname(os.path.abspath(__file__))
DATAPATH = os.path.abspath(
    os.path.join(FILEPATH, os.pardir, "data_definitions")
    )

class FeatureTable():
    """This is the parent class for the ENA formatter.

    This class is designed to keep track of the legal values and settings of a
    'feature table', as described in:
    http://www.insdc.org/files/feature_table.html
    """

    FEATURES_DIR = os.path.join(DATAPATH, "features")
    QUALIFIERS_DIR = os.path.join(DATAPATH, "qualifiers")
    LEGAL_XREF = os.path.join(DATAPATH, "legal_dbxref.json")

    DBXREF_CACHE = {}
    QUALIFIER_CACHE = {}
    FEATURE_CACHE = {}

    def __init__(self, record):
        self.name = None

        self.feature_templates = []
        self.qualifier_templates = []
        self.legal_dbxref = []

        if not self.legal_dbxref:
            self.load_dbxref(FeatureTable.LEGAL_XREF)
        if not self.qualifier_templates:
            self.load_qualifier_definitions(FeatureTable.QUALIFIERS_DIR)
        if not self.feature_templates:
            self.load_feature_definitions(FeatureTable.FEATURES_DIR)

        self.parse_record(record)

    def __repr__(self):
        return "[FeatureTable: {}]".format(self.name)

    def load_dbxref(self, filename):
        """Loads a list of legal database cross reference values.
        """
        if filename not in FeatureTable.DBXREF_CACHE:
            logging.info("Loading legal db xrefs")
            FeatureTable.DBXREF_CACHE[filename] = \
                json.loads(open(filename).read())
        self.legal_dbxref = FeatureTable.DBXREF_CACHE[filename]

    def load_feature_definitions(self, dirname):
        """Loads feature definitions from the supplied directory, and creates
        the self.feature_templates list. This list will be used when creating
        new features for the FeatureTable.
        """
        if dirname not in FeatureTable.FEATURE_CACHE:
            logging.info("Loading Feature definitions")
            logging.debug("Feature definition dir: %s", dirname)
            FeatureTable.FEATURE_CACHE[dirname] = []
            for filename in glob.glob(os.path.join(dirname, "*.json")):
                logging.debug(" - %s", os.path.basename(filename))
                FeatureTable.FEATURE_CACHE[dirname] += \
                    [Feature(json.loads(open(filename).read()))]
        self.feature_templates = FeatureTable.FEATURE_CACHE[dirname]

    def load_qualifier_definitions(self, dirname):
        """Loads qualifier definitions from the supplied directory, and
        creates the self.qualifier_templates list. This list will be used when
        creating new qualifiers for the FeatureTable.
        """
        if dirname not in FeatureTable.QUALIFIER_CACHE:
            logging.info("Loading Qualifier definitions")
            logging.debug("Qualifier definition dir: %s", dirname)
            FeatureTable.QUALIFIER_CACHE[dirname] = []
            for filename in glob.glob(os.path.join(dirname, "*.json")):
                logging.debug(" - %s", os.path.basename(filename))
                FeatureTable.QUALIFIER_CACHE[dirname] += \
                    [Qualifier(json.loads(open(filename).read()))]
        self.qualifier_templates = FeatureTable.QUALIFIER_CACHE[dirname]

    def parse_record(self, record):
        """Parses the record information into Features and Qualifiers of the
        FeatureTable.
        """
        self.name = record.name
