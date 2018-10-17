#!/usr/bin/env python3
"""This file includes the FeatureTable class and supporting functions.
"""

import os
import copy
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
    """This is the parent class for the EMBL formatter.

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

    def __init__(self, record, thread_pool=None, header=None):
        self.name = None
        self.header = header
        self.thread_pool = thread_pool
        self.progress = [0, -1]

        self.features = []
        self.templates = {"qualifiers":{},
                          "features":{}}
        self.legal_dbxref = []

        if not self.legal_dbxref:
            self.load_dbxref(FeatureTable.LEGAL_XREF)
        if not self.templates["qualifiers"]:
            self.load_qualifier_definitions(FeatureTable.QUALIFIERS_DIR)
        if not self.templates["features"]:
            self.load_feature_definitions(FeatureTable.FEATURES_DIR)

        self.parse_record(record)

    def __repr__(self):
        return "[FeatureTable: {}]".format(self.name)

    def get_progress(self):
        """
        Returns the prorgess of the current record import.
        """
        return self.progress[0]/self.progress[1]

    def insert_feature(self, feature, index=None):
        """Creates a new Feature from the feature_templates dictionary,
        and inserts it into the self.features list.
        """

        if feature.type not in self.templates["features"]:
            logging.error("Unknown Feature type: %s", feature.type)
            return

        template = copy.copy(self.templates["features"][feature.type])
        template.update_values(feature)

        if index:
            self.features[index] = [template]
        else:
            self.features += [template]
        self.progress[0] += 1

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
            FeatureTable.FEATURE_CACHE[dirname] = {}

            for filename in glob.glob(os.path.join(dirname, "*.json")):
                feature = Feature(json.loads(open(filename).read()))
                logging.debug(" - %s", feature.name)
                FeatureTable.FEATURE_CACHE[dirname][feature.name] = feature

        self.templates["features"] = FeatureTable.FEATURE_CACHE[dirname]

    def load_qualifier_definitions(self, dirname):
        """Loads qualifier definitions from the supplied directory, and
        creates the self.qualifier_templates list. This list will be used when
        creating new qualifiers for the FeatureTable.
        """
        if dirname not in FeatureTable.QUALIFIER_CACHE:
            logging.info("Loading Qualifier definitions")
            logging.debug("Qualifier definition dir: %s", dirname)
            FeatureTable.QUALIFIER_CACHE[dirname] = {}

            for filename in glob.glob(os.path.join(dirname, "*.json")):
                qualifier = Qualifier(json.loads(open(filename).read()))
                logging.debug(" - %s", qualifier.name)
                FeatureTable.QUALIFIER_CACHE[dirname][qualifier.name] = \
                    qualifier

        self.templates["qualifiers"] = FeatureTable.QUALIFIER_CACHE[dirname]

    def parse_record(self, record):
        """Parses the record information into Features and Qualifiers of the
        FeatureTable.
        """
        self.name = record.name
        self.progress[1] = len(record.features)
        for i, feature in enumerate(record.features):
            if self.thread_pool:
                self.features += [i]
                self.thread_pool.submit(self.insert_feature, (feature, i))
            else:
                self.insert_feature(feature)
