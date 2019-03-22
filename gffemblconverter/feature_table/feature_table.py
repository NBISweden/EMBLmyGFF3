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
TRANSLATIONPATH = os.path.abspath(
    os.path.join(FILEPATH, os.pardir, "translations")
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
    TRANSLATION_DIR = TRANSLATIONPATH

    CACHE_DIRS = []
    DBXREF_CACHE = {}
    QUALIFIER_CACHE = {}
    FEATURE_CACHE = {}
    TRANSLATIONS = {}

    def __init__(self, record, thread_pool=None, header=None):
        self.name = None
        self.header = header
        self.thread_pool = thread_pool
        self.progress = [0, -1]
        self.record = record

        self.features = []
        self.legal_dbxref = []

        if not self.legal_dbxref:
            self.load_dbxref(FeatureTable.LEGAL_XREF)
        self.load_qualifier_definitions(FeatureTable.QUALIFIERS_DIR)
        self.load_feature_definitions(FeatureTable.FEATURES_DIR)

        self._load_translations(FeatureTable.TRANSLATION_DIR)
        self.parse_record(record)

    def __repr__(self):
        return "[FeatureTable: {}]".format(self.name)

    @staticmethod
    def load_translation(translation_file, translations_dir=None):
        """
        Loads a translation file, either from the current working dir, or from
        the translation_dir. Note that a file in the current working dir will
        have priority even if a translations_dir is given.
        """
        try:
            file_path = os.path.join(os.getcwd(), translation_file)
            os.stat(file_path)
        except FileNotFoundError:
            if translations_dir is not None:
                file_path = os.path.join(translations_dir, translation_file)
                os.stat(file_path)
            else:
                raise

        return json.load(open(file_path))

    def _load_translations(self, translations_dir):
        """
        Loads translation files. These will either be loaded from the current
        working dir, or from the translations_dir.

        the loaded files are:

        - translation_gff_feature_to_embl_feature.json
        - translation_gff_attribute_to_embl_qualifier.json
        - translation_gff_other_to_embl_qualifier.json
        """
        if "qualifiers" not in FeatureTable.TRANSLATIONS:
            logging.info("Loading qualifier translations")
            FeatureTable.TRANSLATIONS["qualifiers"] = self.load_translation(
                "translation_gff_attribute_to_embl_qualifier.json",
                translations_dir
                )
            FeatureTable.TRANSLATIONS["qualifiers"].update(
                self.load_translation(
                    "translation_gff_other_to_embl_qualifier.json",
                    translations_dir
                )
            )
            # Also add translations to Feature
            Feature.TRANSLATIONS["qualifiers"] = \
                FeatureTable.TRANSLATIONS["qualifiers"]
        if "features" not in FeatureTable.TRANSLATIONS:
            logging.info("Loading feature translations")
            FeatureTable.TRANSLATIONS["features"] = self.load_translation(
                "translation_gff_feature_to_embl_feature.json", translations_dir
                )
            # Also add translations to Feature
            Feature.TRANSLATIONS["features"] = \
                FeatureTable.TRANSLATIONS["features"]

    def get_progress(self):
        """
        Returns the prorgess of the current record import.
        """
        return self.progress[0]/self.progress[1]

    def new_feature(self, feature):
        """Creates a new Feature from the feature_templates dictionary.
        """
        template = Feature.from_seq_feature(feature, self.header.settings)
        self.progress[0] += 1

        return template

    def load_dbxref(self, filename):
        """Loads a list of legal database cross reference values.
        """
        if filename not in FeatureTable.DBXREF_CACHE:
            logging.info("Loading legal db xrefs")
            FeatureTable.DBXREF_CACHE[filename] = \
                json.loads(open(filename).read())
        self.legal_dbxref = FeatureTable.DBXREF_CACHE[filename]

    @staticmethod
    def load_feature_definitions(dirname):
        """Loads feature definitions from the supplied directory into the class
        variable FeatureTable.FEATURE_CACHE.
        This dictionary will be used when creating new features for the
        FeatureTable.
        """
        if dirname not in FeatureTable.CACHE_DIRS:
            logging.info("Loading feature definitions")
            logging.debug("Feature definition dir: %s", dirname)
            FeatureTable.CACHE_DIRS += [dirname]

            for filename in glob.glob(os.path.join(dirname, "*.json")):
                feature = Feature(json.loads(open(filename).read()))
                FeatureTable.FEATURE_CACHE[feature.name] = feature

        # Also add the feature templates to Feature
        Feature.FEATURE_TEMPLATES = FeatureTable.FEATURE_CACHE

    @staticmethod
    def load_qualifier_definitions(dirname):
        """Loads qualifier definitions from the supplied directory into the
        class variable FeatureTable.QUALIFIER_CACHE.
        This dictionary will be used when creating new qualifiers for the
        FeatureTable, and particularly it's Features.
        """
        if dirname not in FeatureTable.CACHE_DIRS:
            logging.info("Loading qualifier definitions")
            logging.debug("Qualifier definition dir: %s", dirname)
            FeatureTable.CACHE_DIRS += [dirname]

            for filename in glob.glob(os.path.join(dirname, "*.json")):
                qualifier = Qualifier(json.loads(open(filename).read()))
                FeatureTable.QUALIFIER_CACHE[qualifier.name] = qualifier

        # Also add the qualifier templates to Feature
        Feature.QUALIFIER_TEMPLATES = FeatureTable.QUALIFIER_CACHE

    def parse_record(self, record):
        """Parses the record information into Features and Qualifiers of the
        FeatureTable.
        """
        self.name = record.name
        self.progress[1] = len(record.features)
        for i, feature in enumerate(record.features):
            future = self.thread_pool.submit(self.new_feature, feature)
            try:
                self.features += [future.result()]
            except ValueError:
                # still count that we parsed the item
                self.progress[0] += 1
                pass
        logging.info("all jobs submitted")
