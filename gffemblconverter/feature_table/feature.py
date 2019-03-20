#!/usr/bin/env python3
"""Feature object for EMBL feature tables
"""

import logging
from copy import deepcopy
from Bio.SeqIO import InsdcIO

from .embl_utilities import embl_line

class Feature():
    """Feature object for EMBL feature tables.

    From http://www.insdc.org/files/feature_table.html#3.2:
    "Feature keys indicate
    (1) the biological nature of the annotated feature or
    (2) information about changes to or other versions of the sequence.
    The feature key permits a user to quickly find or retrieve similar features
    or features with related functions. "
    """

    QUALIFIER_TEMPLATES = {}
    FEATURE_TEMPLATES = {}
    TRANSLATIONS = {}

    def __init__(self, definition):
        """
        Initializes a Feature from a defintion.
        """
        self.identifier = None
        self.location = None
        self.name = definition['feature_key']
        self.sub_features = []
        self.qualifiers = []
        self.optional_qualifiers = {}
        self.mandatory_qualifiers = {}
        for key, value in definition.items():
            setattr(self, key, value)

    def __repr__(self):
        information = f"{self.name:<16}{self.embl_location(self.location)}"

        output = embl_line("FT", information, add_spacer=False)
        for qualifier in self.qualifiers:
            output += f"{qualifier}"

        for sub_feature in self.sub_features:
            output += str(sub_feature)

        return output

    def add_locus_tag(self, value):
        """
        Checks if this feature needs a locus tag, and adds it if that is the
        case. To make sure that all locus tags are ordered, a placeholder tag
        will be used, which will be updated when the file is written to disk.
        """
        # Only add locus tag if it isn't already set
        if 'locus_tag' in self.qualifiers:
            return

        # Only set valid values
        if value is None:
            logging.warning('No value for locus tag.')

        # Only add if the feature supports locus tag
        for qualifiers in [self.optional_qualifiers, self.mandatory_qualifiers]:
            if 'locus_tag' in qualifiers:
                tag_template = str(value) + '_LOCUS{number}'
                self.set_qualifier('locus_tag', tag_template)

    @staticmethod
    def embl_location(location, rec_length=None):
        """
        Formats a FeatureLocation as an EMBL location by calling the
        InsdcIO._insdc_location_string function in biopython. I'd prefer to not
        call a private function here, but I haven't figured out how to properly
        do this in any other way.
        """
        if rec_length is None:
            rec_length = location.end
        return InsdcIO._insdc_location_string(location, rec_length)

    @staticmethod
    def from_template(template_name):
        """
        This function will create a new Feature from a template in the
        Feature.FEATURE_TEMPLATES dictionary, using the Feature.TRANSLATIONS
        dictionary to try to map to the correct feature type if the type is
        unknown.
        """
        if template_name not in Feature.FEATURE_TEMPLATES:
            # Check if we know a translation for this value
            translations = Feature.TRANSLATIONS.get("features", [])
            if template_name in translations:
                if "target" in translations[template_name]:
                    logging.debug("Translated feature %s to %s", template_name,
                                  translations[template_name]["target"])
                    template_name = translations[template_name]["target"]
                else:
                    logging.info("Feature %s has no translation target",
                                 template_name)
            else:
                logging.error("Unknown Feature type: '%s'", template_name)
                return None

        template = deepcopy(Feature.FEATURE_TEMPLATES[template_name])
        return template

    @staticmethod
    def from_seq_feature(feature, settings=None):
        """
        This template will then be updated with the values contained in feature,
        which should be a SeqFeature, and returned.
        """
        # Create feature from template
        template = Feature.from_template(feature.type)

        # Add values from SeqFeature
        template.update_values(feature)

        # Add values from settings
        if settings is not None:
            # See if we can set locus tag
            locus_tag = settings.get('locus_tag', None)
            if locus_tag is not None:
                template.add_locus_tag(locus_tag)

        return template

    def update_values(self, seq_feature):
        """
        Attempts to update all legal values in the Feature from the
        information in a Bio.SeqFeature.
        """

        # update own identifier to the seq_feature's id
        self.identifier = seq_feature.id

        # update own location to the seq_feature's location
        self.set_location(seq_feature.location)

        # add all legal qualifiers
        for qualifier, value in seq_feature.qualifiers.items():
            self.set_qualifier(qualifier, value)

        # set sub-features
        for sub_feature in seq_feature.sub_features:
            templated_feature = self.from_seq_feature(sub_feature)
            self.sub_features += [templated_feature]

    def set_location(self, seq_location):
        """
        Updates the features location with a seq_location.
        """
        self.location = seq_location

    def set_qualifier(self, qualifier, value, prefix=""):
        """
        Adds a legal qualifier value to the feature, given that it has a
        template in the self.optional_qualifiers or self.mandatory_qualifiers
        dictionaries.

        If the qualifier is unknown, the translations['qualifiers'] dictionary
        will be used to try to find a mapping to a legal qualifier.
        """
        if qualifier not in self.optional_qualifiers and \
           qualifier not in self.mandatory_qualifiers:
            translations = Feature.TRANSLATIONS.get('qualifiers', [])
            if qualifier not in translations:
                logging.warning(("Qualifier %s is neither an optional nor a "
                                 "mandatory qualifier of %s"),
                                qualifier, self.name)
                return
            if translations[qualifier].get("target", ""):
                prefix = translations[qualifier].get("prefix", prefix)
                logging.debug("Translated qualifier %s to %s", qualifier,
                              translations[qualifier]["target"])
                qualifier = translations[qualifier]["target"]
            else:
                logging.info("Qualifier %s has no translation target",
                             qualifier)
                return

        if qualifier not in Feature.QUALIFIER_TEMPLATES:
            logging.error("Legal qualifier %s not found in qualifier cache!",
                          qualifier)
            raise ValueError(f"Legal qualifier {qualifier} not found in cache.")

        if not isinstance(value, list):
            value = [value]

        # add prefixes to values
        for i, val in enumerate(value):
            value[i] = f"{prefix}{val}"

        template = deepcopy(Feature.QUALIFIER_TEMPLATES[qualifier])
        template.set_value(value)

        self.qualifiers += [template]
