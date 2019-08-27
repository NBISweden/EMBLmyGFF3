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
    SINGLETON_TYPES = ['exon']

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
                # if the feature supports locus tags, add them to the
                # subfeatures as well
                for sub_feature in self.sub_features:
                    sub_feature.add_locus_tag(value)

    @staticmethod
    def combine(features):
        """
        Attempts to combine a list of Features into a single item that contains
        all the original information.
        """
        combined_feature = features[0]
        codon_starts = []
        for feature in features[1:]:
            # Add location to combined feature
            combined_feature.location += feature.location

            # Combine all qualifiers except 'codon_start'
            for qualifier in feature.qualifiers:
                # Save all the codon_start's for later parsing
                if qualifier == 'codon_start':
                    codon_starts += [(feature.location, qualifier.value)]
                    continue
                # Otherwise - just add the new qualifier
                combined_feature.set_qualifier(qualifier.name, qualifier.value)

        # If we have more than one codon_start, we want the one corresponding to
        # the lowest feature location if we're on the '+' stand, and the highest
        # if we're on the '-' strand.
        if codon_starts:
            reverse = combined_feature.location.strand == 1
            codon_starts.sort(key=lambda x: (x[0].start, x[0].end)[reverse],
                              reverse=reverse)
            combined_feature.set_qualifier('codon_start', codon_starts[0][1])

        return combined_feature

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
                logging.warning("Unknown Feature type: '%s'", template_name)
                raise ValueError("Unknown Feature type")

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
            # Restructure sub-features for EMBL
            template.restructure_subfeatures(feature.sub_features, settings)

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

    def update_locus_tag(self, number):
        """
        Updates the current feature with a locus tag number (replacing the
        placeholder generated by add_locus_tag()). This functions returns True
        if a locus tag was updated, and False if no locus tag was updated.
        """
        try:
            locus_tag_index = self.qualifiers.index('locus_tag')
        except ValueError:
            return False

        locus_tag = self.qualifiers[locus_tag_index]
        locus_tag.value[0] = locus_tag.value[0].format(number=number)

        # update all sub-features with the same locus tag
        for sub_feature in self.sub_features:
            sub_feature.update_locus_tag(number)

        return True

    def restructure_subfeatures(self, sub_features, settings=None):
        """
        Restructures the sub-features of self into a format that is more
        appropriate for EMBL.

        Some features, such as CDS, is handled by a single EMBL field, but by
        separate SeqFeature fields, which makes some restructuring needed. EMBL
        also requires fields to be oredered by location, so we make sure to
        reorder the fields.

        This function is only called for top-level functions, and there is a
        maximum of three levels if the input GFF is formatted as expected.
        """

        # Parse through all sub-features and create features.
        for feature in sub_features:
            embl_feature = Feature.from_seq_feature(feature)
            # loop through all sub-feature types, and combine them if they are
            # not singleton types. Grouping the sub-features by type is needed
            # for gene features, while all other features should be sorted by
            # location only.
            leaves = [Feature.from_seq_feature(f) for f in feature.sub_features]
            # dict keys will keep the original order in python 3.6+, ant I'm not
            # sure that it's important to have them ordered anyway.
            for feature_type in dict.fromkeys([f.name for f in leaves]):
                type_items = [f for f in leaves if f.name in feature_type]

                # If this is a singleton type, add all items as sub-sub-features
                if feature_type in Feature.SINGLETON_TYPES:
                    embl_feature.sub_features += type_items
                # otherwise combine the items and add as a sub-sub-feature
                else:
                    combined = Feature.combine(type_items)
                    # For the special case that this is a CDS, add translation
                    # table information
                    if combined.name == 'CDS':
                        combined.set_qualifier('transl_table',
                                               settings['translation_table'])
                    embl_feature.sub_features += [combined]

            self.sub_features += [embl_feature]

        # If the feature is _not_ a gene, it needs to be sorted by location
        if self.name != 'gene':
            self.sub_features.sort(key=lambda x: x.location.start)

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
        # check if the qualifier already exists, in that case, just add the
        # value to the existing qualifier
        if qualifier in self.qualifiers:
            index = self.qualifiers.index(qualifier)
            self.qualifiers[index].add_value(value)
            return

        # otherwise, check if we need to translate the qualifier name
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

        # Load the template for the chosen qualifier
        if qualifier not in Feature.QUALIFIER_TEMPLATES:
            logging.error("Legal qualifier %s not found in qualifier cache!",
                          qualifier)
            raise ValueError(f"Legal qualifier {qualifier} not found in cache.")

        # Make sure that the value is a list
        if not isinstance(value, list):
            value = [value]

        # Add prefixes to values
        for i, val in enumerate(value):
            value[i] = f"{prefix}{val}"

        # and finally, set the value to the qualifier
        template = deepcopy(Feature.QUALIFIER_TEMPLATES[qualifier])
        template.set_value(value)

        # and save it
        self.qualifiers += [template]
