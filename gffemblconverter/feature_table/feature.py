#!/usr/bin/env python3
"""Feature object for EMBL feature tables
"""

class Feature():
    """Feature object for EMBL feature tables.

    From http://www.insdc.org/files/feature_table.html#3.2:
    "Feature keys indicate
    (1) the biological nature of the annotated feature or
    (2) information about changes to or other versions of the sequence.
    The feature key permits a user to quickly find or retrieve similar features or
    features with related functions. "
    """

    def __init__(self, definition):
        """
        Initializes a Feature from a defintion.
        """
        self.name = definition['feature_key']
        for key, value in definition.items():
            setattr(self, key, value)
        self.value = []

    # def __repr__(self):
    #     return "[Feature: {}]".format(self.name)

    def set_value(self, value):
        """
        Sets the feature value, and attempts to verify the value according
        if value_format is set.
        """
        self.value = value

    def get_value(self):
        """
        Returns the feature value in the default format.
        """
        return self.value
