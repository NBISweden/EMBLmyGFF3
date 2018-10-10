#!/usr/bin/env python3
"""Qualifier object for EMBL feature tables
"""

class Qualifier():
    """Qualifier object for EMBL feature tables.

    From http://www.insdc.org/files/feature_table.html#3.3:
    "Qualifiers provide a general mechanism for supplying information about
    features in addition to that conveyed by the key and location."
    """

    def __init__(self, definition):
        """
        Initializes a Qualifier from a defintion.
        """
        self.name = definition['qualifier']
        for key, value in definition.items():
            setattr(self, key, value)
        self.value = []

    def __repr__(self):
        return "[Qualifier: {}]".format(self.name)

    def set_value(self, value):
        """
        Sets the qualifier value, and attempts to verify the value according
        if value_format is set.
        """
        self.value = value

    def get_value(self):
        """
        Returns the qualifier value in the default format.
        """
        return self.value
