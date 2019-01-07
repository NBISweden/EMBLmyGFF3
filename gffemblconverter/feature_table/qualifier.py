#!/usr/bin/env python3
"""Qualifier object for EMBL feature tables
"""

from .embl_utilities import embl_line, ensure_quoted

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
        if not isinstance(self.value, list):
            self.value = [self.value]

        output = ""
        for value in self.value:
            if value == []:
                continue

            if isinstance(value, list):
                value = value[0]
            if not value:
                value = "/{self.name}"
            elif isinstance(value, str):
                # quote string values
                value = f"/{self.name}={ensure_quoted(value)}"
            else:
                value = f"/{self.name}={value}"

            output += embl_line("FT", value, False, pad=21)

        return output

    def set_value(self, value):
        """
        Sets the qualifier value, and attempts to verify the value according
        if value_format is set.
        """
        self.value = value if isinstance(value, list) else [value]

    def get_value(self):
        """
        Returns the qualifier value in the default format.
        """
        return self.value
