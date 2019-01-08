#!/usr/bin/env python3
"""Qualifier object for EMBL feature tables
"""

import logging

from .embl_utilities import embl_line, quoted

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
        self.value_format = None
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
                value = f"/{self.name}={quoted(value)}"
            else:
                value = f"/{self.name}={value}"

            output += embl_line("FT", value, False, pad=21)

        return output

    def validate_value(self, value):
        """
        This function attempts to validate qualifier format by using the
        value_format tag of the qualifier definition. Currently, only the
        simpler cases are handled, but given time, more validations will be
        added.
        """
        if isinstance(value, list):
            return [self.validate_value(v) for v in value]

        if self.value_format is None:
            return value

        formatted_value = value
        if self.value_format == "none": # no value taken
            if value:
                logging.warning(("Qualifier '%s' has value '%s', but %s "
                                 "qualifier does not take a value", self.name,
                                 value, self.name))
            return ""
        if self.value_format == "<identifier>":
            pass
        elif self.value_format == "\"text\"":
            formatted_value = quoted(value)
        elif self.value_format == "1 or 2 or 3":
            formatted_value = int(value) + 1
            if formatted_value not in [1, 2, 3]:
                logging.error("Value format '1 or 2 or 3' has value: %i",
                              formatted_value)

        return formatted_value

    def add_value(self, value):
        """
        Runs the qualifier value validator, and adds a value to the current list
        of qualifier values
        """
        value = self.validate_value(value)
        if not isinstance(self.value, list):
            self.set_value(self.value)
        self.value += [value]

    def set_value(self, value):
        """
        Runs the qualifier value validator and sets the qualifier value.
        """
        value = self.validate_value(value)
        self.value = value if isinstance(value, list) else [value]
