#!/usr/bin/env python2.7
"""
Qualifier handler for EMBL feature tables
"""

import os
import re
import json
import logging

LEGAL_DBXREF_FILE="legal_dbxref.json"

class Qualifier( object ):
    
    SINGLE_VALUE_QUALIFIERS = ["gene"]
    
    def __init__(self, qualifier, value = [], mandatory = False, qualifier_definition_dir = "modules/qualifiers"):
        self.name = qualifier
        self.mandatory = mandatory
        self.qualifier_definition_dir = qualifier_definition_dir
        self._load_definition("%s/%s.json" % (qualifier_definition_dir, qualifier.strip()))
        self.value = []
        if value:
            self.add_value( value )
    
    def __repr__(self):
        return self._embl_format()
    
    def _embl_format(self):
        """
        Format lines as EMBL format, limiting lines to 80 characters
        """
        output = ""
        value = self.value
        if type(value) != type([]):
            value = [value]
        for val in value:
            if val == []:
                continue
            if type(val) == type([]) and len(val):
                val = val[0]
            if type(val) == type(""):
                val = "\"%s\"" % val
            if getattr(self, "value_type", None) == "none":
                output += "\nFT                   /%s\n"
                continue
            line = "\nFT                   /%s=%s" % (self.name, val)
            if len(line) <= 79:
                output += line
            else:
                output += line[:79]
                line = line[79:]
                while line:
                    output += "\nFT                   %s" % line[:59]
                    line = line[59:]
            # break if only one value is legal
            if self.name in self.SINGLE_VALUE_QUALIFIERS:
                break
        return output
    
    def _load_definition(self, filename):
        try:
            with open(filename) as data:
                raw = json.load( data )
                for key, value in raw.iteritems():
                    setattr(self, key, value)
        except IOError as e:
            logging.error(e)
    
    def _load_legal_dbxref(self, filename):
        self.legal_dbxref = json.load(open("%s/%s" % (self.qualifier_definition_dir, filename)))
    
    def _by_value_format(self, value):
        """
        This is the main function of the qualifier class, it attempts to validate qualifier 
        format by using the value_format tag of the qualifier definition.
        These definitions are written for humans though, not for scripts...
        """
        formatted_value = value
        
        if not hasattr(self, "value_format"):
            return value
        
        try:
            if self.value_format == "none": # no value taken
                if value:
                    logging.warn("Qualifier '%s' has value '%s', but %s qualifiers does not take a value" % (self.name, value, self.name))
                return ""
            elif self.value_format.startswith("\"<database:identifier>\""): # Handle dbxref's 
                self._load_legal_dbxref( LEGAL_DBXREF_FILE )
                value = value if type(value) == type([]) else [value]
                new_value=[]
                for val in value:
                    if val.split(':')[0].lower() in [v.lower() for v in self.legal_dbxref]:
                        new_value.append(val)
                    else:
                        logging.info("Unknown db_xref '%s' - removing." % val)
                formatted_value=new_value
            elif self.value_format == "<identifier>":
                pass
            elif self.value_format == "\"text\"":
                formatted_value = value
                if formatted_value[0] == "\"":
                    formatted_value = formatted_value[1:]
                if formatted_value[-1] == "\"":
                    formatted_value = formatted_value[:-1]
            elif self.value_format == "1 or 2 or 3":
                formatted_value = int(value) + 1
                if formatted_value not in [1,2,3]:
                    logging.error("Value format '1 or 2 or 3' has value: %i" % formatted_value)
            else:
                logging.debug("No rule to format '%s'" % self.value_format)
        
            # TODO: implement validators for all value_format types
        except Exception as e:
            import traceback
            logging.error(e)
            logging.info("---------+----------------------------------------")
            logging.info("QUALIFIER|  Name: %s" % self.name)
            logging.info("---------+ Value: %s" % value)
            logging.info(" Formatted value: %s" % formatted_value)
            logging.info("       Mandatory: %s" % self.mandatory)
            logging.info("    Value format: %s" % self.value_format)
            logging.info("-"*50 + "\n")
            #logging.info("      dir(self): %s" % dir(self))
            traceback.print_exc(limit=5)
            import sys
            sys.exit(1)
        
        return formatted_value
    
    def add_value(self, value):
        if not value:
            logging.debug("%s - Not adding value (current value: %s)" % (self.name, self.value))
            return
        self.value = [self.value] if type(self.value) != type([]) else self.value
        value = [value] if type(value) != type([]) else value
        
        new_values = [self._by_value_format(v) for v in value]
        logging.debug("%s - Changing value: '%s' to '%s'" % (self.name, value, new_values))
        self.value += new_values
        logging.debug("%s - Current value: '%s'" % (self.name, self.value))
    
    def set_value(self, value):
        value = value if type(value) == type([]) else [value]
        self.value = [self._by_value_format(v) for v in value]

if __name__ == '__main__':
    
    import argparse
    from BCBio import GFF
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("test_gff", help="gff file to load a test qualifier from")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="increase verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="decrease verbosity")
    
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', 
                        level = (5-args.verbose+args.quiet)*10, 
                        datefmt="%H:%M:%S")
    
    for record in GFF.parse(open(args.test_gff)):
        break
    
    for feature in record.features:
        for qualifier, value in feature.qualifiers.iteritems():
            print "%s: %s" % (qualifier, value)
            print Qualifier(qualifier, value, "qualifiers")
        break