#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Qualifier handler for EMBL feature tables
"""

import os
import re
import json
import logging
from utilities import *

LEGAL_DBXREF_FILE="legal_dbxref.json"

class Qualifier( object ):
    
    PREVIOUS_ERRORS = []
    
    def __init__(self, qualifier, value = [], mandatory = False, qualifier_definition_dir = "modules/qualifiers", no_wrap= None):
        """
        Initializes a Qualifier, loads json definition and starts 
        parsing the input data.
        """
        self.name = qualifier
        self.no_wrap = no_wrap
        self.mandatory = mandatory
        self.qualifier_definition_dir = qualifier_definition_dir
        self._load_definition("%s/%s.json" % (qualifier_definition_dir, qualifier.strip()))
        self.value = []
        if value:
            self.add_value( value )
    
    def __repr__(self):
        """
        Prints the current Qualifier in EMBL format
        """
        return self.embl_format()
    
    def _load_definition(self, filename):
        """
        Loads the json definition for a qualifier type and sets its
        own attributes based on the content
        """
        try:
            with open(filename) as data:
                raw = json.load( data )
                for key, value in raw.iteritems():
                    setattr(self, key, value)
        except IOError as e:
            logging.error(e)
    
    def _load_legal_dbxref(self, filename):
        """
        Loads a list of legal external references
        """
        self.legal_dbxref = json.load(open("%s/%s" % (self.qualifier_definition_dir, filename)))
    
    def _by_value_format(self, value):
        """
        This is the main function of the qualifier class, it attempts to validate qualifier 
        format by using the value_format tag of the qualifier definition.
        These definitions are written for humans though, not for scripts...
        TODO: keep implementing these
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

                for val in value:
                    if val.split(':')[0].lower() in [v.lower() for v in self.legal_dbxref]:
                        new_value = val
                    else:
                        msg = "Unknown db_xref '%s' - removing." % (val.split(':')[0])
                        if msg not in Qualifier.PREVIOUS_ERRORS:
                            logging.info(msg)
                            Qualifier.PREVIOUS_ERRORS += [msg]
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

    def embl_format(self, no_wrap=None):
        """
        Format lines as EMBL format, limiting lines to 80 characters
        """
        output = ""
        value = self.value
        if type(value) != type([]):
            value = [value]
        for val in value:
            string=""
            if val == []:
                continue
            if type(val) == type([]) and len(val):
                val = val[0]
            if not val: # if empty => no value
                string = "/{}".format(self.name)
            elif type(val) in [type(""), type(u"")]: #if a string => we add quote to the value
                string = "/{}=\"{}\"".format(self.name, val)
            else:# A value but is an integer => we do not add quote
                 string = "/{}={}".format(self.name, val)
            #if getattr(self, "value_type", None) == "none":
            #    string = "/%s" % (self.name)
            #    continue
            
            output += multiline("FT", string, wrap=59, no_wrap = no_wrap)
            
        return output
         
    def add_value(self, value):
        """
        Adds a value to the current set of values for this qualifier.
        Qualifier wihtout value is possible. e.g environmental_sample
        """
        self.value = [self.value] if type(self.value) != type([]) else self.value
        
        new_value = self._by_value_format(value)
        logging.debug("%s - Changing value: '%s' to '%s'" % (self.name, value, new_value))
        self.value.append(new_value)
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
