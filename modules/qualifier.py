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
    
    def __init__(self, qualifier, value = [], mandatory = False, qualifier_definition_dir = "modules/qualifiers"):
        self.name = qualifier
        self.mandatory = mandatory
        self._load_definition("%s/%s.json" % (qualifier_definition_dir, qualifier.strip()))
        if hasattr(self, "value_format") and value:
            self.value = self.add_value( list(value) )
        else:
            self.value = list(value)
    
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
        self.legal_dbxref = json.load(open(filename))
    
    def _by_value_format(self, value):
        """
        This is the main function of the qualifier class, it attempts to validate qualifier 
        format by using the value_format tag of the qualifier definition.
        These definitions are written for humans though, not for scripts...
        """
        formatted_value = value
        
        if self.value_format == "none": # no value taken
            if value:
                logging.warn("Qualifier '%s' has value '%s', but %s qualifiers does not take a value" % (self.name, value, self.name))
            return ""
        elif self.value_format.startswith("\"<database:identifier>\""): # Handle dbxref's 
            self._load_legal_dbxref( LEGAL_DBXREF_FILE )
            
            
            new_value=[]
            for val in self.value:
                if val.split(':')[0].lower() in [v.lower() for v in self.legal_dbxref]:
                    new_value.append(val)
                else:
                    logging.warn("Unknown db_xref '%s' - removing." % val)
            formatted_value=new_value
        elif self.value_format == "<identifier>":
            pass
        elif self.value_format == "\"text\"":
            formatted_value = value
            if formatted_value[0] == "\"":
                formatted_value = formatted_value[1:]
            if formatted_value[-1] == "\"":
                formatted_value = formatted_value[:-1]
        else:
            logging.debug("No rule to format '%s'" % self.value_format)
        
        # TODO: implement validators for all value_format types
        
        return formatted_value
    
    def add_value(self, value):
        
        self.value = [self.value] if type(self.value) != type([]) else self.value
        value = [value] if type(value) != type([]) else value
        value = [self._by_value_format(v) for v in value]
        self.value += value if type(value) == type([]) else [value]

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