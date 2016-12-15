#!/usr/bin/env python2.7
"""
Script to parse a "raw" text copy of http://www.insdc.org/files/feature_table.html#7.2 into
a set of json feature identifiers.
"""

import re
import gzip
import json
import logging

class Qualifier( object ):
    """
    Simple feature definition to keep track of things in a reasonable way until we
    write to json.
    """
    
    def __init__(self, key):
        self.name = key
        self.known_identifiers = ["sign", "definition", "value_format", "example", "comment"]
        for identifier in self.known_identifiers:
            setattr(self, identifier, None)
    
    def __repr__(self):
        output  = "Qualifier: %s\n" % self.name
        for identifier in self.known_identifiers:
            value = getattr(self, identifier)
            if value:
                output += "%s: %s\n" % (identifier.replace("_", " ").title(), value)
        return output
    
    def save(self):
        with open('%s.json' % self.name, 'w') as outfile:
            data = {"qualifier":self.name}
            for identifier in self.known_identifiers:
                value = getattr(self, identifier)
                if value:
                    data[identifier] = value
            json.dump(data, outfile, indent=True)

def split_qualifiers(text):
    """
    Looks through 'text' and splits the string on all lines that start with '/' and include a '='.
    The function ignores all '/' characters that are within parentheses, and replaces in-parenthesis
    newlines with spaces. The raw text is hard to parse in cases, but it creates a decent baseline that
    can be manually edited.
    """
    out = []
    temp = ""
    counter = 0
    for i, c in enumerate(text):
        if c == "(":
            counter+=1
        elif c == ")":
            counter-=1
        elif c == "/" and counter == 0: # split on forward slashes, as they're used to mark qualifier names
            if temp:
                out += [temp.replace("\n", " ").strip()]
            temp = c
            continue
        elif c == "\n" and counter > 0: # remove newlines in parentheses
            if temp[-1] != " ":
                temp += " "
            continue
        temp += c
    
    if temp:
        if temp[0] == "/":
            out += [temp.replace("\n", " ").strip()]
        else:
            out[-1] += temp.replace("\n", " ").strip()
    return out

def parse_raw_to_json(infile):
    """
    Parse a raw text file into a number of Qualifiers.
    This function will parse everything as long text strings, and then the 
    Qualifier class will parse qualifiers into lists and things like that.
    """
    
    current = None
    identifier = None
    
    for row in infile:
        if row.startswith("Qualifier"):
            if current:
                current.save()
            name = row[len("Qualifier"):].strip().strip("/=")
            logging.info("Parsing qualifier '%s'" % name)
            current = Qualifier( name )
        elif len(row) <= 1 or re.match("^\s*$", row):
            continue
        elif row.startswith(" "): # Extend current identifier
            try:
                base = getattr(current, identifier)
                if base[-1] != " ":
                    base += " "
                extension = row.strip()
                setattr(current, identifier, base + extension)
            except Exception as e:
                print "EXCEPTION: %s" % e
                print "ID: '%s'" % identifier
                print "ROW: '%s'" % row
                print current
                import sys
                sys.exit(0)
        elif current != None:
            split_point = row.find("  ")
            identifier = row[:split_point].lower().replace(" ", "_").strip()
            if identifier[-1] == "s" and identifier[:-1] in current.known_identifiers:
                identifier = identifier[:-1]
            if identifier not in current.known_identifiers:
                logging.warn("Unknown identifier '%s'" % identifier)
                identifier = None
                logging.info(" ++++ row: %s" % row)
                continue
            value = row[len(identifier):].strip()
            setattr(current, identifier, value)
    if current:
        current.save()


if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("raw", help="raw text file qualifier table dump")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="increase verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="decrease verbosity")
    
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', 
                        level = (5-args.verbose+args.quiet)*10, 
                        datefmt="%H:%M:%S")
    
    in_file = gzip.open(args.raw) if args.raw.endswith(".gz") else open(args.raw)
    
    parse_raw_to_json(in_file)
    
    