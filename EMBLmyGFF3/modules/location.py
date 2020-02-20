#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Location examples (http://www.insdc.org/files/feature_table.html)

This used to be a bigger class, but it's been re-written to just be 
a formatting wrapper around Bio.SeqFeature.FeatureLocation.
"""


import logging

class EMBLLocation(object):
    
    def __init__(self, location):
        self.location = location
    
    def __repr__(self):
        
        output = ""
        suffix = ""

        # Store all strand from the location parts to check potential inconsistency 
        strand=[]
        for l in self.location.parts:
            if l.strand != None:
                if l.strand not in strand:
                    strand.append(l.strand)

        if len(strand) == 0:
            logging.debug("No strand stored among the location_parts %s" % self.location.parts)       
        elif len(strand) > 1:
            logging.error("Different strand stored in location_parts (+ strand will be used as default): %s" % self.location.parts)
        elif strand == [1]:
            logging.debug("+ strand")
        else:
            logging.debug("- strand")
            output += "complement("
            suffix += ")"

        # If more than one part let's join the differnt parts together
        if (len(self.location.parts) > 1):
            output += "join("
            suffix += ")"

        output += ",".join(self._format_parts(self.location.parts))
            
        return output + suffix
    
    def _format_parts(self, parts):
        output = []
        for part in parts:
            output += ["%s..%s" % (type(part.start)(part.start+1), type(part.end)(part.end+0))]
        return output
