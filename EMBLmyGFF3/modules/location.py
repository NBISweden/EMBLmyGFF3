#!/usr/bin/env python2.7
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
        complement = True
        if len(self.location.parts) == len([l for l in self.location.parts if l.strand < 0]):
            output += "complement("
            suffix += ")"
            complement = False
        if (len(self.location.parts) > 1):
            output += "join("
            suffix += ")"
        output += ",".join(self._format_parts(self.location.parts, complement=complement))
            
        return output + suffix
    
    def _format_parts(self, parts, complement = True):
        output = []
        for part in parts:
            if part.strand > 0 or complement == False:
                output += ["%s..%s" % (type(part.start)(part.start+1), type(part.end)(part.end+0))]
            else:
                output += ["complement(%s..%s)" % (type(part.start)(part.start+1), type(part.end)(part.end+0))]
        return output
