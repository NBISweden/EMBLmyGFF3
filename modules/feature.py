#!/usr/bin/env python2.7

import sys
import json
from Bio.Seq import Seq
from Bio.Data import CodonTable
from location import Location
from Bio.SeqFeature import FeatureLocation
from qualifier import *

#
# required for next project:
#
#   NcRNA -> ncRNA
#   Non_canonical_five_prime_splice_site
#   Non_canonical_three_prime_splice_site
#   Rfam
#   Transcript -?-> does 'prim_transcript' work?

class Feature(object):
    """
    Super-class for the various feature types defined in:
    http://www.insdc.org/files/feature_table.html.
    """
    
    DEFAULT_TRANSLATION_FILE="translation_gff_to_embl.json"
    
    def __init__(self, feature, seq = None, translation_files = [], feature_definition_dir = "modules/features", qualifier_definition_dir = "modules/qualifiers"):
        logging.debug("Loading feature: %s" % feature.type)
        self.type = feature.type
        self.seq = seq
        self.location = Location(feature, feature.location)
        self.feature_definition_dir = feature_definition_dir
        self.qualifier_definition_dir = qualifier_definition_dir
        self.qualifiers = {}
        self.translation_table = {}
        self.sub_features = []
        self.translation_files = translation_files
        self._load_translations([self.DEFAULT_TRANSLATION_FILE] + translation_files)
        self._load_definition("%s/%s.json" % (feature_definition_dir, feature.type))
        self._load_data(feature)
    
    def __repr__(self):
        """
        Formats the feature as EMBL, limited to 80 character lines.
        """
        output = ""
        line = "\nFT   %s %s" % ("{:15}".format(self.type), self.location)
        if len(line) <= 79:
                output += line
        else: # we have to cut the line (between words)
            output += line[:79]
            line = line[79:]
            while line:
                output += "\nFT                   %s" % line[:59]
                line = line[59:]
            
        for qualifier in self.qualifiers.values():
            if qualifier.value:
                output += str(qualifier)
        for sub_feature in self.sub_features:
            output += str(sub_feature)
        return output
    
    def _from_gff_qualifier(self, qualifier):
        return self.translation_table[qualifier] if qualifier in self.translation_table else qualifier
    
    def _load_data(self, feature):
        for qualifier, value in feature.qualifiers.iteritems():
            logging.debug("Reading qualifier: %s (%s), translating to %s" % (qualifier, value, self._from_gff_qualifier(qualifier)))
            self.add_qualifier( self._from_gff_qualifier(qualifier), value )
        
        for sub_feature in feature.sub_features:
            self.sub_features += [Feature(sub_feature, self.seq, self.translation_files, self.feature_definition_dir, self.qualifier_definition_dir)]
    
    def _load_definition(self, filename):
        """
        Loads a Feature definition json file.
        """
        try:
            with open(filename) as data:
                raw = json.load( data )
                for key, value in raw.iteritems():
                    if "qualifier" in key:
                        for item, definition in value.iteritems():
                            mandatory = "mandatory" in key
                            self.qualifiers[item] = Qualifier(item, mandatory = mandatory, qualifier_definition_dir=self.qualifier_definition_dir)
                    else:
                        # this is not super important, as it just adds comments and 
                        # description from the documentation to the features. I used 
                        # it to have a bit of debugging information.
                        setattr(self, key, value)
        except IOError as e:
            print e
    
    def _load_translations(self, filenames):
        """
        Load translation json files. Files are loaded in order that they are given, 
        thus newer rules can be loaded to replace default rules.
        """
        for filename in filenames:
            logging.debug("Loading translation file: %s/%s" % (self.feature_definition_dir, filename))
            data = json.load( open("%s/%s" % (self.feature_definition_dir, filename)) )
            for gff_feature, info in data.iteritems():
                self.translation_table[gff_feature] = info["target"]
    
    def add_qualifier(self, qualifier, value):
        """
        This is where qualifier values are added to the feature.
        """
        logging.debug("Qualifier: %s - %s" % (qualifier, value))
        
        if not qualifier:
            logging.debug("Skipping empty qualifier with value '%s'" % value)
            return
        
        if qualifier not in self.qualifiers:
            try:
                os.stat( "%s/%s.json" % (self.qualifier_definition_dir, qualifier) )
            except Exception as e:
                logging.warn("Unknown qualifier '%s'" % qualifier)
            else:
                logging.warn("'%s' is not a legal qualifier for feature type '%s'" % (qualifier, self.type))
                
            return
        
        logging.debug("Adding value '%s' to qualifier '%s'" % (value, qualifier))
        
        self.qualifiers[qualifier].add_value(value)

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("gff_file", help="Annotation file in GFF3 format")
    parser.add_argument("--translation_file", default=[], nargs="+", help="GFF to EMBL translation file(s) to load.")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="increase verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="decrease verbosity")
    
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', 
                        level = (5-args.verbose+args.quiet)*10, 
                        datefmt="%H:%M:%S")
    
    try:
        from BCBio import GFF
        for record in GFF.parse( args.gff_file ):
            break
        
        for gff_feature in record.features:
            print gff_feature
            print "_"*80
            feature = Feature( gff_feature, args.translation_file, feature_definition_dir = "features", qualifier_definition_dir="qualifiers" )
            print "_"*80
            print feature
            break
    except Exception as e:
        import traceback
        traceback.print_exc(limit=5)
        sys.stderr.write( str(e) )
        