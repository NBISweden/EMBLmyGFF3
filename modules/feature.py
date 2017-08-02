#!/usr/bin/env python2.7

from __future__ import division

import sys
import json
import logging
from Bio.Seq import Seq
from Bio.Data import CodonTable
from location import EMBLLocation
from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition, AfterPosition
from qualifier import *

def chunk_format(string, chunk_string = None, offset = 0, chunk_size = 3, chunks_per_line = 30, indent = 6):
    offset = offset%chunk_size
    even = True if len(string)%chunk_size == 0 else False
    output = " "*offset
    c_i = 0
    for i, c in enumerate(string):
        output += c
        if (i+offset) and (i+offset) % (chunk_size * chunks_per_line) == (chunk_size * chunks_per_line)-1:
            output += "\n" + " "*indent
            if chunk_string:
                for c in range(chunks_per_line):
                    if c_i >= len(chunk_string):
                        break
                    output += chunk_string[c_i]
                    output += " "*(chunk_size)
                    c_i += 1
                output += "\n" + " "*indent
        elif (i+offset) and (i+offset) % chunk_size == chunk_size-1:
            output += " "
    
    if chunk_string:
        output += "\n" + " "*indent
        for c in range(chunks_per_line):
            if c_i >= len(chunk_string):
                break
            output += chunk_string[c_i]
            output += " "*(chunk_size)
            c_i += 1
    return output

class Feature(object):
    """
    Super-class for the various feature types defined in:
    http://www.insdc.org/files/feature_table.html.
    """
    
    CDS_COUNTER = 0
    OK_COUNTER = 0
    DEFAULT_FEATURE_TRANSLATION_FILE="translation_gff_feature_to_embl_feature.json"
    DEFAULT_QUALIFIER_TRANSLATION_FILE=["translation_gff_attribute_to_embl_qualifier.json", "translation_gff_other_to_embl_qualifier.json"]
    PREVIOUS_ERRORS = []
    
    def __init__(self, feature, seq = None, accessions = [], transl_table = 1, translation_files = [], translate = False, feature_definition_dir = "modules/features", qualifier_definition_dir = "modules/qualifiers", format_data = True, level = 0, reorder_gene_features = True):
        """
        Initializes a Feature, loads json files for feature and 
        qualifiers, and starts parsing the data.
        """
        self.location = feature.location
        self.feature_definition_dir = feature_definition_dir
        self.qualifier_definition_dir = qualifier_definition_dir
        self.qualifiers = {}
        self.qualifier_translation_list = {}
        self.feature_translation_list = {}
        self.singleton_types = ["exon"]
        self.qualifier_prefix = {}
        self.qualifier_suffix = {}
        self.legal_qualifiers = []
        self.sub_features = []
        self.remove = []
        self.translation_files = translation_files
        self._load_qualifier_translations(Feature.DEFAULT_QUALIFIER_TRANSLATION_FILE + translation_files)
        self._load_feature_translations([Feature.DEFAULT_FEATURE_TRANSLATION_FILE])
        self.type = self._from_gff_feature(feature.type)
        self.seq = seq
        self.transl_table = transl_table
        self.translate = translate
        self.level = level
        self.reorder_gene_features = reorder_gene_features

        self._load_definition("%s/%s.json" % (feature_definition_dir, self.type))
        self._load_data(feature, accessions)

        if level == 1:      
            # Parse through subfeatures level2
            featureObj_level2 = None
            for feature_l2 in feature.sub_features:
                featureObj_level2 = Feature(feature_l2, self.seq, accessions, self.transl_table, self.translation_files, self.translate, 
                                                  self.feature_definition_dir, self.qualifier_definition_dir, format_data = True, level=2)
                self.sub_features += [featureObj_level2]

                # Parse through subfeatures level3
                featureObj_level3 = None
                for feature_l3 in feature_l2.sub_features:
                    l3_type = self._from_gff_feature(feature_l3.type)
                    l2_sub_features = [sf.type for sf in featureObj_level2.sub_features]
                    if l3_type in l2_sub_features and l3_type not in self.singleton_types:
                        old_feature = [sf for sf in featureObj_level2.sub_features if sf.type == self._from_gff_feature(feature_l3.type)][0]
                        old_feature.combine(feature_l3)
                    else:
                        featureObj_level3 = Feature(feature_l3, self.seq, accessions, self.transl_table, self.translation_files, self.translate, 
                                                      self.feature_definition_dir, self.qualifier_definition_dir, format_data = False, level=3)
                        featureObj_level2.sub_features += [featureObj_level3]


        if format_data:
            self._format_data(self)

        if self.type == "CDS":
            self.qualifiers['transl_table'].set_value(self.transl_table)
    
    def __repr__(self):
        """
        Formats the feature as EMBL, limited to 80 character lines,
        including sub features.
        """
        output = self._feature_as_EMBL() if self.type not in self.remove else ""
        
        # Sub-features.
        #
        # These need some special formatting - generally features are interleaved, 
        # but for genes they should be printed with sub-features grouped by type
    
        if self.type == "gene" and self.reorder_gene_features:
            
            #print level2
            list_type_l3 = []
            for feature_l2 in self.sub_features:
                output += feature_l2._feature_as_EMBL()
                for feature_l3 in feature_l2.sub_features:
                    if not feature_l3.type in list_type_l3:
                        list_type_l3.append(feature_l3.type)
            #print level3
            for f_type in list_type_l3: 
                for feature_l2 in self.sub_features:
                    for feature_l3 in feature_l2.sub_features:
                        if f_type == feature_l3.type:
                            output += feature_l3._feature_as_EMBL()

        else:
            for sub_feature in self.sub_features:
                output += str(sub_feature)
        
        return output
    
    def _feature_as_EMBL(self):
        """
        Formats the feature as EMBL, limited to 80 character lines.
        """
        
        if self.type == "CDS":
            # with open("feature_%00i.txt" % Feature.CDS_COUNTER, "w") as out:
            #     self.CDS_report(out)
            if self.translate:
                self.qualifiers["translation"].set_value(self.translation())
            Feature.CDS_COUNTER += 1
        
        # Print the feature line (type and location)
        output = ""
        line = "\nFT   %s %s" % ("{:15}".format(self.type), EMBLLocation(self.location))
        if len(line) <= 79:
                output += line
        else: # we have to cut the line (between words)
            output += line[:79]
            line = line[79:]
            while line:
                output += "\nFT                   %s" % line[:59]
                line = line[59:]
        
        # Print qualifiers for the feature
        for qualifier in self.qualifiers.values():
            if qualifier.value:
                output += str(qualifier)
        
        return output
    
    def _format_data(self, feature):
        """
        Reformats the data somewhat to better map to the expected EMBL
        structure
        """
        # according to Jacques, EMBL files shouldn't have an mRNA feature but use the 
        # exon information together with the mRNA features as an mRNA feature
        self._reformat_exons()
        
        # EMBL files are quite picky with complete reading frames, so we check the 
        # features for correct start and stop codons, as well as phase to avoid 
        # errors later.
        self._infer_ORFs(feature)
    
    def _from_gff_feature(self, feature):
        """
        Returns the EMBL feature name from the translation list based
        on the GFF feature name.
        """
        return self.feature_translation_list[feature] if feature in self.feature_translation_list else feature
    
    def _from_gff_qualifier(self, qualifier):
        """
        Returns the EMBL qualifier name from the translation list based
        on the GFF qualifier name.
        """
        return self.qualifier_translation_list[qualifier] if qualifier in self.qualifier_translation_list else qualifier
    
    def _infer_ORFs(self, feature):
        """
        Checks a CDS feature to see if it has a start codon and a stop codon,
        and adjusts the location after that.
        """
        if self.type in ['CDS']:
            seq = self.sequence()
            start_codon = seq[:3]
            stop_codon = seq[-3:]
            
            # load the current codon table
            codon_table = CodonTable.unambiguous_dna_by_id[self.transl_table]
            
            # basic info
            strand = self.location.strand
            
            if start_codon.upper() not in codon_table.start_codons:
                self.location = self._set_before(self.location)
            if stop_codon.upper() not in codon_table.stop_codons:
                self.location = self._set_after(self.location)
            if start_codon.upper() in codon_table.start_codons and stop_codon.upper() in codon_table.stop_codons:
                Feature.OK_COUNTER += 1
            
        for sub_feature in self.sub_features:
            sub_feature._infer_ORFs(feature)
    
    def _load_data(self, feature, accessions):
        """
        Parses a GFF feature and stores the data in the current Feature
        """
        for qualifier, value in feature.qualifiers.iteritems():
            logging.debug("Reading qualifier: %s (%s), translating to %s" % (qualifier, value, self._from_gff_qualifier(qualifier)))
            self.add_qualifier( qualifier, value )
        
        if 'locus_tag' in self.qualifiers:
            self.qualifiers['locus_tag'].set_value( accessions )
    
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
                            self.legal_qualifiers += [key]
                            mandatory = "mandatory" in key
                            self.qualifiers[item] = Qualifier(item, mandatory = mandatory, qualifier_definition_dir=self.qualifier_definition_dir)
                    else:
                        # this is not super important, as it just adds comments and 
                        # description from the documentation to the features. I used 
                        # it to have a bit of debugging information.
                        setattr(self, key, value)
        except IOError as e:
            logging.error(e)
    
    def _load_feature_translations(self, filenames):
        """
        Load translation json files. Files are loaded in order that they are given, 
        thus newer rules can be loaded to replace default rules.
        """
        module_dir = os.path.dirname(os.path.abspath(sys.modules[Feature.__module__].__file__))
        
        for filename in filenames:
            logging.debug("Loading feature translation file: %s/%s" % (module_dir, filename))
            data = json.load( open("%s/%s" % (module_dir, filename)) )
            for gff_feature, info in data.iteritems():
                if info.get("remove", False):
                    self.remove += [gff_feature]
                if "target" in info:
                    self.feature_translation_list[gff_feature] = info["target"]
    
    def _load_qualifier_translations(self, filenames):
        """
        Load translation json files. Files are loaded in order that they are given, 
        thus newer rules can be loaded to replace default rules.
        """
        module_dir = os.path.dirname(os.path.abspath(sys.modules[Feature.__module__].__file__))
        
        for filename in filenames:
            logging.debug("Loading qualifier translation file: %s/%s" % (module_dir, filename))
            data = json.load( open("%s/%s" % (module_dir, filename)) )
            for gff_feature, info in data.iteritems():
                if "target" in info:
                    self.qualifier_translation_list[gff_feature] = info["target"]
                if "prefix" in info:
                    self.qualifier_prefix[gff_feature] = info["prefix"]
                if "suffix" in info:
                    self.qualifier_suffix[gff_feature] = info["suffix"]
    
    def _reformat_exons(self):
        """
        Reformats mRNA features to have the location of its exon sub-features,
        and removes the exon sub-features.
        """
        
        if self.level == 2: # level 2 means e.g: mRNA, tRNA, etc.
            first = True
            for i, sf in enumerate(self.sub_features):
                if sf.type != 'exon':
                    continue
                # replace mRNA location with exon location(s)
                if first:
                    self.location = sf.location
                    first = False
                else:
                    self.location += sf.location
        
        for sf in self.sub_features:
            sf._reformat_exons()
    
    def _set_before(self, location):
        """
        Changes a FeatureLocation to include a "BeforePosition" or 
        "AfterPosition" to indicate that the mRNA does not include 
        start codon.
        """
        if location.strand >= 0: # forward strand
            if len(location.parts) > 1:
                location.parts[0] = FeatureLocation( BeforePosition(location.parts[0].start), location.parts[0].end, strand = location.parts[0].strand )
            else:
                location = FeatureLocation( BeforePosition(location.start), location.end, strand = location.strand)
        else:
            if len(location.parts) > 1:
                location.parts[-1] = FeatureLocation( location.parts[-1].start, AfterPosition(location.parts[-1].end), strand = location.parts[-1].strand )
            else:
                location = FeatureLocation( location.start, AfterPosition(location.end), strand = location.strand)
        return location
    
    def _set_after(self, location):
        """
        Changes a FeatureLocation to include a "BeforePosition" or 
        "AfterPosition" to indicate that the mRNA does not include 
        stop codon.
        """
        if location.strand >= 0: # forward strand
            if len(location.parts) > 1:
                location.parts[-1] = FeatureLocation( location.parts[-1].start, AfterPosition(location.parts[-1].end), strand = location.parts[-1].strand )
            else:
                location = FeatureLocation( location.start, AfterPosition(location.end), strand = location.strand)
        else:
            if len(location.parts) > 1:
                location.parts[0] = FeatureLocation( BeforePosition(location.parts[0].start), location.parts[0].end, strand = location.parts[0].strand )
            else:
                location = FeatureLocation( BeforePosition(location.start), location.end, strand = location.strand)
        return location
    
    def add_qualifier(self, gff_qualifier, value):
        """
        This is where qualifier values are added to the feature.
        """
        qualifier = self._from_gff_qualifier(gff_qualifier)
        logging.debug("Qualifier: %s - %s" % (qualifier, value))
        
        if not qualifier:
            logging.debug("Skipping empty qualifier with value '%s'" % value)
            return
        
        if qualifier not in self.qualifiers:
            try:
                os.stat( "%s/%s.json" % (self.qualifier_definition_dir, qualifier) )
            except Exception as e:
                msg = "Unknown qualifier '%s'" % qualifier
                if msg not in Feature.PREVIOUS_ERRORS:
                    logging.warn(msg)
                    Feature.PREVIOUS_ERRORS += [msg]
            else:
                logging.debug("'%s' is not a legal qualifier for feature type '%s'" % (qualifier, self.type))
                
            return
        
        logging.debug("Adding value '%s' to qualifier '%s'" % (value, qualifier))
        
        if self.qualifier_prefix.get(gff_qualifier, None):
            value = ["%s%s" % (self.qualifier_prefix[gff_qualifier], v) for v in value]
        
        if self.qualifier_suffix.get(gff_qualifier, None):
            value = ["%s%s" % (v, self.qualifier_suffix[gff_qualifier]) for v in value]
        
        self.qualifiers[qualifier].add_value(value)
    
    def combine(self, other):
        """
        Attempt to combine all features from another feature into this one.
        """
        
        # add new location
        self.location += other.location
        
        # combine qualifiers
        for name, qualifier in self.qualifiers.iteritems():
            other_qualifier = other.qualifiers.get(name, None)
            for val in getattr(other_qualifier, "value", []):
                if val not in qualifier.value:
                    self.qualifiers[name].add_value(other_value)
        
        # Sort out phase
        current_phase = int(self.qualifiers.get("phase", [0])[0])
        other_phase = int(other.qualifiers.get("phase", [0])[0])
        
        phase = current_phase if self.location.start < other.location.start else other_phase
        if "codon_start" in self.legal_qualifiers:
            if not "codon_start" in self.qualifiers:
                self.qualifiers["codon_start"] = Qualifier("codon_start", phase, qualifier_definition_dir = self.qualifier_definition_dir)
            else:
                self.qualifiers["codon_start"].set_value(phase)
    
    def CDS_report(self, out = sys.stdout, parts = False, codon_info = True):
        """
        Writes a short report about a CDS to a file, used for debugging.
        """
        seq = "%s" % self.sequence()
        aa = self.translation()
        
        start_codon = seq[:3]
        stop_codon = seq[-3:]
        
        out.write("Name: %s\n" % self.qualifiers.get("gene").value[0])
        out.write("Location: %s\n" % EMBLLocation(self.location))
        
        if parts:
            offset = 0
            aa_offset = 0
            for i, part in enumerate(self.location.parts) if self.location.strand > 0 else enumerate(reversed(self.location.parts)):
                codon_table = CodonTable.unambiguous_dna_by_id[self.transl_table]
                part_seq = SeqFeature(location = part).extract(self.seq)
                aa_len = (len(part)+offset)//3
                part_aa  = aa[aa_offset:aa_offset+aa_len]
                aa_offset += aa_len
            
                out.write("Part %02i %s  " % (i, "(+)" if part.strand > 0 else "(-)"))
                out.write("%s\n" % chunk_format(part_seq, part_aa, offset, 3, 20, 13))
                offset += len(SeqFeature(location = part).extract(self.seq)) % 3
                offset %= 3
        else:
            out.write("Sequence:    %s\n" % chunk_format(seq, indent = 13))
            out.write("Translation: %s\n" % chunk_format(aa, None, 0, 8, 6, 13))
        
        if codon_info:
            codon_table = CodonTable.unambiguous_dna_by_id[self.transl_table]
            out.write("Start codon: %s (%s) \n" % (start_codon, ", ".join(codon_table.start_codons)))
            out.write("Stop codon: %s (%s) \n" % (stop_codon, ", ".join(codon_table.stop_codons)))
    
    def sequence(self):
        """
        Returns the nucleotide sequence of self
        """
        if self.location.strand > 0:
            return SeqFeature(location = self.location).extract(self.seq)
        
        seq = Seq("")
        for part in reversed(self.location.parts):
            seq += SeqFeature(location = part).extract(self.seq)
        
        return seq
    
    def translation(self):
        """
        Returns the amino acid sequence of self
        """
        codon_table = CodonTable.unambiguous_dna_by_id[self.transl_table]
        seq = self.sequence()
        
        output = "%s" % seq.translate(codon_table)
        return output

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
            feature = Feature( gff_feature, args.translation_file, 1, feature_definition_dir = "features", qualifier_definition_dir="qualifiers" )
            print "_"*80
            print feature
            break
    except Exception as e:
        import traceback
        traceback.print_exc(limit=5)
        sys.stderr.write( str(e) )
        
