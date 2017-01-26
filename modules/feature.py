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

#
# required for next project:
#
#   NcRNA -> ncRNA
#   Non_canonical_five_prime_splice_site
#   Non_canonical_three_prime_splice_site
#   Rfam
#   Transcript -?-> does 'prim_transcript' work?


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
    DEFAULT_TRANSLATION_FILE="translation_gff_to_embl.json"
    
    def __init__(self, feature, seq = None, accessions = [], transl_table = 1, translation_files = [], feature_definition_dir = "modules/features", qualifier_definition_dir = "modules/qualifiers", format_data = True):
        self.type = feature.type
        self.seq = seq
        self.transl_table = transl_table
        
        self.location = feature.location
        self.feature_definition_dir = feature_definition_dir
        self.qualifier_definition_dir = qualifier_definition_dir
        self.qualifiers = {}
        self.translation_list = {}
        self.sub_features = []
        self.translation_files = translation_files
        self._load_translations([self.DEFAULT_TRANSLATION_FILE] + translation_files)
        self._load_definition("%s/%s.json" % (feature_definition_dir, feature.type))
        self._load_data(feature, accessions)
        if format_data:
            self._format_data(feature)
    
    def __repr__(self):
        """
        Formats the feature as EMBL, limited to 80 character lines.
        """
        
        if self.type == "CDS":
            with open("feature_%00i.txt" % Feature.CDS_COUNTER, "w") as out:
                self.CDS_report(out)
            Feature.CDS_COUNTER += 1
        
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
            
        for qualifier in self.qualifiers.values():
            if qualifier.value:
                output += str(qualifier)
        for sub_feature in self.sub_features:
            output += str(sub_feature)
        return output
    
    def _from_gff_qualifier(self, qualifier):
        return self.translation_list[qualifier] if qualifier in self.translation_list else qualifier
    
    def _load_data(self, feature, accessions):
        for qualifier, value in feature.qualifiers.iteritems():
            logging.debug("Reading qualifier: %s (%s), translating to %s" % (qualifier, value, self._from_gff_qualifier(qualifier)))
            self.add_qualifier( self._from_gff_qualifier(qualifier), value )
        
        # create locus tag from feature ID and accessions
        output_accession="|".join(accessions if type(accessions) == type([]) else [accessions])
        
        if 'locus_tag' in self.qualifiers:
            locus_tag = "%s_%s" % (output_accession, "_".join(feature.qualifiers['ID']))
            self.qualifiers['locus_tag'].set_value( locus_tag )
        
        # Parse through subfeatures
        sub_feature_types = []
        for sub_feature in feature.sub_features:
            if sub_feature.type in [sf.type for sf in self.sub_features]:
                old_feature = [sf for sf in self.sub_features if sf.type == sub_feature.type][0]
                old_feature.combine(sub_feature)
            else:
                self.sub_features += [Feature(sub_feature, self.seq, accessions, self.transl_table, self.translation_files, 
                                              self.feature_definition_dir, self.qualifier_definition_dir, format_data = False)]
    
    def _format_data(self, feature):
        # according to Jacques, EMBL files shouldn't have an mRNA feature but use the 
        # exon information together with the mRNA features as an mRNA feature
        self._reformat_exons()
        
        # EMBL files are quite picky with complete reading frames, so we check the 
        # features for correct start and stop codons, as well as phase to avoid 
        # errors later.
        self._infer_ORFs(feature)
    
    def _reformat_exons(self):
        if self.type == "mRNA":
            for i, sf in enumerate(self.sub_features):
                if sf.type != 'exon':
                    continue
                # replace mRNA location with exon location(s)
                self.location = sf.location
                del self.sub_features[i]
                break
        
        for sf in self.sub_features:
            sf._reformat_exons()
    
    def _infer_ORFs(self, feature):
        if self.type in ['CDS']:
            # create a sequence feature using the location of the current EMBLFeature
            # seq = SeqFeature(location = self.location).extract(self.seq)
            # start_codon = seq[:3].upper()
            
            offset = len(self.location)%3 if len(self.location) %3 != 0 else 0
            if  self.location.strand < 0:
              start_codon=self.seq[ self.location.end - 3 : self.location.end].reverse_complement()
              stop_codon=self.seq[ self.location.start - 3 + offset: self.location.start + offset].reverse_complement()
            else:
              start_codon=self.seq[ self.location.start : self.location.start + 3]
              stop_codon=self.seq[ self.location.end - offset: self.location.end + 3 - offset]
            
            #logging.error("%s - %s %s" % (seq_start_codon, start_codon, "<----" if seq_start_codon == start_codon else ""))
            
            # load the current codon table
            codon_table = CodonTable.unambiguous_dna_by_id[self.transl_table]
            
            if start_codon not in codon_table.start_codons:
                strand = self.location.parts[0].strand
                start = BeforePosition(self.location.parts[0].start)
                end = self.location.parts[0].end
                self.location.parts[0] = FeatureLocation( start, end, strand = strand )
            if stop_codon not in codon_table.stop_codons:
                strand = self.location.parts[-1].strand
                start = self.location.parts[-1].start
                end = AfterPosition(self.location.parts[-1].end)
                self.location.parts[-1] = FeatureLocation( start, end, strand = strand )
            
        for sub_feature in self.sub_features:
            sub_feature._infer_ORFs(feature)
    
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
        module_dir = os.path.dirname(os.path.abspath(sys.modules[Feature.__module__].__file__))
        
        for filename in filenames:
            logging.debug("Loading translation file: %s/%s" % (module_dir, filename))
            data = json.load( open("%s/%s" % (module_dir, filename)) )
            for gff_feature, info in data.iteritems():
                self.translation_list[gff_feature] = info["target"]
    
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
                logging.debug("'%s' is not a legal qualifier for feature type '%s'" % (qualifier, self.type))
                
            return
        
        logging.debug("Adding value '%s' to qualifier '%s'" % (value, qualifier))
        
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
        if not "codon_start" in self.qualifiers:
            self.qualifiers["codon_start"] = Qualifier("codon_start", phase, qualifier_definition_dir = self.qualifier_definition_dir)
        else:
            self.qualifiers["codon_start"].set_value(phase)
    
    def CDS_report(self, out = sys.stdout, parts = False, codon_info = True):
        seq = "%s" % self.sequence()
        aa = self.translation()
        
        offset = len(self.location)%3
        if  self.location.strand < 0:
            start_codon=self.seq[ self.location.end - 3 : self.location.end].reverse_complement()
            stop_codon=self.seq[ self.location.start - 3 + offset : self.location.start + offset].reverse_complement()
        else:
            start_codon=self.seq[ self.location.start : self.location.start + 3]
            stop_codon=self.seq[ self.location.end - offset: self.location.end + 3 - offset]
        
        out.write("Name: %s\n" % self.type)
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
            out.write("Start codon: %s (%s)\n" % (start_codon, ", ".join(codon_table.start_codons)))
            out.write("Stop codon: %s (%s)\n" % (stop_codon, ", ".join(codon_table.stop_codons)))
    
    def sequence(self):
        seq = Seq("", self.seq.alphabet)
        if self.location.strand > 0:
            for part in self.location.parts:
                seq += SeqFeature(location = part).extract(self.seq)
        else:
            for part in reversed(self.location.parts):
                seq += SeqFeature(location = part).extract(self.seq)
        
        return seq
    
    def translation(self):
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
        