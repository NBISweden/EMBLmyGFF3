#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .qualifier import *
from .location import EMBLLocation
from .utilities import *

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet.IUPAC import *
from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition, AfterPosition

import os
import sys
import json
import logging
from operator import attrgetter

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

    Often the stucture is complex. We build a the Feature object based on a bcbio Record. A record is in fact a group of features linked together by relationship. e.g: A gene has mRNA as child which has exon as child. The exon has mRNA as parent and mRNA has gene has parent.
    Would be nice to make create a Record objcet made from Feature object. We are mixing concepts here.
    Some information are udeful only a the top feature. To make the child feature lighter no need to attach tehm all the information.
    """

    CDS_COUNTER = 0
    OK_COUNTER = 0
    DEFAULT_FEATURE_TRANSLATION_FILE=["translation_gff_feature_to_embl_feature.json"]
    DEFAULT_QUALIFIER_TRANSLATION_FILE=["translation_gff_attribute_to_embl_qualifier.json", "translation_gff_other_to_embl_qualifier.json"]
    PREVIOUS_ERRORS = {}
    SCRIPT_DIR = os.path.dirname(__file__)
    DEFAULT_FEATURE_DIR = os.path.join(SCRIPT_DIR, "modules/features")
    DEFAULT_QUALIFIER_DIR = os.path.join(SCRIPT_DIR, "modules/qualifiers/")

    def __init__(self, feature = None, seq = None, accessions = [], transl_table = 1, translation_files = [], translate = False,
                feature_definition_dir = DEFAULT_FEATURE_DIR, qualifier_definition_dir = DEFAULT_QUALIFIER_DIR, format_data = True,
                level = 0, reorder_gene_features = True, skip_feature = False, force_unknown_features = False,
                force_uncomplete_features = False, uncompressed_log = None, no_wrap_qualifier = False):
        """
        Initializes a Feature, loads json files for feature and
        qualifiers, and starts parsing the data.
        """
        self.feature = feature
        self.legal_qualifiers = []
        self.level = level
        self.location = feature.location
        self.qualifier_definition_dir = qualifier_definition_dir
        self.qualifiers = {}
        self.qualifier_translation_list = {}
        self.qualifier_prefix = {}
        self.qualifier_suffix = {}
        self.feature_definition_dir = feature_definition_dir
        self.feature_translation_list = {}
        self.force_uncomplete_features = force_uncomplete_features
        self.force_unknown_features = force_unknown_features
        self.no_wrap_qualifier = no_wrap_qualifier
        self.reorder_gene_features = reorder_gene_features
        self.remove = []
        self.seq = seq
        self.combine_types = ["CDS","3'UTR","5'UTR"]
        self.skip_feature = skip_feature
        self.sub_features = []
        self.translate = translate
        self.translation_files = translation_files
        self.transl_table = transl_table
        self.uncompressed_log = uncompressed_log

        self._load_qualifier_translations(Feature.DEFAULT_QUALIFIER_TRANSLATION_FILE + translation_files)
        self._load_feature_translations(Feature.DEFAULT_FEATURE_TRANSLATION_FILE)
        self.type = self._from_gff_feature(feature.type)
        self._load_definition("%s/%s.json" % (feature_definition_dir, self.type))
        self._load_data(feature, accessions)
        self._check_qualifier(feature)

        if level == 1 :
            # Parse through subfeatures level2
            featureObj_level2 = None
            for feature_l2 in feature.sub_features:
                featureObj_level2 = Feature(feature_l2, self.seq, accessions, self.transl_table, self.translation_files, self.translate,
                                                      self.feature_definition_dir, self.qualifier_definition_dir, format_data = True, level=2)
                self.sub_features += [featureObj_level2]


                # Sort the sub-features in case they were not ordered into the gff3 file (issue 1)
                feature_l2.sub_features.sort(key=lambda x: x.location.start)

                # Parse through subfeatures level3
                featureObj_level3 = None
                for feature_l3 in feature_l2.sub_features:
                    l3_type = self._from_gff_feature(feature_l3.type)
                    l2_sub_features = [sf.type for sf in featureObj_level2.sub_features]
                    if l3_type in l2_sub_features and l3_type in self.combine_types:
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
        output=str("")

        if self.skip_feature is False or self.force_unknown_features or self.force_uncomplete_features:
            output = self._feature_as_EMBL(self.no_wrap_qualifier) if self.type not in self.remove else ""

        # Sub-features.
        #
        # These need some special formatting - generally features are interleaved,
        # but for genes they should be printed with sub-features grouped by type

        if self.type == "gene" and self.reorder_gene_features:

            #print level2
            list_type_l3 = []
            for feature_l2 in self.sub_features:

                for feature_l3 in feature_l2.sub_features:
                    if not feature_l3.type in list_type_l3:
                        list_type_l3.append(feature_l3.type)

                if feature_l2.skip_feature is False or self.force_unknown_features or self.force_uncomplete_features:
                    output += feature_l2._feature_as_EMBL(self.no_wrap_qualifier) if feature_l2.type not in feature_l2.remove else ""
                #else:
                    #check if CDS exist in subfeature. It could be helpful to create a mRNA feature instead to skip stupidly the L2 feature ! But it's not the philosophy of the tool. It should be done by using the json mapping file.
                #    if "CDS" in list_type_l3:
                #        feature_l2.type = "mRNA"
                #        output += feature_l2._feature_as_EMBL()

            #print level3
            for f_type in list_type_l3:
                for feature_l2 in self.sub_features:
                    for feature_l3 in feature_l2.sub_features:
                        if f_type == feature_l3.type:
                            if feature_l3.skip_feature is False or self.force_unknown_features or self.force_uncomplete_features:
                                output += feature_l3._feature_as_EMBL(self.no_wrap_qualifier) if feature_l3.type not in feature_l3.remove else ""

        else:
            for sub_feature in self.sub_features:
                if sub_feature.skip_feature is False or self.force_unknown_features or self.force_uncomplete_features:
                    output += str(sub_feature) if sub_feature.type not in sub_feature.remove else ""

        return output

    def _feature_as_EMBL(self, no_wrap_qualifier):
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
        string = "%s" % (EMBLLocation(self.location))
        output = multiline("FT", string, featureType=self.type, wrap=59, split_char=",")

        # Print qualifiers for the feature
        for qualifier in sorted(self.qualifiers): # sort by qualifier name

            # continue if qualifier has a value
            if self.qualifiers[qualifier].value:
                # sort by value
                self.qualifiers[qualifier].value = sorted(self.qualifiers[qualifier].value)
                output += self.qualifiers[qualifier].embl_format(no_wrap_qualifier)

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
            codon_table = CodonTable.ambiguous_dna_by_id[self.transl_table]

            # basic info
            strand = self.location.strand
            # raise an error if no strand for the CDS. Strand is not mandatory (can be a dot) except for CDS where it has an
            # impact on the translation, and to check where is start and stop codon...
            if strand == None:
                ID=''
                for qualifier in self.feature.qualifiers:
                    if 'id' == qualifier.lower():
                        ID =  "%s" % " ".join(self.feature.qualifiers[qualifier])
                        break
                logging.error('CDS %s does not have any strand! Please check your gff file.'  %  ID)
                sys.exit()

            if start_codon.upper() not in codon_table.start_codons:
                self.location = self._set_before(self.location)
            if stop_codon.upper() not in codon_table.stop_codons:
                self.location = self._set_after(self.location)
            if start_codon.upper() in codon_table.start_codons and stop_codon.upper() in codon_table.stop_codons:
                Feature.OK_COUNTER += 1

        for sub_feature in self.sub_features:
            sub_feature._infer_ORFs(feature)

    def _check_qualifier(self, feature):
        for qualifier, value in self.qualifiers.items():

            # Check presence of mandatory qualifier
            if self.qualifiers[qualifier].mandatory:# Check if the qualifier is mandatory
                if not self.qualifiers[qualifier].value: # No value for this mandatory qualifier

                    msg = "The qualifier >%s< is mandatory for the feature >%s<. We will not report the feature." % (qualifier, self.type)
                    self.handle_message('warning', msg, msg, None)

                    self.skip_feature = True

    def _load_data(self, feature, accessions):
        """
        Parses a GFF feature and stores the data in the current Feature
        """
        for qualifier, value in feature.qualifiers.items():
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
                for key, value in raw.items():
                    #logging.error("key:%s value:%s",key,value)
                    if "qualifier" in key:
                        for item, definition in value.items():
                            #logging.error("item:%s definition:%s",item,definition)
                            self.legal_qualifiers += [item]
                            mandatory = "mandatory" in key
                            self.qualifiers[item] = Qualifier(item, mandatory = mandatory, qualifier_definition_dir=self.qualifier_definition_dir)
                    else:
                        # this is not super important, as it just adds comments and
                        # description from the documentation to the features. I used
                        # it to have a bit of debugging information.
                        setattr(self, key, value)
        except IOError as e:
            logging.debug("%s" % e)
            msg = ">>%s<< is not a valid EMBL feature type. You can ignore this message if you don't need the feature.\nOtherwise tell me which EMBL feature it corresponds to by adding the information within the json mapping file." % (self.type)
            self.handle_message("error", msg, msg, 1)

            self.skip_feature=True

    def _load_feature_translations(self, filenames):
        """
        Load translation json files. Files are loaded in order that they are given,
        thus newer rules can be loaded to replace default rules.
        """
        module_dir = os.path.dirname(os.path.abspath(sys.modules[Feature.__module__].__file__))
        local_dir = os.getcwd()

        for filename in filenames:
            try:
                data = json.load( open("%s/%s" % (local_dir, filename)) )
            except IOError:
                data = json.load( open("%s/%s" % (module_dir, filename)) )
            for gff_feature, info in data.items():
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
        local_dir = os.getcwd()

        for filename in filenames:
            try:
                data = json.load( open("%s/%s" % (local_dir, filename)) )
            except IOError:
                data = json.load( open("%s/%s" % (module_dir, filename)) )
            for gff_feature, info in data.items():
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
                msg = "Unknown qualifier '%s' - skipped" % qualifier
                self.handle_message('warning', msg, msg, 1)
            else:
                logging.debug("'%s' is not a legal qualifier for feature type '%s'" % (qualifier, self.type))

            return

        # Check if qualifier follow rules (i.e regex)
        if self.qualifiers[qualifier].value_format.startswith("regex:"):
            error_regex=False
            regex= self.qualifiers[qualifier].value_format[6:]
            pattern = re.compile(regex)

            newListValue=[]
            if type(value) == type([]):
                for val in value:
                    if pattern.search(val):
                        newListValue+=value

                if not newListValue:
                    error_regex = True
                else:
                    value = newListValue

            elif not pattern.search(value):
                error_regex = True

            if error_regex:

                msg_type = "The value(s) is(are) invalid for the qualifier %s of the feature %s. We will not report the qualifier. (Here is the regex expected: %s)"  % (qualifier, self.type, regex)
                msg = "The value(s) %s is(are) invalid for the qualifier %s of the feature %s. We will not report the qualifier. (Here is the regex expected: %s)"  % (value, qualifier, self.type, regex)
                self.handle_message("warning", msg_type, msg, None)

                return


        logging.debug("Adding value '%s' to qualifier '%s'" % (value, qualifier))

        if self.qualifier_prefix.get(gff_qualifier, None):
            value = ["%s%s" % (self.qualifier_prefix[gff_qualifier], v) for v in value]

        if self.qualifier_suffix.get(gff_qualifier, None):
            value = ["%s%s" % (v, self.qualifier_suffix[gff_qualifier]) for v in value]

        ###########################################
        # add the value only if not already present
        # List case
        if isinstance(value, list):
            for val in value:
                if val not in self.qualifiers[qualifier].value:
                    self.qualifiers[qualifier].add_value(val)
                else:
                    logging.debug('val %s alredy exist (list case)' %  val)
        # Scalar case
        else:
            if value not in self.qualifiers[qualifier].value:
                self.qualifiers[qualifier].add_value(value)
            else:
                logging.debug('val %s alredy exist (scalar case)' %  val)

    def combine(self, other):
        """
        Attempt to combine all features from another feature into this one.
        """

        # add new location
        self.location += other.location

        # combine qualifier except codon start
        for gff_qualifier, list_val_other in other.qualifiers.items():
            other_qualifier = self._from_gff_qualifier(gff_qualifier) # get the real qualifier name in EMBL format to be able to compare with the one alredy saved
            if other_qualifier != "codon_start":
                self.add_qualifier(gff_qualifier, list_val_other)
            else:
                # as the feature are sorted by increasing order location if + strand the first CDS codon_start qualifier was the good one
                # Of we are in a minus strand case we have to replace the start_codon, only hte last one will left
                if self.location.strand < 0:
                    # get phase of the last CDS
                    phase = int(other.qualifiers.get("phase", [0])[0])

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
                codon_table = CodonTable.ambiguous_dna_by_id[self.transl_table]
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
            codon_table = CodonTable.ambiguous_dna_by_id[self.transl_table]
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
        B = "Asx";  Aspartic acid (R) or Asparagine (N)
        X = "Xxx";  Unknown or 'other' amino acid
        Z = "Glx";  Glutamic acid (E) or Glutamine (Q)
        J = "Xle";  Leucine (L) or Isoleucine (I), used in mass-spec (NMR)
        U = "Sec";  Selenocysteine
        O = "Pyl";  Pyrrolysine
        """
        codon_table = CodonTable.ambiguous_dna_by_id[self.transl_table]
        seq = Seq(str(self.sequence()),IUPACAmbiguousDNA())

        #start translation according to the phase. Phase and codon_start are not the same coordinate system. It is why we have to remove 1
        phase = int(self.qualifiers.get('codon_start').value[0]) - 1
        if phase != 0:
            seq = seq[phase:]

        #check if multiple of three
        remaining_nuc = len(seq)%3
        if remaining_nuc != 0:
            #create warning
            ID=''
            for qualifier in self.feature.qualifiers:
                if 'id' == qualifier.lower():
                    ID =  "%s" % " ".join(self.feature.qualifiers[qualifier])
                    break
            logging.warning('Partial CDS. The CDS with ID = %s not a multiple of three.' %  ID)

        #translate the sequence in AA with normal frame even if stop inside
        translated_seq = seq.translate(codon_table).tostring().replace('B','X').replace('Z','X').replace('J','X')

        #Extra check about stop codon in CDS
        if '*' in translated_seq[:-1]: # check if premature stop codon in the translation
            ID=''
            for qualifier in self.feature.qualifiers:
                if 'id' == qualifier.lower():
                    ID =  "%s" % " ".join(self.feature.qualifiers[qualifier])
                    break
            logging.error('Stop codon found within the CDS (ID = %s). It will rise an error submiting the data to ENA. Please fix your gff file (check the phase Column 8).' %  ID)

        # remove the stop character. It's not accepted by embl
        if translated_seq[-1:] == "*":
            translated_seq = translated_seq[:-1]

        return translated_seq

    def handle_message(self, type, msg_type, msg, value):

        if msg_type in Feature.PREVIOUS_ERRORS:
            Feature.PREVIOUS_ERRORS[msg_type] += 1

        level = eval("logging.%s" % type.upper())

        if self.uncompressed_log:
            logging.log(level, msg)
        else:
            if not value:   # number of line accepted to display (defaut or given to the method)
                value = 5
            if msg_type not in Feature.PREVIOUS_ERRORS or Feature.PREVIOUS_ERRORS[msg_type] < value:
                logging.log(level, msg)
                Feature.PREVIOUS_ERRORS.setdefault(msg_type,1)
            elif Feature.PREVIOUS_ERRORS[msg_type] == value:
                logging.log(level, msg)
                final_message = 'We will not display anymore this %s. Please use the --uncompressed_log parameter if you wish having all of them.' % type
                logging.log(level, final_message)


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
            print(gff_feature)
            print("_"*80)
            feature = Feature( gff_feature, args.translation_file, 1, feature_definition_dir = "features", qualifier_definition_dir="qualifiers" )
            print("_"*80)
            print(feature)
            break
    except Exception as e:
        import traceback
        traceback.print_exc(limit=5)
        sys.stderr.write( str(e) )
