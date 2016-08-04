#!/usr/bin/python2.7
"""
EMBL writer for ENA data submission. Note that this implementation is basically 
just the documentation at ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt 
in python form - the implementation could be a lot more efficient!


TODO: find list of previous ENA release dates and numbers
TODO: find way to retrieve current release date
TODO: implement all features
TODO: implement proper data types for all feature qualifiers
TODO: make better guessing system
TODO: add reasonable way to add references
TODO: add helpful way to get classification via tax_id or scientific name
TODO: add way to handle mandatory features and feature qualifiers (especially contingent dependencies)


"""

import sys
import gzip
import time
import argparse
import curses.ascii
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from modules.feature import *

class EMBL( object ):
    """
    The basic structure of an EMBL file is like this:
    
    ID - identification             (begins each entry; 1 per entry)
    AC - accession number           (>=1 per entry)
    PR - project identifier         (0 or 1 per entry)
    DT - date                       (2 per entry)
    DE - description                (>=1 per entry)
    KW - keyword                    (>=1 per entry)
    OS - organism species           (>=1 per entry)
    OC - organism classification    (>=1 per entry)
    OG - organelle                  (0 or 1 per entry)
    RN - reference number           (>=1 per entry)
    RC - reference comment          (>=0 per entry)
    RP - reference positions        (>=1 per entry)
    RX - reference cross-reference  (>=0 per entry)
    RG - reference group            (>=0 per entry)
    RA - reference author(s)        (>=0 per entry)
    RT - reference title            (>=1 per entry)
    RL - reference location         (>=1 per entry)
    DR - database cross-reference   (>=0 per entry)
    CC - comments or notes          (>=0 per entry)
    AH - assembly header            (0 or 1 per entry)   
    AS - assembly information       (0 or >=1 per entry)
    FH - feature table header       (2 per entry)
    FT - feature table data         (>=2 per entry)    
    XX - spacer line                (many per entry)
    SQ - sequence header            (1 per entry)
    CO - contig/construct line      (0 or >=1 per entry) 
    bb - (blanks) sequence data     (>=1 per entry)
    // - termination line           (ends each entry; 1 per entry)
    
    """
    
    legal_values = {'data_class':{"CON":"Entry constructed from segment entry sequences; if unannotated, annotation may be drawn from segment entries",
                                  "PAT":"Patent",
                                  "EST":"Expressed Sequence Tag",
                                  "GSS":"Genome Survey Sequence",
                                  "HTC":"High Thoughput CDNA sequencing",
                                  "HTG":"High Thoughput Genome sequencing",
                                  "MGA":"Mass Genome Annotation",
                                  "WGS":"Whole Genome Shotgun",
                                  "TSA":"Transcriptome Shotgun Assembly",
                                  "STS":"Sequence Tagged Site",
                                  "STD":"Standard (all entries not classified as above)",
                                 },
                    'taxonomy':{"PHG":"Bacteriophage",
                                "ENV":"Environmental Sample",
                                "FUN":"Fungal",
                                "HUM":"Human",
                                "INV":"Invertebrate",
                                "MAM":"Other Mammal",
                                "VRT":"Other Vertebrate",
                                "MUS":"Mus musculus", 
                                "PLN":"Plant",
                                "PRO":"Prokaryote",
                                "ROD":"Other Rodent",
                                "SYN":"Synthetic",
                                "TGN":"Transgenic",
                                "UNC":"Unclassified",
                                "VRL":"Viral",
                               },
                    'topology':['linear', 'circular'],
                    'molecule_class':["genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", "other RNA", "other DNA", 
                                      "transcribed RNA", "viral cRNA", "unassigned DNA", "unassigned RNA"],
                    'organelle':["chromatophore", "hydrogenosome", "mitochondrion", "nucleomorph", "plastid", 
                                 "mitochondrion:kinetoplast", "plastid:chloroplast", "plastid:apicoplast", 
                                 "plastid:chromoplast", "plastid:cyanelle", "plastid:leucoplast", "plastid:proplastid"],
                    'transl_table':[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25],
#                    'transl_table':["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"],
                    }
    
    release_dates = {125:time.strptime("2015-09-23", "%Y-%m-%d"),
                     124:time.strptime("2015-07-01", "%Y-%m-%d"), 
                     123:time.strptime("2015-03-23", "%Y-%m-%d"),
                     122:time.strptime("2014-12-09", "%Y-%m-%d"),
                     121:time.strptime("2014-09-24", "%Y-%m-%d"),
                     120:time.strptime("2014-07-01", "%Y-%m-%d"),
                     119:time.strptime("2014-03-17", "%Y-%m-%d"),
                     118:time.strptime("2013-12-17", "%Y-%m-%d"),
                     117:time.strptime("2013-09-12", "%Y-%m-%d"),
                     116:time.strptime("2013-06-27", "%Y-%m-%d"),
                     114:time.strptime("2012-12-21", "%Y-%m-%d"),
                     }
    
    spacer = "\nXX"
    termination = "\n//\n"
    
    def __init__(self, record = None, verify = False, guess = True):
        self.record = record
        self.verify = verify
        self.guess = guess
        self.refs = []
        self.dbxref = []
        self.assembly_information = []
        self.construct_information = []
        self.comment = ""
        super(EMBL, self).__init__()
    
    def _add_mandatory(self):
        # Make sure that there's at least one source qualifier
        if not [f for f in self.record.features if f.type == 'source']:
            source_location = FeatureLocation(ExactPosition(0), ExactPosition(len(self.record.seq)))
            source_location.strand = 1
            source_feature = SeqFeature( source_location )
            source_feature.qualifiers["mol_type"] = self.molecule_class
            source_feature.qualifiers["organism"] = self.species
            source_feature.type = "source"
            self.record.features[0:0] = [source_feature]
        
        # Make sure that there's a gap feature for every span of n's
        start = None
        for i, c in enumerate(self.record.seq):
            if c in ['n', 'N']:
                if start == None:
                    start = i
            else:
                if start != None:
                    found = False
                    for f in [f for f in self.record.features if f.type == 'gap']:

                        if f.location.start == start and f.location.end == i:
                            found = True
                    if not found:
                        gap_location = FeatureLocation(ExactPosition(start), ExactPosition(i))
                        gap_location.strand = 1
                        gap_feature = SeqFeature( gap_location )
                        gap_feature.qualifiers["estimated_length"] = i-start
                        gap_feature.type = "gap"
                        self.record.features += [gap_feature]
                    
                start = None

    def _get_release(self, date):
        """
        Tries to find the correct release number for a given date.
        """
        
        previous = None
        for release in sorted(self.release_dates):
            if not previous:
                previous = self.release_dates[release]
                continue
            if date > previous and date < self.release_dates[release]:
                return release
        
        #TODO: find way to get latest release number and date!
        
        return max(self.release_dates.keys())+1
    
    @staticmethod
    def guess_qualifier(qualifier):
        if qualifier in ['Dbxref']:
            return "db_xref"
        if qualifier in ['description']:
            return "note"
        return False
    
    def _multiline(self, prefix, data, sep=";", suffix="", indent = 3, quoted = False):
        """
        Creates a multiline entry.
        
        If data is a list, the entries are listed with "sep" as separator, if data is
        a string, it's quoted over as many lines as needed.
        """
        
        #particular case when RT come empty. We must print ; wihtout quotes
        if(prefix == "RT" and data == ";"):
            output = "%s%s" % (prefix, " "*indent)
            output += str(data)
            return "\n" + output + suffix

        output = "%s%s" % (prefix, " "*indent)
        if type(data) == type([]):
            for item in data:
                if len(output) - output.rfind("\n") + len(item) > 80:
                    output += "\n%s%s" % (prefix, " "*indent)
                output += item + "%s " % sep
        else:
            string = " ".join(data.split("\n"))
            output += "\"" if quoted else ""
            while string:
                line = string[:80-(len(output) - output.rfind("\n"))]
                output += line
                string = string[len(line):]
                if string or (len(line)+indent+2) == 80:
                    output += "\n%s%s" % (prefix, " "*indent)
            output += "\"" if quoted else ""
                
        return "\n" + output.strip().strip(sep) + suffix
    
    def _print_feature(self, feature):
        output = "FT   %s %i..%i\n" % ("{:15}".format(feature.type), feature.location.start+1, feature.location.end)
        
        for key, values in feature.qualifiers.iteritems():
            key, values = EMBL.legal_qualifier(feature.type, key, values, self.guess)
            if not key:
                continue
            
            for value in values:
                line = "FT                   /%s=" % key
                if type(value) == type(0):
                    line += str(value) + "\n"
                else:
                    while value:
                        linebreak = 80-len(line)
                        line += value[:linebreak]
                        value = value[linebreak:]
                        if value:
                            output += line + "\n"
                            line = "FT                   "
                        else:
                            line += "\n"
                output += line
        
        for sub_feature in feature.sub_features:
            sub_feature.type = EMBL.legal_feature(sub_feature.type)
            if not sub_feature.type:
                continue
            output += self._print_feature(sub_feature)
        
        return output
    
    def _set_all(self):
        self.set_accession()
        self.set_classification()
        self.set_created()
        self.set_description()
        self.set_sample_id()
        self.set_project_id
        self.set_version()
        self.set_topology()
        self.set_transl_table()
        self.set_molecule_class()
        self.set_data_class()
        self.set_taxonomy()
    
    def _verify(self, key, key_type):
        if key_type not in self.legal_values:
            sys.stderr.write("Can't verify value for %s, legal values unknown." % key_type)
            return key
        legal_values = self.legal_values[key_type]
        while key not in self.legal_values[key_type]:
            sys.stderr.write("\n'%s' is not a legal value for %s.\n" % (key, key_type))
            sys.stderr.write("Legal values are:\n")
            if type(self.legal_values[key_type]) == type({}):
                for value, description in self.legal_values[key_type].iteritems():
                    sys.stderr.write("  - %s\t%s\n" % (value, description))
            else:
                for value in self.legal_values[key_type]:
                    sys.stderr.write("  - %s\n" % value)
            sys.stderr.write("Please enter new value: ")
            key = raw_input()
            
        return key
    
    def add_xref(self, xref):
        self.dbxref += [xref]
    
    def add_reference(self, title, positions = "all", location = "", comment = "", xrefs = [], group = [], authors = []):
        self.refs += [{'title':title,
                       'positions':positions if positions != 'all' else [(1,len(self.record.seq))],
                       'location':location,
                       'comment':comment,
                       'xrefs':xrefs,
                       'group':group,
                       'authors':authors}]
    
    @staticmethod
    def legal_feature(feature_type):
        while feature_type not in FeatureTable.features:
            sys.stderr.write("\n%s is not a legal feature type.\n" % feature_type)
            sys.stderr.write("Please enter new feature type: ")
            string = raw_input()
            feature_type = string
        return feature_type
    
    @staticmethod
    def legal_qualifier(feature_type, qualifier, values, guess = False):
        if not FeatureTable.features[feature_type]:
            return qualifier, values
        while qualifier and qualifier not in FeatureTable.features[feature_type]['qualifiers']:
            if guess:
                temp = EMBL.guess_qualifier(qualifier)
                if temp:
                    qualifier = temp
                else:
                    sys.stderr.write( "\nqualifier: '%s' with value(s) '%s' is illegal for feature '%s'\n" % (qualifier, ", ".join(values), feature_type) )
                    sys.stderr.write( "  Legal qualifiers are:\n" )
                    for key, value in FeatureTable.features[feature_type]['qualifiers'].iteritems():
                        sys.stderr.write( "  -  %s\n" % key )
                    sys.stderr.write( "Please enter new qualifier name (or blank to remove): ")
                    qualifier = raw_input()
        
        if not qualifier:
            return qualifier, values
        
        out_values = []
        if type(values) != type([]):
            values = [values]
        for value in values:
            if str in FeatureTable.features[feature_type]['qualifiers'][qualifier]:
                out_values += ["\"%s\"" % value]
            elif int in FeatureTable.features[feature_type]['qualifiers'][qualifier]:
                out_values += [int(value)]
            else:
                try:
                    value = int(value)
                except:
                    value = "\"%s\"" % value
                out_values += [value]
                
                #sys.stderr.write( "\nqualifier: '%s' is illegal for feature '%s'\n" % (qualifier, feature_type) )
                #sys.stderr.write( "Please enter new qualifier name (or blank to remove): ")
                #qualifier = raw_input()
        
        return qualifier, out_values
    
    def ID(self):
        """
        from ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt:
        
        The ID (IDentification) line is always the first line of an entry. The
        format of the ID line is:
        ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
        The tokens represent:
           1. Primary accession number
           2. Sequence version number
           3. Topology: 'circular' or 'linear'
           4. Molecule type (see note 1 below)
           5. Data class (see section 3.1)
           6. Taxonomic division (see section 3.2)
           7. Sequence length (see note 2 below)
        
        """
        
        if self.verify:
            self.topology       = self._verify( self.topology,       "topology")
            self.molecule_class = self._verify( self.molecule_class, "molecule_class")
            self.data_class     = self._verify( self.data_class,     "data_class")
            self.taxonomy       = self._verify( self.taxonomy,       "taxonomy")
        
        return "ID   %s; SV %s; %s; %s; %s; %s; %i BP." % (self.accessions[0], self.version, self.topology, 
                                                           self.molecule_class, self.data_class, self.taxonomy, 
                                                           len(self.record.seq) ) + self.spacer
    
    def AC(self):
        """
        from ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt:
        
        The AC (ACcession number) line lists the accession numbers associated with 
        the entry. Each accession number, or range of accession numbers, is terminated by a
        semicolon. Where necessary, more than one AC line is used.
        """
        
        output = "AC   "
        
        if len(self.accessions) == 0 and self.verify:
            sys.stderr.write("At least one accession number is needed: ")
            self.accessions += [raw_input()]
        
        for accession in self.accessions:
            if len(output) + len(accession) > 80:
                output += "\nAC   "
            output += accession + "; "
            
        return "\n" + output.strip() + self.spacer
    
    def PR(self):
        """
        The PR (PRoject) line shows the International Nucleotide Sequence Database
        Collaboration (INSDC) Project Identifier that has been assigned to the entry.
        Full details of INSDC Project are available at
        http://www.ebi.ac.uk/ena/about/page.php?page=project_guidelines.
        """
        return "\nPR   Project:%s;" % self.project_id + self.spacer
    
    def DT(self):
        """
        The DT (DaTe) line shows when an entry first appeared in the database and
        when it was last updated.  Each entry contains two DT lines, formatted
        as follows:
        DT   DD-MON-YYYY (Rel. #, Created)
        DT   DD-MON-YYYY (Rel. #, Last updated, Version #)
        
        the Release number (Rel.) indicates the first quarterly release made *after* 
        the entry was created or last updated.
        
        The Release number is in places like the header of ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
        but I have no idea where to get it most easily... 
        """
        
        updated = time.localtime() # this should be the latest update...
        
        output  = "\nDT   %s (Rel. %s, Created)" % (time.strftime("%d-%b-%Y", self.created), self._get_release(self.created))
        output += "\nDT   %s (Rel. %s, Last updated, Version %i)" % (time.strftime("%d-%b-%Y", updated), 
                                                                     self._get_release(updated), self.version)
        return output + self.spacer
    
    def DE(self):
        """
        The DE (Description) lines contain general descriptive information about the
        sequence stored. This may include the designations of genes for which the
        sequence codes, the region of the genome from which it is derived, or other
        information which helps to identify the sequence.
        """
        output = ""
        temp = str(self.description)
        while temp:
            output += "\nDE   %s" % temp[:75]
            temp = temp[75:]
        return output + self.spacer
    
    def KW(self):
        """
        The KW (KeyWord) lines provide information which can be used to generate
        cross-reference indexes of the sequence entries based on functional,
        structural, or other categories deemed important.
        """
        
        output = "KW   "
        
        if len(self.keywords) == 0 and self.verify:
            sys.stderr.write("At least one keyword is needed: ")
            self.keywords += [raw_input()]
        
        return self._multiline("KW", self.keywords, suffix=".") + self.spacer
    
    def OS(self):
        """
        The OS (Organism Species) line specifies the preferred scientific name of
        the organism which was the source of the stored sequence. In most 
        cases this is done by giving the Latin genus and species designations, 
        followed (in parentheses) by the preferred common name in English where
        known. The preferred format is:
             OS   Genus species (name)
        """
        return "\nOS   %s" % self.species + self.spacer
    
    def OC(self):
        """
        The OC (Organism Classification) lines contain the taxonomic classification
        of the source organism as described in Section 2.2 of 
        ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
        """
        
        output = "OC   "
        
        if len(self.classification) == 0 and self.verify:
            sys.stderr.write("At least one classification level is needed: ")
            self.classification += [raw_input()]
        
        return self._multiline("OC", self.classification, suffix=";") + self.spacer
    
    def OG(self):
        """
        The OG (OrGanelle) linetype indicates the sub-cellular location of non-nuclear
        sequences.  It is only present in entries containing non-nuclear sequences
        and appears after the last OC line in such entries.
        
        The OG line contains
        a) one data item (title cased) from the controlled list detailed under the
        /organelle qualifier definition in the Feature Table Definition document
        that accompanies this release or
        b) a plasmid name.
        Examples include "Mitochondrion", "Plastid:Chloroplast" and "Plasmid pBR322".
        
        legal values from http://www.insdc.org/controlled-vocabulary-organelle-qualifier
        """
        
        if self.organelle:
            return ("OG   %s %s" % (self.organelle, self.organelle_name)).strip() + "\n" + self.spacer
        return ""
    
    def RF(self):
        """
        The Reference (RN, RC, RP, RX, RG, RA, RT, RL) Lines
        These lines comprise the literature citations within the database.
        The citations provide access to the papers from which the data has been 
        abstracted. The reference lines for a given citation occur in a block, and
        are always in the order RN, RC, RP, RX, RG, RA, RT, RL. 
        
        >>>> Within each such reference block <<<<<

        the RN line occurs once, the RC, RP and RX lines occur zero
        or more times, and the RA, RT, RL lines each occur one or more times. 
        If several references are given, there will be a reference block for each.
        """
        
        output = ""
        
        for i, ref in enumerate(self.refs):
            output += "\nRN   [%i]" % (i+1)                         # RN - reference number           (>=1 per entry)
            if ref['comment']:                                      # RC - reference comment          (>=0 per entry)
                output += self._multiline("RC", ref['comment'])
                                                                    # RP - reference positions        (>=1 per entry)
            output += self._multiline("RP", ["%i-%i" % pos for pos in ref['positions']])   
            for xref in ref['xrefs']:                               # RX - reference cross-reference  (>=0 per entry)
                output += self._multiline("RX", xref, suffix=".")
            if ref['group']:                                        # RG - reference group            (>=0 per entry)
                output += self._multiline("RG", ref['group'])
            if ref['authors']:                                      # RA - reference author(s)        (>=0 per entry)
                output += self._multiline("RA", ref['authors'], sep=", ", suffix=";")
             
            if ref['title'] == ";":                                 # RT - reference title            (>=1 per entry)
                output += self._multiline("RT", ref['title'], quoted=False)
            else:
                output += self._multiline("RT", ref['title'], quoted=True)
            # TODO: There are lots of recommended formatting for the references,
            #       but I won't bother implementing them right now.
            output += self._multiline("RL", ref['location'])        # RL - reference location         (>=1 per entry)
        
        if output:
            return output + self.spacer
        
        return ""
    
    def DR(self):
        """
        The DR (Database Cross-reference) line cross-references other databases which
        contain information related to the entry in which the DR line appears. For
        example, if an annotated/assembled sequence in ENA is cited in the IMGT/LIGM
        database there will be a DR line pointing to the relevant IMGT/LIGM entry.
        """
        if self.dbxref:
            output = ""
            for xref in self.dbxref:
                output += self._multiline("DR", xref, suffix=".")
            return output + self.spacer
        return ""
    
    def CC(self):
        return self._multiline("CC", self.comment, quoted=True) + self.spacer if self.comment else ""
    
    def AH(self):
        """
        Third Party Annotation (TPA) and Transcriptome Shotgun Assembly (TSA) records
        may include information on the composition of their sequences to show
        which spans originated from which contributing primary sequences. The AH
        (Assembly Header) line provides column headings for the assembly information.
        The lines contain no data and may be ignored by computer programs.
        """
        return "AH   LOCAL_SPAN     PRIMARY_IDENTIFIER     PRIMARY_SPAN     COMP"
    
    def AS(self):
        """
        The AS (ASsembly Information) lines provide information on the composition of 
        a TPA or TSA sequence. These lines include information on local sequence spans
        (those spans seen in the sequence of the entry showing the AS lines) plus
        identifiers and base spans of contributing primary sequences (for ENA
        primary entries only).
    
        a) LOCAL_SPAN               base span on local sequence shown in entry  
        b) PRIMARY_IDENTIFIER       acc.version of contributing ENA sequence(s)
                                    or trace identifier for ENA read(s)
        c) PRIMARY_SPAN             base span on contributing ENA primary
                                    sequence or not_available for ENA read(s)
                                   
        d) COMP                     'c' is used to indicate that contributing sequence
                                    originates from complementary strand in primary
                                    entry
        """
        output = ""
        for assembly in self.assembly_information:
            output += "AS   %s%s%s%s" % ("{:16}".format(assembly['local_span']),
                                         "{:24}".format(assembly['identifier']),
                                         "{:18}".format(assembly['primary_span']),
                                         assembly['complementary'])
        return output
    
    def FH(self):
        """
        The FH (Feature Header) lines are present only to improve readability of
        an entry when it is printed or displayed on a terminal screen. The lines 
        contain no data and may be ignored by computer programs. The format of these
        lines is always the same.
        """
        return "\nFH   Key             Location/Qualifiers\nFH"
    
    def FT(self):
        """
        The FT (Feature Table) lines provide a mechanism for the annotation of the
        sequence data. Regions or sites in the sequence which are of interest are
        listed in the table. In general, the features in the feature table represent
        signals or other characteristics reported in the cited references. In some
        cases, ambiguities or features noted in the course of data preparation have 
        been included.  The feature table is subject to expansion or change as more
        becomes known about a given sequence.
        
        Feature Table Definition Document:
        A complete and definitive description of the feature table is given 
        in the document "The DDBJ/ENA/GenBank Feature Table:  Definition". 
        URL: ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.txt
        """
        #add extra information to CDS features
        for feature in self.record.features:
            if not feature.type.lower() == "source":
                for feature_l2 in feature.sub_features:
                    for l3_feature in feature_l2.sub_features:
                        if(l3_feature.type.lower() == "cds"):
                            l3_feature.qualifiers['transl_table']=self.transl_table

        #process features
        output = ""
        for feature in self.record.features:
             for feat in parse_gff_feature(self.accessions, feature, self.record.seq):
                 output += str(feat)
        return output + self.spacer
    
    def CO(self):
        """
        Con(structed) sequences in the CON data classes represent complete
        chromosomes, genomes and other long sequences constructed from segment entries.
        CON data class entries do not contain sequence data per se, but rather the
        assembly information on all accession.versions and sequence locations relevant
        to building the constructed sequence. The assembly information is represented in
        the CO lines.
        """
    
    def SQ(self, out = None):
        """
        The SQ (SeQuence header) line marks the beginning of the sequence data and 
        Gives a summary of its content.
        
        This is followed by sequence data lines, which have a line code consisting 
        of two blanks. The sequence is written 60 bases per line, in groups of 10 
        bases separated by a blank character, beginning at position 6 of the line. 
        The direction listed is always 5' to 3', and wherever possible the 
        non-coding strand (homologous to the message) has been stored.
        
        This function can be passed a streamhandler to write directly to, instead of 
        buffering the entire sequence.
        """
        seq = str(self.record.seq) if self.record else ""
        num_a = seq.count("a")+seq.count("A")
        num_c = seq.count("c")+seq.count("C")
        num_g = seq.count("g")+seq.count("G")
        num_t = seq.count("t")+seq.count("T")
        num_o = len(seq) - (num_a + num_c + num_g + num_t)
        
        output = "\nSQ   Sequence %i BP; %i A; %i C; %i G; %i T; %i other;" % (len(seq), num_a, num_c, num_g, num_t, num_o)
        
        if out:
            out.write(output)
            output = ""
        
        seq_len = 0
        while seq:
            current_line = " ".join([seq[i*10:(i+1)*10] for i in range(0, 6)])
            seq_len += min(60, len(seq))
            formatted_line = "\n     %s %s" % ("{:65}".format(current_line), "{:>9}".format(str(seq_len)))
            
            if out:
                out.write(formatted_line)
            else:
                output += formatted_line
            seq = seq[60:]
        return output
    
    def set_accession(self, accessions = []):
        """
        Sets the entry accession numbers, or parses it from the current record
        """
        if not hasattr(self, "accessions"):
            self.accessions = []
        if accessions:
            self.accessions = accessions
        if hasattr(self.record, "dbxrefs"):
            self.accessions += self.record.dbxrefs
        if not getattr(self, "accessions", False):
            self.accessions = ["Unknown"]
    
    def set_classification(self, classification = []):
        """
        Sets the entry phylogenetic classification, or parses it from the current record
        """
        if classification:
            self.classification = classification
        elif hasattr(self.record, "classification"):
            self.classification += self.record.classification
        if not getattr(self, "classification", False):
            self.classification = [""]
    
    def set_created(self, timestamp = None):
        """
        Sets the creation time of the original entry.
        """
        if timestamp:
            self.created = timestamp
        elif hasattr(self.record, "created"):
            self.created = self.record.created
        elif not hasattr(self, "created"):
            self.created = time.localtime()
    
    def set_data_class(self, data_class = None):
        """
        Sets the sample data class, or parses it from the current record.
        """
        if data_class:
            self.data_class = data_class
        elif hasattr(self.record, "data_class"):
            self.data_class = self.record.data_class
        elif not hasattr(self, "data_class"):
            self.data_class = ""
        
        if self.verify:
            self.data_class = self._verify( self.data_class, "data_class")
    
    def set_description(self, description = None):

        """
        Sets the sample description, or parses it from the current record.
        """
        if description:
            self.description = description
        elif hasattr(self.record, "description"):
            self.description = self.record.description
        elif not hasattr(self, "description"):
            self.description = ""
    
    def set_keywords(self, keywords = []):
        """
        Sets the entry keywords, and parses those of the current record
        """
        if keywords:
            self.keywords = keywords
        if hasattr(self.record, "keywords"):
            self.keywords += self.record.keywords
        if not getattr(self, "keywords", False):
            self.keywords = [""]
    
    def set_molecule_class(self, molecule_class = None):
        """
        Sets the sample molecule class, or parses it from the current record.
        """
        if molecule_class:
            self.molecule_class = molecule_class
        elif hasattr(self.record, "molecule_class"):
            self.molecule_class = self.record.molecule_class
        elif not hasattr(self, "molecule_class"):
            self.molecule_class = ""
        
        if self.verify:
            self.molecule_class = self._verify( self.molecule_class, "molecule_class")
    
    def set_organelle(self, organelle = None, organelle_name = None):
        """
        Sets the sample organelle, and organelle name, or parses them from the current record.
        """
        if organelle:
            self.organelle = organelle
        elif hasattr(self.record, "organelle"):
            self.organelle = self.record.organelle
        elif not hasattr(self, "organelle"):
            self.organelle = ""
        
        if organelle_name:
            self.organelle_name = organelle_name
        elif hasattr(self.record, "organelle name"):
            self.organelle_name = organelle_name
        elif not hasattr(self, "organelle_name"):
            self.organelle_name = ""
        
        if self.verify and self.organelle:
            self.organelle = self._verify( self.organelle, "organelle")
    
    def set_project_id(self, project_id = None):
        """
        Sets the project id, or parses it from the current record
        """
        if project_id:
            self.project_id = project_id
        elif hasattr(self.record, "project_id"):
            self.project_id = self.record.project_id
        elif not hasattr(self, "project_id"):
            self.project_id = "Unknown"
    
    def set_record(self, record):
        self.record = record
    
    def set_sample_id(self, sample_id = None):
        """
        Sets the sample id, or parses it from the current record
        """
        if sample_id:
            self.sample_id = sample_id
        elif hasattr(self.record, "id"):
            self.sample_id = self.record.id
        elif not hasattr(self, "sample_id"):
            self.sample_id = "Unknown"
    
    def set_species(self, species = None):
        """
        Sets the species, or parses it from the current record
        """
        if species:
            self.species = species
        elif hasattr(self.record, "species"):
            self.species = self.record.species
        if not getattr(self, "species", False):
            self.species = "Unknown (Unknown)"
        
    def set_taxonomy(self, taxonomy = None):
        """
        Sets the sample taxonomy, or parses it from the current record.
        """
        if taxonomy:
            self.taxonomy = taxonomy
        elif hasattr(self.record, "taxonomy"):
            self.taxonomy = self.record.taxonomy
        elif not hasattr(self, "taxonomy"):
            self.taxonomy = ""
        
        if self.verify:
            self.taxonomy = self._verify( self.taxonomy, "taxonomy")
    
    def set_topology(self, topology = None):
        """
        Sets the sample topology, or parses it from the current record.
        """
        if topology:
            self.topology = topology
        elif hasattr(self.record, "topology"):
            self.topology = self.record.topology
        elif not hasattr(self, "topology"):
            self.topology = ""
        
        if self.verify:
            self.topology = self._verify( self.topology,       "topology")

    def set_transl_table(self, transl_table = None):
        """
        Sets the translation table, or parses it from the current record.
        """
        if transl_table:
            self.transl_table = transl_table
        elif hasattr(self.record, "transl_table"):
            self.transl_table = self.record.transl_table
        elif not hasattr(self, "transl_table"):
            self.transl_table = ""
        
        if self.verify:
            self.transl_table = self._verify( self.transl_table,       "transl_table")   
    def set_version(self, version = None):
        """
        Sets the release version, or parses it from the current record.
        """
        if version:
            self.version = version
        elif hasattr(self.record, "version"):
            self.version = self.record.version
        elif not hasattr(self, "version"):
            self.version = 1
    
    def write_all(self, out = sys.stdout):
        
        self._set_all()
        
        if type(out) == type(""):
            out = open(out, 'w')
        
        # Add missing mandatory features:
        self._add_mandatory()
        
        out.write( self.ID() ) # ID - identification             (begins each entry; 1 per entry)
        out.write( self.AC() ) # AC - accession number           (>=1 per entry)
        out.write( self.PR() ) # PR - project identifier         (0 or 1 per entry)
        out.write( self.DT() ) # DT - date                       (2 per entry)
        out.write( self.DE() ) # DE - description                (>=1 per entry)
        out.write( self.KW() ) # KW - keyword                    (>=1 per entry)
        out.write( self.OS() ) # OS - organism species           (>=1 per entry)
        out.write( self.OC() ) # OC - organism classification    (>=1 per entry)
        out.write( self.OG() ) # OG - organelle                  (0 or 1 per entry)
        out.write( self.RF() ) # References, including RN, RC, RP, RX, RG, RA, RT and RL
        out.write( self.DR() ) # DR - database cross-reference   (>=0 per entry)
        out.write( self.CC() ) # CC - comments or notes          (>=0 per entry)
        if self.assembly_information:
            out.write( self.AH() ) # AH - assembly header            (0 or 1 per entry)
            out.write( self.AS() ) # AS - assembly information       (0 or >=1 per entry)
        
        if self.record and self.record.features:
            out.write( self.FH() ) # FH - feature table header       (2 per entry)
        out.write( self.FT() ) # FT - feature table data         (>=2 per entry)
        if self.construct_information:
            out.write( self.CO() )
        
        self.SQ( out )         # SQ - sequence header            (1 per entry)
                               # + sequence...
        
        out.write( self.termination ) # // - termination line    (ends each entry; 1 per entry)
        
        
        # CO - contig/construct line      (0 or >=1 per entry) ?????
    
##########################
#        MAIN            #
##########################

if __name__ == '__main__':
    
    sys.stderr.write( """
    ########################################################
    # NBIS 2016 - Sweden                                   #  
    # Authors: Martin Norling, Jacques Dainat              #
    # Please cite NBIS (www.nbis.se) when using this tool. #
    ########################################################\n\n""")

    parser = argparse.ArgumentParser( description = __doc__ )
    
    parser.add_argument("gff_file", help="input gff-file")
    parser.add_argument("fasta", help="input fasta sequence")
    parser.add_argument("-a", "--accession", default=[], nargs="+", help="Accession number(s) for the entry")
    parser.add_argument("-c", "--created", default=None, help="Creation time of the original entry")
    parser.add_argument("-d", "--data_class", default="STD", help="Data class of the sample.", choices=["CON", "PAT", "EST", "GSS", "HTC", "HTG", "MGA", "WGS", "TSA", "STS", "STD"])
    parser.add_argument("-g", "--organelle", default=None, help="Sample organelle.")
    
    parser.add_argument("-i", "--sample_id", default=None, help="Sample ID, as returned from ENA.")
    parser.add_argument("-k", "--keyword", default=[], nargs="+", help="Keywords for the entry")
    parser.add_argument("-l", "--classification", default=["Life"], nargs="+", help="Phylogenetic classification")
    parser.add_argument("-m", "--molecule_class", default="genomic DNA", help="Molecule class of the sample.", choices=["genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", "other RNA", "other DNA", "transcribed RNA", "viral cRNA", "unassigned DNA", "unassigned RNA"])
    parser.add_argument("-o", "--output", default=None, help="output filename.")
    parser.add_argument("-p", "--project_id", default=None, help="Project ID (optional).")
#    parser.add_argument("-r", "--table", default="1", help="Translation table.", choices=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24",25"])
    parser.add_argument("-r", "--table", type=int, default=1, help="Translation table.", choices=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])
    parser.add_argument("-s", "--species", default=None, help="Sample Species, formatted as 'Genus species (english name)'.")
    parser.add_argument("-t", "--topology", default="linear", help="Sequence topology.", choices=[None, "linear", "circular"])
    
    parser.add_argument("--rc", default="", help="Reference Comment.")
    parser.add_argument("--rx", default="", help="Reference cross-reference.")
    parser.add_argument("--rg", default="", help="Reference Group.")
    parser.add_argument("--ra", "--author", nargs="+", default="", help="Author for the reference")
    parser.add_argument("--rt", default=";", help="Reference Title.")
    parser.add_argument("--rl", default="Unpublished", help="Reference publishing location.")
    
    
    parser.add_argument("-v", "--version", default=1, type=int, help="Sequence version number")
    parser.add_argument("-x", "--taxonomy", default=None, help="Source taxonomy.", choices=["PHG", "ENV", "FUN", "HUM", "INV", "MAM", "VRT", "MUS", "PLN", "PRO", "ROD", "SYN", "TGN", "UNC", "VRL"])
    
    parser.add_argument("--organelle_name", default=None, help="Sample organelle name.")
    
    args = parser.parse_args()
    
    if args.output:
        outfile = args.output
        if not outfile.endswith(".embl.gz"):
            outfile += ".gz" if outfile.endswith(".embl") else ".embl.gz"
        outfile = gzip.open(args.output, "wb")
    else:
        outfile = sys.stdout
    
    infile = gzip.open(args.gff_file) if args.gff_file.endswith(".gz") else open(args.gff_file)
    infasta = gzip.open(args.fasta) if args.fasta.endswith(".gz") else open(args.fasta)
    seq_dict = SeqIO.to_dict( SeqIO.parse(infasta, "fasta") )
    
    #Manage reference(s)
    referenceExist=""
    if not args.rg :
        sys.stderr.write( """It is not mandatory to add a reference. Yes I know it'obvious when it's an unpublished work and you don't plan to. 
So we will fill all infromation expected by ENA for unpublished work.
BUT it's always interesting to know who produced the record. So, please enter a Reference Group (or leave it empty) and press ENTER:\n""")
        key = raw_input()
        if key == "":
            print "It's a pity... Let's go anyway !"
        else:
            args.rg="key"

    for record in GFF.parse(infile, base_dict=seq_dict):
        
        writer = EMBL( record, True )
        writer.set_accession( args.accession )
        writer.set_created( args.created )
        writer.set_classification( args.classification )
        writer.set_data_class( args.data_class )
        writer.set_keywords( args.keyword )
        writer.set_molecule_class( args.molecule_class )
        writer.set_organelle( args.organelle, args.organelle_name )
        writer.set_project_id( args.project_id )
        writer.set_sample_id( args.sample_id )
        writer.set_species( args.species )
        writer.set_taxonomy( args.taxonomy )
        writer.set_topology( args.topology )
        writer.set_transl_table( args.table )
        writer.set_version( args.version )
        writer.add_reference(args.rt, location = args.rl, comment = args.rc, xrefs = args.rx, group = args.rg, authors = args.ra)

        writer.write_all( outfile )
     
    sys.stderr.write( """Well done\n""")