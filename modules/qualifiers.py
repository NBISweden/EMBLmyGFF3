#!/usr/bin/env python2.7

class Qualifier( object ):
    
    comment = None
    
    example = None
    
    regexp = None
    
    def __init__(self, value):
        self.value = value
    
    def add(self, value):
        if type(self.value) != type([]):
            self.value = [self.value]
        self.value += value if type(value) == type([]) else [value]
    
    # length of a line:79 characters
    def format(self, qualifier, value):
        output = ""
        if type(value) != type([]):
            value = [value]
        for val in value:
            if type(val) == type(""):
                val = "\"%s\"" % val
            line = "\nFT                   /%s=%s" % (qualifier, val)
            if len(line) <= 79:
                output += line
            else:
                output += line[:79]
                line = line[79:]
                while line:
                    output += "\nFT                   %s" % line[:59]
                    line = line[59:]
        return output

class AlleleQualifier( Qualifier ):
    
    definition = "name of the allele for the given gene."
    
    value_format = "text"
    
    comment = """all gene-related features (exon, CDS etc) for a given 
                 gene should share the same /allele qualifier value; 
                 the /allele qualifier value must, by definition, be 
                 different from the /gene qualifier value; when used with 
                 the variation feature key, the allele qualifier value 
                 should be that of the variant."""
    
    def __init__(self, value):
        super(AlleleQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("allele", self.value)

class AltitudeQualifier( Qualifier ):
    
    definition = "geographical altitude of the location from which the sample was collected."
    
    value_format = "text"
    
    comment = """Values indicate altitudes above or below nominal sea level 
                 provided in metres."""
    
    regexp = "-?[0-9]+\.?[0-9]* m"
    
    def __init__(self, value):
        super(AltitudeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("altitude", self.value)

class AnticodonQualifier( Qualifier ):
    
    definition = "location of the anticodon of tRNA and the amino acid for which it codes."
    
    value_format = """(pos:<location>,aa:<amino_acid>,seq:<text>) where location
                      is the position of the anticodon and amino_acid is the abbreviation for the
                      amino acid encoded and seq is the sequence of the anticodon"""
    
    def __init__(self, value):
        super(AnticodonQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("anticodon", self.value)

class Artificial_locationQualifier( Qualifier ):
    
    definition = """ indicates that location of the CDS or mRNA is modified to adjust
                     for the presence of a frameshift or internal stop codon and not
                     because of biological processing between the regions.."""
    
    value_format = """\"heterogeneous population sequenced\", \"low-quality sequence region\""""
    
    comment = """expected to be used only for genome-scale annotation."""
    
    def __init__(self, value):
        super(Artificial_locationQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("artificial_location", self.value)

class Bio_materialQualifier( Qualifier ):
    
    definition = """identifier for the biological material from which the nucleic
                    acid sequenced was obtained, with optional institution code and
                    collection code for the place where it is currently stored."""
    
    value_format = "[<institution-code>:[<collection-code>:]]<material_id>"
    
    comment = """the bio_material qualifier should be used to annotate the
                 identifiers of material in biological collections that are not
                 appropriate to annotate as either /specimen_voucher or
                 /culture_collection; these include zoos and aquaria, stock
                 centres, seed banks, germplasm repositories and DNA banks;
                 material_id is mandatory, institution_code and collection_code
                 are optional; institution code is mandatory where collection
                 code is present; institution code and collection code are taken
                 from a controlled vocabulary maintained by the INSDC."""
    
    def __init__(self, value):
        super(Bio_materialQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("bio_material", self.value)

class Bound_moietyQualifier( Qualifier ):
    
    definition = """name of the molecule/complex that may bind to the given feature."""
    
    value_format = "text"
    
    comment = """A single /bound_moiety qualifier is legal on the "misc_binding", 
                 "oriT" and "protein_bind" features."""
    
    def __init__(self, value):
        super(Bound_moietyQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("bound_moiety", self.value)

class Cell_lineQualifier( Qualifier ):
    
    definition = """cell line from which the sequence was obtained."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(Cell_lineQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("cell_line", self.value)

class Cell_typeQualifier( Qualifier ):
    
    definition = """cell type from which the sequence was obtained."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(Cell_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("cell_type", self.value)

class ChromosomeQualifier( Qualifier ):
    
    definition = """chromosome (e.g. Chromosome number) from which
                    the sequence was obtained."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(ChromosomeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("chromosome", self.value)

class CitationQualifier( Qualifier ):
    
    definition = """reference to a citation listed in the entry reference field."""
    
    value_format = """[integer-number] where integer-number is the number of the
                      reference as enumerated in the reference field"""
    
    comment = """used to indicate the citation providing the claim of and/or
                 evidence for a feature; brackets are used for conformity."""
    
    def __init__(self, value):
        super(CitationQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("citation", self.value)

class CloneQualifier( Qualifier ):
    
    definition = """clone from which the sequence was obtained."""
    
    value_format = "text"
    
    comment = """not more than one clone should be specified for a given source 
                 feature;  to indicate that the sequence was obtained from
                 multiple clones, multiple source features should be given.."""
    
    def __init__(self, value):
        super(CloneQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("clone", self.value)

class Clone_libQualifier( Qualifier ):
    
    definition = """clone library from which the sequence was obtained."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(Clone_libQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("clone_lib", self.value)

class Codon_startQualifier( Qualifier ):
    
    definition = """clone library from which the sequence was obtained."""
    
    value_format = "text"
    
    comment = """indicates the offset at which the first complete codon of a
                 coding feature can be found, relative to the first base of that
                 feature."""
    
    regexp = "[123]"
    
    def __init__(self, value):
        if type(value) == type([]):
            for i, val in enumerate(value):
                if type(val) == type(""):
                    val = int(val)
                val += 1
                value[i] = val
        super(Codon_startQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("codon_start", self.value)

class Collected_byQualifier( Qualifier ):
    
    definition = """name of persons or institute who collected the specimen."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(Collected_byQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("collected_by", self.value)

class Collection_dateQualifier( Qualifier ):
    
    definition = """The date on which the specimen was collected.
                    Date/time ranges are supported by providing two collection dates from among the supported value
                    formats, delimited by a forward-slash character.
                    Collection times are supported by adding "T", then the hour and minute and seconds, after the date.
                    Collection times must be in Coordinated Universal Time (UTC), otherwise known as "Zulu Time" (Z)."""
    
    value_format = """\"DD-Mmm-YYYY\", \"Mmm-YYYY\", \"YYYY\", \"YYYY-MM-DDThh:mmZ\", 
                      \"YYYY-MM-DDThh:mm:ssZ\", \"YYYY-MM-DDThhZ\", \"YYYY-MM-DD\", or \"YYYY-MM\""""
    
    comment = """'Mmm' represents a three-letter month abbreviation, and can be one of the following:
                 Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec

                 'YYYY' is a four-digit value representing the year. 'MM' is a two-digit value representing
                 the month. 'DD' is a two-digit value representing the day of the month.

                 'hh' is a two-digit value representing the hour of the day (00 to 23)
                 'mm' is a two-digit value representing the minute of the hour (00 to 59)
                 'ss' is a two-digit value representing the second of the hour (00 to 59)

                 Within a date range, value formats that make use of 'Mmm' (month abbreviations) cannot be
                 combined with value formats that make use of 'MM' (two-digit month number)

                 Collection dates that are specified to at least the month, day, and year (DD-Mmm-YYYY or YYYY-MM-DD)
                 are strongly encouraged. If the day and/or month of the collection date are not known, 
                 Mmm-YYYY or YYYY-MM or YYYY may be used.

                 Within a collection date range, the first date (possibly including time) must be
                 prior to the second date (possibly including time).

                 Within a collection date range for which the day, month, and year are identical, the first time value
                 must be prior to the second time value."""
    
    def __init__(self, value):
        super(Collection_dateQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("collection_date", self.value)

class CompareQualifier( Qualifier ):
    
    definition = """Reference details of an existing public INSD entry 
                    to which a comparison is made."""
    
    value_format = "[accession-number.sequence-version]"
    
    comment = """This qualifier may be used on the following features:
                 misc_difference, unsure, old_sequence and variation.
                 The feature "old_sequence" must have either a
                 /citation or a /compare qualifier. Multiple /compare
                 qualifiers with different contents are allowed within a 
                 single feature. 
                 This qualifier is not intended for large-scale annotation 
                 of variations, such as SNPs."""
    
    def __init__(self, value):
        super(CompareQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("compare", self.value)

class CountryQualifier( Qualifier ):
    
    definition = """locality of isolation of the sequenced organism indicated in
                    terms of political names for nations, oceans or seas, followed
                    by regions and localities."""
    
    value_format = """"\"<country_value>[:<region>][, <locality>]\" where 
                       country_value is any value from the controlled vocabulary at 
                       http://www.insdc.org/documents/country-qualifier-vocabulary"""
    
    comment = """Intended to provide a reference to the site where the source
                 organism was isolated or sampled. Regions and localities should
                 be indicated where possible. Note that the physical geography of
                 the isolation or sampling site should be represented in
                 /isolation_source."""
    
    def __init__(self, value):
        super(CountryQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("country", self.value)

class CultivarQualifier( Qualifier ):
    
    definition = """cultivar (cultivated variety) of plant from which sequence was 
                    obtained."""
    
    value_format = "text"
    
    comment = """'cultivar' is applied solely to products of artificial 
                 selection;  use the variety qualifier for natural, named 
                 plant and fungal varieties;"""
    
    def __init__(self, value):
        super(CultivarQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("cultivar", self.value)

class Culture_collectionQualifier( Qualifier ):
    
    definition = """institution code and identifier for the culture from which the
                    nucleic acid sequenced was obtained, with optional collection
                    code."""
    
    value_format = "<institution-code>:[<collection-code>:]<culture_id>"
    
    comment = """the /culture_collection qualifier should be used to annotate
                 live microbial and viral cultures, and cell lines that have been
                 deposited in curated culture collections; microbial cultures in
                 personal or laboratory collections should be annotated in strain
                 qualifiers;
                 annotation with a culture_collection qualifier implies that the
                 sequence was obtained from a sample retrieved (by the submitter
                 or a collaborator) from the indicated culture collection, or
                 that the sequence was obtained from a sample that was deposited
                 (by the submitter or a collaborator) in the indicated culture
                 collection; annotation with more than one culture_collection
                 qualifier indicates that the sequence was obtained from a sample
                 that was deposited (by the submitter or a collaborator) in more
                 than one culture collection.
                 culture_id and institution_code are mandatory, collection_code
                 is optional; institution code and collection code are taken
                 from a controlled vocabulary maintained by the INSDC.
                 http://www.insdc.org/controlled-vocabulary-culturecollection-qualifier"""
    
    def __init__(self, value):
        super(Culture_collectionQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("culture_collection", self.value)

class Db_xrefQualifier( Qualifier ):
    
    definition = """database cross-reference: pointer to related information in 
                    another database."""
    
    value_format = """"<database:identifier>" where database is
                   the name of the database containing related information, and 
                   identifier is the internal identifier of the related information
                   according to the naming conventions of the cross-referenced 
                   database."""
    
    comment = """the complete list of allowed database types is kept at 
                 http://www.insdc.org/db_xref.html"""
    
    def __init__(self, value):
        super(Db_xrefQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("db_xref", self.value)

class Dev_stageQualifier( Qualifier ):
    
    definition = """if the sequence was obtained from an organism in a specific 
                    developmental stage, it is specified with this qualifier."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(Dev_stageQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("dev_stage", self.value)

class DirectionQualifier( Qualifier ):
    
    definition = "direction of DNA replication"
    
    value_format = """left, right, or both where left indicates toward the 5' end of
                      the entry sequence (as presented) and right indicates toward
                      the 3' end"""
    
    def __init__(self, value):
        super(DirectionQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("direction", self.value)

class EC_numberQualifier( Qualifier ):
    
    definition = "Enzyme Commission number for enzyme product of sequence"
    
    value_format = "text"
    
    comment = """valid values for EC numbers are defined in the list prepared by the 
                 Nomenclature Committee of the International Union of Biochemistry and
                 Molecular Biology (NC-IUBMB) (published in Enzyme Nomenclature 1992,
                 Academic Press, San Diego, or a more recent revision thereof). 
                 The format represents a string of four numbers separated by full
                 stops; up to three numbers starting from the end of the string can 
                 be replaced by dash "." to indicate uncertain assignment. 
                 Symbol "n" can be used in the last position instead of a number 
                 where the EC number is awaiting assignment. Please note that such
                 incomplete EC numbers are not approved by NC-IUBMB."""
    
    def __init__(self, value):
        super(EC_numberQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("EC_number", self.value)

class EcotypeQualifier( Qualifier ):
    
    definition = """a population within a given species displaying genetically 
                    based, phenotypic traits that reflect adaptation to a local habitat."""
    
    value_format = "text"
    
    comment = """an example of such a population is one that has adapted hairier
                 than normal leaves as a response to an especially sunny habitat.
                 'Ecotype' is often applied to standard genetic stocks of
                 Arabidopsis thaliana, but it can be applied to any sessile 
                 organism."""
    
    def __init__(self, value):
        super(EcotypeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("ecotype", self.value)

class Environmental_sampleQualifier( Qualifier ):
    
    definition = """identifies sequences derived by direct molecular
                    isolation from a bulk environmental DNA sample
                    (by PCR with or without subsequent cloning of the
                    product, DGGE, or other anonymous methods) with no
                    reliable identification of the source organism.
                    Environmental samples include clinical samples,
                    gut contents, and other sequences from anonymous
                    organisms that may be associated with a particular
                    host. They do not include endosymbionts that can be
                    reliably recovered from a particular host, organisms
                    from a readily identifiable but uncultured field sample
                    (e.g., many cyanobacteria), or phytoplasmas that can be 
                    reliably recovered from diseased plants (even though 
                    these cannot be grown in axenic culture)."""
    
    value_format = None
    
    comment = """used only with the source feature key; source feature 
                 keys containing the /environmental_sample qualifier 
                 should also contain the /isolation_source qualifier.
                 entries including /environmental_sample must not include 
                 the /strain qualifier."""
    
    def __init__(self, value):
        super(Environmental_sampleQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("environmental_sample")

class Estimated_lengthQualifier( Qualifier ):
    
    definition = "estimated length of the gap in the sequence"
    
    value_format = "unknown or <integer>"
    
    def __init__(self, value):
        super(Estimated_lengthQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("estimated_length", self.value)

class ExceptionQualifier( Qualifier ):
    
    definition = """indicates that the coding region cannot be translated using
                    standard biological rules"""
    
    value_format = """"RNA editing", "reasons given in citation",
                "rearrangement required for product", "annotated by transcript
                or proteomic data\""""
    
    comment = """only to be used to describe biological mechanisms such 
                 as RNA editing;  where the exception cannot easily be described 
                 a published citation must be referred to; protein translation of
                 /exception CDS will be different from the according conceptual 
                 translation; 
                 - An /inference qualifier should accompany any use of
                 /exception="annotated by transcript or proteomic data", to
                 provide support for the existence of the transcript/protein.
                 - must not be used where transl_except would be adequate,
                   e.g. in case of stop codon completion use:
                 /transl_except=(pos:6883,aa:TERM)
                 /note="TAA stop codon is completed by addition of 3' A residues to   
                 mRNA".
                 - must not be used for ribosomal slippage, instead use join operator, 
                   e.g.: CDS   join(486..1784,1787..4810)
                               /note="ribosomal slip on tttt sequence at 1784..1787\""""
    
    def __init__(self, value):
        super(ExceptionQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("exception", self.value)

class ExperimentQualifier( Qualifier ):
    
    definition = """a brief description of the nature of the experimental 
                    evidence that supports the feature identification or assignment."""
    
    value_format = """"[CATEGORY:]text"
                       where CATEGORY is one of the following:
                       "COORDINATES" support for the annotated coordinates
                       "DESCRIPTION" support for a broad concept of function such as that
                       based on phenotype, genetic approach, biochemical function, pathway
                       information, etc.
                       "EXISTENCE" support for the known or inferred existence of the product
                       where text is free text (see examples)"""
    
    comment = """detailed experimental details should not be included, and would
                 normally be found in the cited publications; PMID, DOI and any 
                 experimental database ID is allowed to be used in /experiment
                 qualifier; value "experimental evidence, no additional details
                 recorded" was used to  replace instances of /evidence=EXPERIMENTAL in
                 December 2005"""
    
    def __init__(self, value):
        super(ExperimentQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("experiment", self.value)

class FocusQualifier( Qualifier ):
    
    definition = """identifies the source feature of primary biological
                    interest for records that have multiple source features
                    originating from different organisms and that are not
                    transgenic."""
    
    value_format = None
    
    comment = """the source feature carrying the /focus qualifier
                 identifies the main organism of the entry, this
                 determines: a) the name displayed in the organism
                 lines, b) if no translation table is specified, the
                 translation table, c) the DDBJ/EMBL/GenBank taxonomic
                 division in which the entry will appear; only one
                 source feature with /focus is allowed in an entry; the
                 /focus and /transgenic qualifiers are mutually exclusive
                 in an entry."""
    
    def __init__(self, value):
        super(FocusQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("focus")

class FrequencyQualifier( Qualifier ):
    
    definition = "frequency of the occurrence of a feature"
    
    value_format = """text representing the proportion of a population carrying the
                      feature expressed as a fraction"""
        
    def __init__(self, value):
        super(FrequencyQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("frequency", self.value)

class FunctionQualifier( Qualifier ):
    
    definition = "function attributed to a sequence"
    
    value_format = "text"
    
    comment = """/function is used when the gene name and/or product name do not 
                 convey the function attributable to a sequence."""
    
    def __init__(self, value):
        super(FunctionQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("function", self.value)

class Gap_typeQualifier( Qualifier ):
    
    definition = """type of gap connecting components in records of a genome assembly, 
                    or the type of biological gap in a record that is part of a genome 
                    assembly;"""
    
    value_format = """"between scaffolds", "within scaffold", "telomere", "centromere",
                      "short arm", "heterochromatin", "repeat within scaffold", 
                      "repeat between scaffolds", "unknown\""""
    
    comment = """This qualifier is used only for assembly_gap features and its values
                 are controlled by the AGP Specification version 2.0:
                 https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
                 Please also visit: http://www.insdc.org/controlled-vocabulary-gaptype-qualifier"""
    
    def __init__(self, value):
        super(Gap_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("gap_type", self.value)

class GeneQualifier( Qualifier ):
    
    definition = "symbol of the gene corresponding to a sequence region"
    
    value_format = "text"
    
    def __init__(self, value):
        super(GeneQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("gene", self.value)

class Gene_synonymQualifier( Qualifier ):
    
    definition = "synonymous, replaced, obsolete or former gene symbol"
    
    value_format = "text"
    
    comment = """used where it is helpful to indicate a gene symbol
                 synonym; when used, a primary gene symbol must always be
                 indicated in /gene or a /locus_tag must be used."""
    
    def __init__(self, value):
        super(Gene_synonymQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("gene_synonym", self.value)

class GermlineQualifier( Qualifier ):
    
    definition = """the sequence presented in the entry has not undergone somatic
                    rearrangement as part of an adaptive immune response; it is the
                    unrearranged sequence that was inherited from the parental
                    germline"""
    
    value_format = None
    
    comment = """/germline should not be used to indicate that the source of
                 the sequence is a gamete or germ cell;
                 /germline and /rearranged cannot be used in the same source
                 feature;
                 /germline and /rearranged should only be used for molecules that
                 can undergo somatic rearrangements as part of an adaptive immune 
                 response; these are the T-cell receptor (TCR) and immunoglobulin
                 loci in the jawed vertebrates, and the unrelated variable 
                 lymphocyte receptor (VLR) locus in the jawless fish (lampreys
                 and hagfish);
                 /germline and /rearranged should not be used outside of the
                 Craniata (taxid=89593)"""
    
    def __init__(self, value = None):
        super(GermlineQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("germline")

class HaplogroupQualifier( Qualifier ):
    
    definition = "synonymous, replaced, obsolete or former gene symbol"
    
    value_format = "text"
    
    def __init__(self, value):
        super(HaplogroupQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("haplogroup", self.value)

class HaplotypeQualifier( Qualifier ):
    
    definition = """name for a combination of alleles that are linked together
                    on the same physical chromosome. In the absence of
                    recombination, each haplotype is inherited as a unit, and may
                    be used to track gene flow in populations."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(HaplotypeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("haplotype", self.value)

class HostQualifier( Qualifier ):
    
    definition = """natural (as opposed to laboratory) host to the organism from
                which sequenced molecule was obtained."""
    
    value_format = "text"
    
    def __init__(self, value):
        super(HostQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("host", self.value)

class Identified_byQualifier( Qualifier ):
    
    definition = "name of the expert who identified the specimen taxonomically"
    
    value_format = "text"
    
    def __init__(self, value):
        super(Identified_byQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("identified_by", self.value)

class InferenceQualifier( Qualifier ):
    
    definition = """a structured description of non-experimental evidence that supports
                    the feature identification or assignment."""
    
    value_format = """"[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
  
                      where CATEGORY is one of the following:
                      "COORDINATES" support for the annotated coordinates
                      "DESCRIPTION" support for a broad concept of function such as that
                      based on phenotype, genetic approach, biochemical function, pathway
                      information, etc.
                      "EXISTENCE" support for the known or inferred existence of the product
  
                      where TYPE is one of the following:
                      "non-experimental evidence, no additional details recorded"
                         "similar to sequence"
                            "similar to AA sequence"
                            "similar to DNA sequence"
                            "similar to RNA sequence"
                            "similar to RNA sequence, mRNA"
                            "similar to RNA sequence, EST"
                            "similar to RNA sequence, other RNA"
                         "profile"
                            "nucleotide motif"
                            "protein motif"
                            "ab initio prediction"
                         "alignment"
  
                      where the optional text "(same species)" is included when the
                      inference comes from the same species as the entry.
  
                      where the optional "EVIDENCE_BASIS" is either a reference to a
                      database entry (including accession and version) or an algorithm
                      (including version) , eg 'INSD:AACN010222672.1', 'InterPro:IPR001900',
                      'ProDom:PD000600', 'Genscan:2.0', etc. and is structured 
                      "[ALGORITHM][:EVIDENCE_DBREF[,EVIDENCE_DBREF]*[,...]]\""""
    
    comment = """/inference="non-experimental evidence, no additional details 
                 recorded" was used to replace instances of 
                 /evidence=NOT_EXPERIMENTAL in December 2005; any database ID can be
                 used in /inference= qualifier; recommentations for choice of resource 
                 acronym for[EVIDENCE_BASIS] are provided in the /inference qualifier
                 vocabulary recommendation document (http://www.insdc.org/inference.html);"""
    
    def __init__(self, value):
        super(InferenceQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("inference", self.value)

class IsolateQualifier( Qualifier ):
    
    definition = "individual isolate from which the sequence was obtained"
    
    value_format = "text"
    
    def __init__(self, value):
        super(IsolateQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("isolate", self.value)

class Isolation_sourceQualifier( Qualifier ):
    
    definition = """describes the physical, environmental and/or local
                    geographical source of the biological sample from which
                    the sequence was derived"""
    
    value_format = "text"
    
    comment = """used only with the source feature key;
                 source feature keys containing an /environmental_sample
                 qualifier should also contain an /isolation_source
                 qualifier; the /country qualifier should be used to 
                 describe the country and major geographical sub-region."""
    
    def __init__(self, value):
        super(Isolation_sourceQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("isolation_source", self.value)

class Lab_hostQualifier( Qualifier ):
    
    definition = """scientific name of the laboratory host used to propagate the
                    source organism from which the sequenced molecule was obtained"""
    
    value_format = "text"
    
    comment = """the full binomial scientific name of the host organism should
                 be used when known; extra conditional information relating to
                 the host may also be included."""
    
    def __init__(self, value):
        super(Lab_hostQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("lab_host", self.value)

class Lat_lonQualifier( Qualifier ):
    
    definition = """geographical coordinates of the location where the specimen was
                    collected."""
    
    value_format = "text"
    
    example = """/lat_lon="47.94 N 28.12 W" 
                 /lat_lon="45.0123 S 4.1234 E\""""
    
    comment = """degrees latitude and longitude in format "d[d.dddd] N|S d[dd.dddd] W|E"
                 (see the examples)"""
    
    def __init__(self, value):
        super(Lat_lonQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("lat_lon", self.value)

class Linkage_evidenceQualifier( Qualifier ):
    
    definition = """type of evidence establishing linkage across an 
                    assembly_gap. Only allowed to be used with assembly_gap features that 
                    have a /gap_type value of "within scaffold"or "repeat within scaffold";"""
    
    value_format = """"pcr", "paired-ends", "align genus", "align xgenus", "align trnscpt", "within clone", 
                      "clone contig", "map", "strobe", "unspecified\""""
    
    comment = """This qualifier is used only for assembly_gap features and its values are
                 controlled by the AGP Specification version 2.0:
                 https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
                 Please also visit: http://www.insdc.org/controlled-vocabulary-linkageevidence-qualifier"""
    
    def __init__(self, value):
        super(Linkage_evidenceQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("linkage_evidence", self.value)

class Locus_tagQualifier( Qualifier ):
    
    definition = """a submitter-supplied, systematic, stable identifier for a gene
                    and its associated features, used for tracking purposes"""
    
    value_format = """"text"(single token) but not "<1-5 letters><5-9 digit integer>[.<integer>]\""""
    
    comment = """/locus_tag can be used with any feature that /gene can be used with;  
                identical /locus_tag values may be used within an entry/record, 
                but only if the identical /locus_tag values are associated 
                with the same gene; in all other circumstances the /locus_tag 
                value must be unique within that entry/record. Multiple /locus_tag 
                values are not allowed within one feature for entries created 
                after 15-OCT-2004. 
                If a /locus_tag needs to be re-assigned the /old_locus_tag qualifier 
                should be used to store the old value. The /locus_tag value should
                not be in a format which resembles INSD accession numbers,                 
                accession.version, or /protein_id identifiers."""
    
    def __init__(self, value):
        super(Locus_tagQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("locus_tag", self.value)

class MacronuclearQualifier( Qualifier ):
    
    definition = """if the sequence shown is DNA and from an organism which 
                    undergoes chromosomal differentiation between macronuclear and
                    micronuclear stages, this qualifier is used to denote that the 
                    sequence is from macronuclear DNA. """
    
    value_format = None
    
    def __init__(self, value = None):
        super(MacronuclearQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("macronuclear")

class MapQualifier( Qualifier ):
    
    definition = "genomic map position of feature"
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(MapQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("map", self.value)

class Mating_typeQualifier( Qualifier ):
    
    definition = """mating type of the organism from which the sequence was
                    obtained; mating type is used for prokaryotes, and for
                    eukaryotes that undergo meiosis without sexually dimorphic
                    gametes"""
    
    value_format = "text"
    
    comment = """/mating_type="male" and /mating_type="female" are
                 valid in the prokaryotes, but not in the eukaryotes;
                 for more information, see the entry for /sex."""
    
    def __init__(self, value = None):
        super(Mating_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("mating_type", self.value)

class Mobile_element_typeQualifier( Qualifier ):
    
    definition = """type and name or identifier of the mobile element which is
                    described by the parent feature"""
    
    value_format = """"<mobile_element_type>[:<mobile_element_name>]" where
                mobile_element_type is one of the following:
                "transposon", "retrotransposon", "integron", 
                "insertion sequence", "non-LTR retrotransposon", 
                "SINE", "MITE", "LINE", "other"."""
    
    comment = """/mobile_element_type is legal on mobile_element feature key only.  
                 Mobile element should be used to represent both elements which 
                 are currently mobile, and those which were mobile in the past.  
                 Value "other" requires a mobile_element_name."""
    
    def __init__(self, value = None):
        super(Mobile_element_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("mobile_element_type", self.value)

class Mod_baseQualifier( Qualifier ):
    
    definition = "abbreviation for a modified nucleotide base"
    
    value_format = "modified_base"
    
    comment = """ modified nucleotides not found in the restricted vocabulary
                 list can be annotated by entering '/mod_base=OTHER' with
                 '/note="name of modified base"'"""
    
    def __init__(self, value = None):
        super(Mod_baseQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("mod_base", self.value)

class Mol_typeQualifier( Qualifier ):
    
    definition = "in vivo molecule type of sequence"
    
    value_format = """"genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", "other
                      RNA", "other DNA", "transcribed RNA", "viral cRNA", "unassigned
                      DNA", "unassigned RNA\""""
    
    comment = """all values refer to the in vivo or synthetic molecule for
                 primary entries and the hypothetical molecule in Third Party
                 Annotation entries; the value "genomic DNA" does not imply that
                 the molecule is nuclear (e.g. organelle and plasmid DNA should
                 be described using "genomic DNA"); ribosomal RNA genes should be
                 described using "genomic DNA"; "rRNA" should only be used if the
                 ribosomal RNA molecule itself has been sequenced; /mol_type is
                 mandatory on every source feature key; all /mol_type values
                 within one entry/record must be the same; values "other RNA" and
                 "other DNA" should be applied to synthetic molecules, values
                 "unassigned DNA", "unassigned RNA" should be applied where in
                 vivo molecule is unknown
                 Please also visit:
                 http://www.insdc.org/controlled-vocabulary-moltype-qualifier"""
    
    def __init__(self, value = None):
        super(Mol_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("mol_type", self.value)

class NcRNA_classQualifier( Qualifier ):
    
    definition = """a structured description of the classification of the
                    non-coding RNA described by the ncRNA parent key"""
    
    value_format = "TYPE"
    
    comment = """TYPE is a term taken from the INSDC controlled vocabulary for ncRNA
                 classes (http://www.insdc.org/rna_vocab.html); on
                 15-Oct-2013, the following terms were valid:

                       "antisense_RNA"
                       "autocatalytically_spliced_intron" 
                       "ribozyme"
                       "hammerhead_ribozyme"
                       "lncRNA"
                       "RNase_P_RNA"
                       "RNase_MRP_RNA"
                       "telomerase_RNA"
                       "guide_RNA"
                       "rasiRNA"
                       "scRNA"
                       "siRNA"
                       "miRNA"
                       "piRNA"
                       "snoRNA"
                       "snRNA"
                       "SRP_RNA"
                       "vault_RNA"
                       "Y_RNA"
                       "other"

                 ncRNA classes not yet in the INSDC /ncRNA_class controlled
                 vocabulary can be annotated by entering
                 '/ncRNA_class="other"' with '/note="[brief explanation of
                 novel ncRNA_class]"';"""
    
    def __init__(self, value = None):
        super(NcRNA_classQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("nc", self.value)

class NoteQualifier( Qualifier ):
    
    definition = "any comment or additional information"
    
    value_format = "text"
    
    def __init__(self, value = None):
        
        # if the value is a list, join it by commas
        # TODO: evaluate if this is ultra-dumb.
        if type(value) == type([]):
            value = ", ".join(value)
        
        super(NoteQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("note", self.value)

class NumberQualifier( Qualifier ):
    
    definition = """a number to indicate the order of genetic elements (e.g.,
                    exons or introns) in the 5' to 3' direction"""
    
    value_format = "unquoted text (single token)"
    
    comment = """text limited to integers, letters or combination of integers and/or 
                 letters represented as an unquoted single token (e.g. 5a, XIIb);
                 any additional terms should be included in /standard_name.
                 Example:  /number=2A
                           /standard_name="long\""""
    
    def __init__(self, value = None):
        super(NumberQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("number", self.value)
    
    def add(self, value):
        raise Exception("A feature can only have one number qualifier.")

class Old_locus_tagQualifier( Qualifier ):
    
    definition = "feature tag assigned for tracking purposes "
    
    value_format = """"text" (single token)"""
    
    comment = """/old_locus_tag can be used with any feature where /gene is valid and 
                 where a /locus_tag qualifier is present.  
                 Identical /old_locus_tag values may be used within an entry/record, 
                 but only if the identical /old_locus_tag values are associated 
                 with the same gene; in all other circumstances the /old_locus_tag 
                 value must be unique within that entry/record. 
                 Multiple/old_locus_tag qualifiers with distinct values are 
                 allowed within a single feature; /old_locus_tag and /locus_tag 
                 values must not be identical within a single feature."""
    
    def __init__(self, value = None):
        super(Old_locus_tagQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("old_locus_tag", self.value)
    
    def add(self, value):
        raise Exception("A feature can only have one 'old locus tag' qualifier.")

class OperonQualifier( Qualifier ):
    
    definition = """name of the group of contiguous genes transcribed into a 
                    single transcript to which that feature belongs."""
    
    value_format = "text"
    
    comment = """currently valid only on Prokaryota-specific features"""
    
    def __init__(self, value = None):
        super(OperonQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("operon", self.value)

class OrganelleQualifier( Qualifier ):
    
    definition = """type of membrane-bound intracellular structure from which the 
                    sequence was obtained."""
    
    value_format = """chromatophore, hydrogenosome, mitochondrion, nucleomorph, plastid,
                      mitochondrion:kinetoplast, plastid:chloroplast, plastid:apicoplast,
                      plastid:chromoplast, plastid:cyanelle, plastid:leucoplast, plastid:proplastid"""
    
    comment = """modifier text limited to values from controlled list
                 Please also visit: http://www.insdc.org/controlled-vocabulary-organelle-qualifier"""
    
    def __init__(self, value = None):
        super(OrganelleQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("organelle", self.value)

class OrganismQualifier( Qualifier ):
    
    definition = """scientific name of the organism that provided the 
                    sequenced genetic material."""
    
    value_format = "text"
    
    comment = """the organism name which appears on the OS or ORGANISM line 
                 will match the value of the /organism qualifier of the 
                 source key in the simplest case of a one-source sequence."""
    
    def __init__(self, value = None):
        super(OrganismQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("organism", self.value)

class PartialQualifier( Qualifier ):
    
    definition = """differentiates between complete regions and partial ones."""
    
    value_format = None
    
    comment = """not to be used for new entries from 15-DEC-2001;
                 use '<' and '>' signs in the location descriptors to
                 indicate that the sequence is partial."""
    
    def __init__(self, value = None):
        super(PartialQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("partial")

class PCR_conditionsQualifier( Qualifier ):
    
    definition = """description of reaction conditions and components for PCR."""
    
    value_format = "text"
    
    comment = """used with primer_bind key."""
    
    def __init__(self, value = None):
        super(PCR_conditionsQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("PCR_conditions=\"%s\"" % self.value)

class PCR_primersQualifier( Qualifier ):
    
    definition = """PCR primers that were used to amplify the sequence.
                    A single /PCR_primers qualifier should contain all the primers used  
                    for a single PCR reaction. If multiple forward or reverse primers are                   
                    present in a  single PCR reaction, multiple sets of fwd_name/fwd_seq 
                    or rev_name/rev_seq values will be  present."""
    
    value_format = """/PCR_primers="[fwd_name: XXX1, ]fwd_seq: xxxxx1,[fwd_name: XXX2,]
                      fwd_seq: xxxxx2, [rev_name: YYY1, ]rev_seq: yyyyy1, 
                      [rev_name: YYY2, ]rev_seq: yyyyy2\""""
    
    comment = """fwd_seq and rev_seq are both mandatory; fwd_name and rev_name are
                 both optional. Both sequences should be presented in 5'>3' order. 
                 The sequences should be given in the IUPAC degenerate-base alphabet,
                 except for the modified bases; those must be enclosed within angle
                 brackets <>"""
    
    def __init__(self, value = None):
        super(PCR_conditionsQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("PCR_conditions=\"%s\"" % self.value)

class PhenotypeQualifier( Qualifier ):
    
    definition = """phenotype conferred by the feature, where phenotype is defined as a 
                    physical, biochemical or behavioural characteristic or set of 
                    characteristics."""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(PhenotypeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("phenotype", self.value)

class PlasmidQualifier( Qualifier ):
    
    definition = """name of naturally occurring plasmid from which the sequence was 
                    obtained, where plasmid is defined as an independently replicating
                    genetic unit that cannot be described by /chromosome or /segment"""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(PlasmidQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("Plasmid=\"%s\"" % self.value)

class Pop_variantQualifier( Qualifier ):
    
    definition = """name of subpopulation or phenotype of the sample from which the sequence
                    was derived."""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(Pop_variantQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("pop_variant", self.value)

class ProductQualifier( Qualifier ):
    
    definition = """name of the product associated with the feature, e.g. the mRNA of an 
                    mRNA feature, the polypeptide of a CDS, the mature peptide of a 
                    mat_peptide, etc."""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(ProductQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("product", self.value)

class Protein_idQualifier( Qualifier ):
    
    definition = """protein identifier, issued by International collaborators.
                    this qualifier consists of a stable ID portion (3+5 format
                    with 3 position letters and 5 numbers) plus a version number
                    after the decimal point."""
    
    value_format = "<identifier>"
    
    comment = """when the protein sequence encoded by the CDS changes, only 
                 the version number of the /protein_id value is incremented; 
                 the stable part of the /protein_id remains unchanged and as a
                 result will permanently be associated with a given protein;
                 this qualifier is valid only on CDS features which translate
                 into a valid protein."""
    
    def __init__(self, value = None):
        super(Protein_idQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("protein_id", self.value)

class ProviralQualifier( Qualifier ):
    
    definition = """this qualifier is used to flag sequence obtained from a virus or
                    phage that is integrated into the genome of another organism"""
    
    value_format = None
    
    def __init__(self, value = None):
        super(ProviralQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("proviral")

class PseudoQualifier( Qualifier ):
    
    definition = """indicates that this feature is a non-functional version of the
                    element named by the feature key."""
    
    value_format = None
    
    comment = """The qualifier /pseudo should be used to describe non-functional 
                 genes that are not formally described as pseudogenes, e.g. CDS 
                 has no translation due to other reasons than pseudogenisation events.
                 Other reasons may include sequencing or assembly errors.
                 In order to annotate pseudogenes the qualifier /pseudogene= must be
                 used indicating the TYPE which can be taken from the INSDC controlled vocabulary 
                 for pseudogenes."""
    
    def __init__(self, value = None):
        super(PseudoQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("pseudo")

class PseudogeneQualifier( Qualifier ):
    
    definition = """indicates that this feature is a pseudogene of the element named
                    by the feature key."""
    
    value_format = """"TYPE"
                      where TYPE is one of the following:
                      processed, unprocessed, unitary, allelic, unknown"""
    
    comment = """TYPE is a term taken from the INSDC controlled vocabulary for pseudogenes
                 (http://www.insdc.org/documents/pseudogene-qualifier-vocabulary):

                 processed: the pseudogene has arisen by reverse transcription of a 
                 mRNA into cDNA, followed by reintegration into the genome. Therefore,
                 it has lost any intron/exon structure, and it might have a pseudo-polyA-tail.

                 unprocessed: the pseudogene has arisen from a copy of the parent gene by duplication
                 followed by accumulation of random mutations. The changes, compared to their
                 functional homolog, include insertions, deletions, premature stop codons, frameshifts
                 and a higher proportion of non-synonymous versus synonymous substitutions.

                 unitary: the pseudogene has no parent. It is the original gene, which is
                 functional is some species but disrupted in some way (indels, mutation, 
                 recombination) in another species or strain.

                 allelic: a (unitary) pseudogene that is stable in the population but
                 importantly it has a functional alternative allele also in the population. i.e.,
                 one strain may have the gene, another strain may have the pseudogene.
                 MHC haplotypes have allelic pseudogenes.

                 unknown: the submitter does not know the method of pseudogenisation."""
    
    def __init__(self, value = None):
        super(PseudogeneQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("pseudogene")

class RearrangedQualifier( Qualifier ):
    
    definition = """the sequence presented in the entry has undergone somatic
                    rearrangement as part of an adaptive immune response; it is not
                    the unrearranged sequence that was inherited from the parental
                    germline."""
    
    value_format = None
    
    comment = """/rearranged should not be used to annotate chromosome
                 rearrangements that are not involved in an adaptive immune
                 response;
                 /germline and /rearranged cannot be used in the same source
                 feature;
                 /germline and /rearranged should only be used for molecules that
                 can undergo somatic rearrangements as part of an adaptive immune 
                 response; these are the T-cell receptor (TCR) and immunoglobulin
                 loci in the jawed vertebrates, and the unrelated variable 
                 lymphocyte receptor (VLR) locus in the jawless fish (lampreys
                 and hagfish);
                 /germline and /rearranged should not be used outside of the
                 Craniata (taxid=89593)"""
    
    def __init__(self, value = None):
        super(RearrangedQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("rearranged")

class Regulatory_classQualifier( Qualifier ):
    
    definition = """a structured description of the classification of transcriptional 
                    and translational regulatory elements in a sequence."""
    
    value_format = "TYPE"
    
    comment = """TYPE is a term taken from the INSDC controlled vocabulary for regulatory classes
                 (http://www.insdc.org/controlled-vocabulary-regulatoryclass); on  15-DEC-2014,
                 the following terms were valid:

                 "attenuator"  (replaces attenuator Feature Key)
                 "CAAT_signal"  (replaces CAAT_signal Feature Key)
                 "enhancer"  (replaces enhancer Feature Key)
                 "enhancer_blocking_element"
                 "GC_signal"  (replaces GC_signal Feature Key)
                 "imprinting_control_region"
                 "insulator"
                 "locus_control_region"
                 "minus_35_signal" (replaces -35_signal Feature Key)
                 "minus_10_signal"  (replaces -10_signal Feature Key)
                 "response_element"
                 "polyA_signal_sequence"  (replaces polyA_signal Feature Key)
                 "promoter" (replaces promoter Feature Key)
                 "ribosome_binding_site" (replaces RBS Feature Key)
                 "riboswitch"
                 "silencer"
                 "TATA_box"  (replaces TATA_signal Feature Key)
                 "terminator" (replaces terminator Feature Key)
                 "other"
                 
                 regulatory classes not yet in the INSDC /regulatory_class controlled vocabulary
                 can be annotated by entering /regulatory_class="other" with
                 /note="[brief explanation of novel regulatory_class]";"""
    
    def __init__(self, value = None):
        super(Regulatory_classQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("regulatory_class", self.value)

class ReplaceQualifier( Qualifier ):
    
    definition = """indicates that the sequence identified a feature's intervals is  
                    replaced by the sequence shown in "text"; if no sequence is 
                    contained within the qualifier, this indicates a deletion."""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(ReplaceQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("replace", self.value)

class Ribosomal_slippageQualifier( Qualifier ):
    
    definition = """during protein translation, certain sequences can program
                    ribosomes to change to an alternative reading frame by a 
                    mechanism known as ribosomal slippage."""
    
    value_format = None
    
    comment = """a join operator,e.g.: [join(486..1784,1787..4810)] should be used 
                 in the CDS spans to indicate the location of ribosomal_slippage"""
    
    def __init__(self, value = None):
        super(Ribosomal_slippageQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("ribosomal_slippage")

class Rpt_familyQualifier( Qualifier ):
    
    definition = """type of repeated sequence; "Alu" or "Kpn", for example."""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(Rpt_familyQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("rpt_family", self.value)

class Rpt_typeQualifier( Qualifier ):
    
    definition = """structure and distribution of repeated sequence."""
    
    value_format = """tandem, direct, inverted, flanking, nested, dispersed, terminal, 
                      long_terminal_repeat, non_ltr_retrotransposon_polymeric_tract, 
                      centromeric_repeat, telomeric_repeat, x_element_combinatorial_repeat,
                      y_prime_element and other"""
    
    comment = """the values are case-insensitive, i.e. both "INVERTED" and "inverted" 
                 are valid; For the most current list of allowed values and their definitions please visit:
                 http://www.insdc.org/controlled-vocabulary-rpttype-qualifier"""
    
    def __init__(self, value = None):
        super(Rpt_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("rpt_type", self.value)

class Rpt_unit_rangeQualifier( Qualifier ):
    
    definition = """identity of a repeat range."""
    
    value_format = """<base_range>"""
    
    comment = """used to indicate the base range of the sequence that constitutes 
                 a repeated sequence specified by the feature keys oriT and
                 repeat_region; qualifiers /rpt_unit_range and /rpt_unit_seq
                 replaced qualifier /rpt_unit in December 2005"""
    
    def __init__(self, value = None):
        super(Rpt_unit_rangeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("rpt_unit_range", self.value)

class Rpt_unit_seqQualifier( Qualifier ):
    
    definition = """identity of a repeat sequence."""
    
    value_format = "text"
    
    comment = """used to indicate the literal sequence that constitutes a
                 repeated sequence specified by the feature keys oriT and
                 repeat_region; qualifiers /rpt_unit_range and /rpt_unit_seq
                 replaced qualifier /rpt_unit in December 2005"""
    
    def __init__(self, value = None):
        super(Rpt_unit_seqQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("rpt_unit_seq", self.value)

class SatelliteQualifier( Qualifier ):
    
    definition = """identifier for a satellite DNA marker, compose of many tandem
                    repeats (identical or related) of a short basic repeated unit;"""
    
    value_format = """"<satellite_type>[:<class>][ <identifier>]"
                      where satellite_type is one of the following 
                          "satellite", "microsatellite", "minisatellite\""""
    
    comment = """many satellites have base composition or other properties
                 that differ from those of the rest of the genome that allows
                 them to be identified.
                 Please also visit: http://www.insdc.org/controlled-vocabulary-satellite-qualifier"""
    
    def __init__(self, value = None):
        super(SatelliteQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("satellite", self.value)

class SegmentQualifier( Qualifier ):
    
    definition = """name of viral or phage segment sequenced"""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(SegmentQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("segment", self.value)

class SerotypeQualifier( Qualifier ):
    
    definition = """serological variety of a species characterized by its
                    antigenic properties"""
    
    value_format = "text"
    
    comment = """used only with the source feature key;
                 the Bacteriological Code recommends the use of the
                 term 'serovar' instead of 'serotype' for the 
                 prokaryotes; see the International Code of Nomenclature
                 of Bacteria (1990 Revision) Appendix 10.B "Infraspecific
                 Terms"."""
    
    def __init__(self, value = None):
        super(SerotypeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("serotype", self.value)

class SerovarQualifier( Qualifier ):
    
    definition = """serological variety of a species (usually a prokaryote)
                    characterized by its antigenic properties"""
    
    value_format = "text"
    
    comment = """used only with the source feature key;
                 the Bacteriological Code recommends the use of the
                 term 'serovar' instead of 'serotype' for prokaryotes;
                 see the International Code of Nomenclature of Bacteria
                 (1990 Revision) Appendix 10.B "Infraspecific Terms"."""
    
    def __init__(self, value = None):
        super(SerovarQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("serovar", self.value)

class SexQualifier( Qualifier ):
    
    definition = """sex of the organism from which the sequence was obtained;
                    sex is used for eukaryotic organisms that undergo meiosis
                    and have sexually dimorphic gametes"""
    
    value_format = "text"
    
    comment = """/sex should be used (instead of /mating_type)
                 in the Metazoa, Embryophyta, Rhodophyta & Phaeophyceae;
                 /mating_type should be used (instead of /sex)
                 in the Bacteria, Archaea & Fungi;
                 neither /sex nor /mating_type should be used
                 in the viruses;
                 outside of the taxa listed above, /mating_type
                 should be used unless the value of the qualifier
                 is taken from the vocabulary given in the examples
                 above."""
    
    def __init__(self, value = None):
        super(SexQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("sex", self.value)

class Specimen_voucherQualifier( Qualifier ):
    
    definition = """identifier for the specimen from which the nucleic acid
                    sequenced was obtained"""
    
    value_format = "text"
    
    comment = """the /specimen_voucher qualifier is intended to annotate a
                 reference to the physical specimen that remains after the
                 sequence has been obtained;
                 if the specimen was destroyed in the process of sequencing,
                 electronic images (e-vouchers) are an adequate substitute for a
                 physical voucher specimen; ideally the specimens will be
                 deposited in a curated museum, herbarium, or frozen tissue
                 collection, but often they will remain in a personal or
                 laboratory collection for some time before they are deposited in
                 a curated collection;
                 there are three forms of specimen_voucher qualifiers; if the
                 text of the qualifier includes one or more colons it is a
                 'structured voucher'; structured vouchers include
                 institution-codes (and optional collection-codes) taken from a
                 controlled vocabulary maintained by the INSDC that denotes the
                 museum or herbarium collection where the specimen resides;
                 Please also visit: http://www.insdc.org/controlled-vocabulary-specimenvoucher-qualifier"""
    
    def __init__(self, value = None):
        super(Specimen_voucherQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("specimen_voucher", self.value)

class Standard_nameQualifier( Qualifier ):
    
    definition = """accepted standard name for this feature"""
    
    value_format = "text"
    
    comment = """use /standard_name to give full gene name, but use /gene to
                 give gene symbol (in the above example /gene="Dt")."""
    
    def __init__(self, value = None):
        super(Standard_nameQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("standard_name", self.value)

class StrainQualifier( Qualifier ):
    
    definition = """strain from which sequence was obtained"""
    
    value_format = "text"
    
    comment = """entries including /strain must not include
                 the /environmental_sample qualifier"""
    
    def __init__(self, value = None):
        super(StrainQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("strain", self.value)

class Sub_cloneQualifier( Qualifier ):
    
    definition = """sub-clone from which sequence was obtained"""
    
    value_format = "text"
    
    comment = """the comments on /clone apply to /sub_clone"""
    
    def __init__(self, value = None):
        super(Sub_cloneQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("sub_clone", self.value)

class Sub_speciesQualifier( Qualifier ):
    
    definition = """name of sub-species of organism from which sequence was
                obtained"""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(Sub_speciesQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("sub_species", self.value)

class Sub_strainQualifier( Qualifier ):
    
    definition = """name or identifier of a genetically or otherwise modified 
                    strain from which sequence was obtained, derived from a 
                    parental strain (which should be annotated in the /strain 
                    qualifier).sub_strain from which sequence was obtained"""
    
    value_format = "text"
    
    comment = """If the parental strain is not given, this should
                 be annotated in the strain qualifier instead of sub_strain.
                 Either:
                 /strain="K-12"
                 /sub_strain="MG1655"
                 or:
                 /strain="MG1655\""""
    
    def __init__(self, value = None):
        super(Sub_strainQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("sub_strain", self.value)

class Tag_peptideQualifier( Qualifier ):
    
    definition = """base location encoding the polypeptide for proteolysis tag of 
                    tmRNA and its termination codon;"""
    
    value_format = "<base_range>"
    
    comment = """it is recommended that the amino acid sequence corresponding
                 to the /tag_peptide be annotated by describing a 5' partial 
                 CDS feature; e.g. CDS    <90..122;"""
    
    def __init__(self, value = None):
        super(Tag_peptideQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("tag_peptide", self.value)

class Tissue_libQualifier( Qualifier ):
    
    definition = """tissue library from which sequence was obtained"""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(Tissue_libQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("tissue_lib", self.value)

class Tissue_typeQualifier( Qualifier ):
    
    definition = """tissue type from which the sequence was obtained"""
    
    value_format = "text"
    
    def __init__(self, value = None):
        super(Tissue_typeQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("tissue_type", self.value)

class TransgenicQualifier( Qualifier ):
    
    definition = """identifies the source feature of the organism which was 
                    the recipient of transgenic DNA."""
    
    value_format = None
    
    comment = """transgenic sequences must have at least two source feature keys; 
                 the source feature key having the /transgenic qualifier must 
                 span the whole sequence; the source feature carrying the 
                 /transgenic qualifier identifies the main organism of the entry, 
                 this determines: a) the name displayed in the organism lines, 
                 b) if no translation table is specified, the translation table;
                 only one source feature with /transgenic is allowed in an entry; 
                 the /focus and /transgenic qualifiers are mutually exclusive in 
                 an entry."""
    
    def __init__(self, value = None):
        super(TransgenicQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("transgenic")

class TranslationQualifier( Qualifier ):
    
    definition = """automatically generated one-letter abbreviated amino acid
                    sequence derived from either the universal genetic code or the
                    table as specified in /transl_table and as determined by an
                    exception in the /transl_except qualifier"""
    
    value_format = """IUPAC one-letter amino acid abbreviation, "X" is to be used
                      for AA exceptions."""
    
    comment = """to be used with CDS feature only; this is a mandatory qualifier 
                 in the CDS feature key except where /pseudogene="TYPE" or /pseudo
                 is shown; see /transl_table for definition and location of genetic
                 code tables. """
    
    def __init__(self, value = None):
        super(TranslationQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("translation", self.value)

class Transl_exceptQualifier( Qualifier ):
    
    definition = """translational exception: single codon the translation of which
                    does not conform to genetic code defined by /organism or 
                    /transl_table."""
    
    value_format = """(pos:location,aa:<amino_acid>) where amino_acid is the
                      amino acid coded by the codon at the base_range position."""
    
    comment = """if the amino acid is not on the restricted vocabulary list use
                 e.g., '/transl_except=(pos:213..215,aa:OTHER)' with
                 '/note="name of unusual amino acid"';
                 for modified amino-acid selenocysteine use three letter code
                 'Sec'  (one letter code 'U' in amino-acid sequence)
                 /transl_except=(pos:1002..1004,aa:Sec);
                 for partial termination codons where TAA stop codon is
                 completed by the addition of 3' A residues to the mRNA
                 either a single base_position or a base_range is used, e.g.
                 if partial stop codon is a single base:
                 /transl_except=(pos:1017,aa:TERM)
                 if partial stop codon consists of two bases:
                 /transl_except=(pos:2000..2001,aa:TERM) with
                 '/note='stop codon completed by the addition of 3' A residues 
                 to the mRNA'."""
    
    def __init__(self, value = None):
        super(Transl_exceptQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("transl_except", self.value)

class Transl_tableQualifier( Qualifier ):
    
    definition = """definition of genetic code table used if other than universal
                    genetic code table. Tables used are described in appendix IV."""
    
    value_format = """<integer; 1=universal table 1;2=non-universal table 2;..."""
    
    comment = """genetic code exceptions outside range of specified tables are
                 reported in /transl_except qualifier."""
    
    def __init__(self, value = None):
        super(Transl_tableQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("transl_table", self.value)

class Trans_splicingQualifier( Qualifier ):
    
    definition = """indicates that exons from two RNA molecules are ligated in
                    intermolecular reaction to form mature RNA."""
    
    value_format = None
    
    comment = """should be used on features such as CDS, mRNA and other features
                 that are produced as a result of a trans-splicing event. This
                 qualifier should be used only when the splice event is indicated in
                 the "join" operator, eg join(complement(69611..69724),139856..140087)."""
    
    def __init__(self, value = None):
        super(Trans_splicingQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("trans_splicing")

class Type_materialQualifier( Qualifier ):
    
    definition = """indicates that the organism from which this sequence was obtained is
                    a nomenclatural type of the species (or subspecies) corresponding with
                    the /organism identified in the sequence entry."""
    
    value_format = """"<type-of-type> of <organism name>"
                      where type-of-type is one of the following:
                      type strain, neotype strain, holotype, paratype, neotype, allotype, hapanotype,
                      syntype, lectotype, paralectotype, isotype, epitype, isosyntype, ex-type,
                      reference strain, type material;"""
    
    comment = """<type-of-type> is taken from a controlled vocabularly, listed above.
                 <organism name> should be listed as the scientific name 
                 (or as a synonym) at the species (or subsopecies) node in the taxonomy database.
                 Usage of /type_material will start in the second half of 2014."""
    
    def __init__(self, value = None):
        super(Type_materialQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("type_material", self.value)

class VarietyQualifier( Qualifier ):
    
    definition = """variety (= varietas, a formal Linnaean rank) of organism 
                    from which sequence was derived."""
    
    value_format = "text"
    
    comment = """use the cultivar qualifier for cultivated plant 
                 varieties, i.e., products of artificial selection;
                 varieties other than plant and fungal variatas should be            
                 annotated via /note, e.g. /note="breed:Cukorova\""""
    
    def __init__(self, value = None):
        super(VarietyQualifier, self).__init__(value)
    
    def __repr__(self):
        return self.format("variety", self.value)
