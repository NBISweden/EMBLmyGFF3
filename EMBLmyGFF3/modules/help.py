#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import logging
import sys

def Help(string):
  output=""

  if(string == "a" or string == "accession" or string == "all"):   
    output += string+""":
EMBL specific
This option is used to set up the accession number of the AC line (see above) and the first token of the ID line.
Default value = XXX
Advise: Let the default value as it is, ENA will automatically replace it during the submission by a unique accession number they will assign.

3.4.2  The AC Line
The AC (ACcession number) line lists the accession numbers associated with 
the entry.

Examples of accession number lines are shown below:
 AC   X56734; S46826;
 AC   Y00001; X00001-X00005; X00008; Z00001-Z00005;
Each accession number, or range of accession numbers, is terminated by a
semicolon. Where necessary, more than one AC line is used. Consecutive
secondary accession numbers in ENA flatfiles are shown in the form of 
inclusive accession number ranges.
Accession numbers are the primary means of identifying sequences providing 
a stable way of identifying entries from release to release. An accession
number, however, always remains in the accession number list of the latest
version of the entry in which it first appeared.  Accession numbers allow
unambiguous citation of database entries. Researchers who wish to cite entries
in their publications should always cite the first accession number in the
list (the "primary" accession number) to ensure that readers can find the
relevant data in a subsequent release. Readers wishing to find the data thus
cited must look at all the accession numbers in each entry's list.
Secondary accession numbers: One reason for allowing the existence of several
accession numbers is to allow tracking of data when entries are merged
or split. For example, when two entries are merged into one, a "primary" 
accession number goes at the start of the list, and those from the 
merged entries are added after this one as "secondary" numbers.  

Example:  AC   X56734; S46826;

Similarly, if an existing entry is split into two or more entries (a rare 
occurrence), the original accession number list is retained in all the derived
entries.
An accession number is dropped from the database only when the data to
which it was assigned have been completely removed from the database.
"""
  if(string == "i" or string == "locus_tag"):   
    output += string+""":
This option is used to set up the prefix of the locus_tag qualifier.
Mandatory - Default value = XXX

Locus tags are identifiers that are systematically applied to every gene in a genome within the context of sequencing projects. These tags have become surrogate gene names by the biological community. If two submitters of two different genomes use the same systematic names to describe two very different genes in two very different genomes, it can be very confusing. In order to prevent this from happening INSDC has created a registry of locus tag prefixes. Submitters of eukaryotic and prokaryotic genomes should register a locus tag prefix prior to submitting their genome annotation into ENA.

Locus tags can be registered in Webin at the time of project registration. All genome assembly and annotation projects are required to register a sequencing project. If during the submission process the project is said to contain functional annotation, then the user will be prompted to register a locus tag prefix. Users can opt to select their preferred locus tag prefix (subject to availability) or to have ENA assign one automatically from one of the following ranges: BN1-BN9999, BQ1-BQ9999 or CZ1-CZ9999.

The locus tag prefix is to be separated from the tag value by an underscore ‘_’ (e.g. /locus_tag='BN5_00001').

Locus tags should be assigned to all protein coding and non-coding genes such as structural RNAs. /locus_tag should appear on gene, mRNA, CDS, 5'UTR, 3'UTR, intron, exon, tRNA, rRNA, misc_RNA, etc. within a genome project submission. We discourage the use of the /locus_tag qualifier on repeat_region and misc_feature features in the context of complete genome annotation. The same /locus_tag should be used for all components of a single gene. For example, all of the exons, CDS, mRNA and gene features for a particular gene would have the same /locus_tag. There should only be one /locus_tag associated with one /gene, i.e. if a /locus_tag is associated with a /gene symbol in any feature, that gene symbols (and only that /gene symbol) must also be present on every other feature that contains that /locus_tag.

Locus tags are systematically added to genes within a genome. They are generally in sequential order on the genome. If a genome center were to update a genome and provide additional annotation, the new genes could either be assigned the next sequential available /locus_tag or the submitter can leave gaps when initially assigning /locus_tags and fill in new annotation with tag values that are between the gaps.
For more information please visit: https://www.ebi.ac.uk/ena/submit/locus-tags
"""

  if(string == "c" or string == "created" or string == "all"):   
    output += string+""":
EMBL specific
This option is used to set up the DT line (see above).
The default value will be the current date.

3.4.4  The DT Line
The DT (DaTe) line shows when an entry first appeared in the database and
when it was last updated.  Each entry contains two DT lines, formatted
as follows:
DT   DD-MON-YYYY (Rel. #, Created)
DT   DD-MON-YYYY (Rel. #, Last updated, Version #)
The DT lines from the above example are:
DT   12-SEP-1991 (Rel. 29, Created)
DT   13-SEP-1993 (Rel. 37, Last updated, Version 8)
The date supplied on each DT line indicates when the entry was created or 
Last updated; that will usually also be the date when the new or modified 
Entry became publicly visible via the EBI network servers. The release 
number indicates the first quarterly release made *after* the entry was 
created or last updated. The version number appears only on the "Last 
updated" DT line.
The absolute value of the version number is of no particular significance; its
purpose is to allow users to determine easily if the version of an entry 
which they already have is still the most up to date version. Version numbers
are incremented by one every time an entry is updated; since an entry may be
updated several times before its first appearance in a quarterly release, the
version number at the time of its first release appearance may be greater than
one. Note that because an entry may also be updated several times between
two quarterly releases, there may be gaps in the sequence of version numbers 
which appear in consecutive releases.
If an entry has not been updated since it was created, it will still have 
two DT lines and the "Last updated" line will have the same date (and 
release number) as the "Created" line.
"""
  if(string == "d" or string == "data_class" or string == "all"):   
    output += string+""":
EMBL specific
This option is used to set up the 5th token of the ID line.
Default value = XXX

3.1  Data Class
The data class of each entry, representing a methodological approach to the
generation of the data or a type of data, is indicated on the first (ID) line
of the entry. Here an example of ID line:

ID   XXX; SV 1; linear; genomic DNA; STD; FUN; 3422859 BP.

Each entry belongs to exactly one data class.
  Class    Definition
  -----------    -----------------------------------------------------------
  CON     Entry constructed from segment entry sequences; if unannotated,
     annotation may be drawn from segment entries
  PATPatent
  ESTExpressed Sequence Tag
  GSSGenome Survey Sequence
  HTCHigh Thoughput CDNA sequencing
  HTGHigh Thoughput Genome sequencing
  MGAMass Genome Annotation
  WGSWhole Genome Shotgun
  TSATranscriptome Shotgun AssEMBLy
  STSSequence Tagged Site
  STDStandard (all entries not classified as above)
"""
  if(string == "k" or string == "keyword" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the KW line.
No default value.

3.4.6  The KW Line
The KW (KeyWord) lines provide information which can be used to generate
cross-reference indexes of the sequence entries based on functional,
structural, or other categories deemed important.
The format for a KW line is:
     KW   keyword[; keyword ...].
More than one keyword may be listed on each KW line; the keywords are 
separated by semicolons, and the last keyword is followed by a full
stop. Keywords may consist of more than one word, and they may contain
embedded blanks and stops. A keyword is never split between lines. 
An example of a keyword line is:
     KW   beta-glucosidase.
The keywords are ordered alphabetically; the ordering implies no hierarchy
of importance or function.  If an entry has no keywords assigned to it,
it will contain a single KW line like this:
     KW   .
"""
  if(string == "l" or string == "classification" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the OC line.
Default value = Life

3.4.8  The OC Line
The OC (Organism Classification) lines contain the taxonomic classification
Of the source organism as described in Section 2.2 above. 
The classification is listed top-down as nodes in a taxonomic tree in which 
the most general grouping is given first.  The classification may be 
distributed over several OC lines, but nodes are not split or hyphenated 
between lines. The individual items are separated by semicolons and the
list is terminated by a full stop. The format for the OC line is:
     OC   Node[; Node...].
     
Example classification lines:
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   euphyllophytes; Spermatophyta; Magnoliophyta; eudicotyledons; Rosidae;
OC   Fabales; Fabaceae; Papilionoideae; Trifolium.
"""
  if(string == "m" or string == "molecule_type" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the Molecule type that is the 4th token of the ID line.
Mandatory - No default value

Note 1 - Molecule type: this represents the type of molecule as stored and can
be any value from the list of current values for the mandatory mol_type source
qualifier. This item should be the same as the value in the mol_type
qualifier(s) in a given entry.
possible value: ["genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", "other RNA", "other DNA", 
    "transcribed RNA", "viral cRNA", "unassigned DNA", "unassigned RNA"]
"""
  if(string == "p" or string == "project_id" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the PR line. The value should have been provided to you by EMBL.
Mandatory - Default value = XXX

3.4.3  The PR Line
The PR (PRoject) line shows the International Nucleotide Sequence Database
Collaboration (INSDC) Project Identifier that has been assigned to the entry.
Full details of INSDC Project are available at
http://www.ebi.ac.uk/ena/about/page.php?page=project_guidelines.
Example:  PR   Project:17285;
"""
  if(string == "r" or string == "transl_table" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the translation table qualifier transl_table of the CDS features.
Mandatory and no default
"""
  if(string == "s" or string == "species" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the OS line.
Mandatory - No default value

3.4.7  The OS Line
The OS (Organism Species) line specifies the preferred scientific name of
the organism which was the source of the stored sequence. In most 
cases this is done by giving the Latin genus and species designations, 
followed (in parentheses) by the preferred common name in English where
known. The format is:
     OS   Genus species (name)
In some cases, particularly for viruses and genetic elements, the only
accepted designation is a simple name such as "Canine adenovirus type 2".
In these cases only this designation is given. The species line from the 
example is:
     OS   Trifolium repens (white clover)
Hybrid organisms are classified in their own right. A rat/mouse hybrid,
for example, would appear as follows:
     OS   Mus musculus x Rattus norvegicus
     OC   (OC for mouse)
 
If the source organism is unknown but has been/will be cultured, the OS
line will contain a unique name derived from the what is known of the
classification. The unique name serves to identify the database entry,
which will be updated once the full classification is known. In the
case of an unknown bacterium, for example:
     OS   unidentified bacterium B8
     OC   Bacteria.
For environmental samples where there is no intention to culture the
organism and complete taxonomy cannot be determined, collective names
are used in the OS line and the classification given extends down to
the most resolved taxonomic node possible, for example:
     OS   uncultured proteobacterium
     OC   Bacteria; Proteobacteria; environmental samples.
 
For naturally occurring plasmids the OS/OC lines will contain the 
source organism and the plasmid name will appear on the OG line. 
For example:
     OS   Escherichia coli
     OC   Prokaryota; ... Enterobacteriaceae.
     XX
     OG   Plasmid colE1
For artificial plasmids the OS line will be "OS Cloning vector" and the
sequence will be classified as an artificial sequence. For example:
     OS   Cloning vector M13plex17 
     OC   Artificial sequences; vectors.
 
Where only a naturally occurring part of a plasmid is reported, the plasmid
name will appear on the OG line and the OS/OC lines will describe the natural
source.
For example:
     OS   Escherichia coli
     OC   Prokaryota; ... Enterobacteriaceae.
     XX
     OG   Plasmid pUC8

P.S: The organism name which appears on the OS or ORGANISM line will match the value of
     the /organism qualifier of the source key in the simplest case of a one-source sequence.
"""
  if(string == "t" or string == "topology" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set up the Topology that is the 3th token of the ID line.
Mandatory - No default value

The possible choice are 'circular' or 'linear'.
"""
  if(string == "rc" or string == "all"):
    output += string+""":
EMBL specific
3.4.10.2  The RC Line
The RC (Reference Comment) linetype is an optional linetype which appears if 
The reference has a comment. The comment is in English and as many RC lines as
are required to display the comment will appear. They are formatted thus:
RC   comment
"""
  if(string == "rx" or string == "all"):
    output += string+""":
EMBL specific
3.4.10.4  The RX Line
The RX (reference cross-reference) linetype is an optional linetype which
contains a cross-reference to an external citation or abstract resource.
For example, if a journal citation exists in the PUBMED database, there will
be an RX line pointing to the relevant PUBMED identifier.
The format of the RX line is as follows:
     RX  resource_identifier; identifier.   
The first item on the RX line, the resource identifier, is the abbreviated 
name of the data collection to which reference is made. The current
set of cross-referenced resources is:
     Resource ID    Fullname
     -----------    ------------------------------------
     PUBMED   PUBMED bibliographic database (NLM)
     DOIDigital Object Identifier (International DOI Foundation)
     AGRICOLA US National Agriculture Library (NAL) of the US Department
  of Agriculture (USDA)
The second item on the RX line, the identifier, is a pointer to the entry in
the external resource to which reference is being made. The data item used as
the primary identifier depends on the resource being referenced.
For example:
RX   DOI; 10.1016/0024-3205(83)90010-3.
RX   PUBMED; 264242.
Note that further details of DOI are available at http://www.doi.org/. URLs
formulated in the following way are resolved to the correct full text URLs:
     http://dx.doi.org/<doi>
     eg. http:/dx.doi.org/10.1016/0024-3205(83)90010-3
"""
  if(string == "rg" or string == "all"):
    output += string+""":
EMBL specific
Default value = XXX

3.4.10.5  The RG Line
The RG (Reference Group) lines list the working groups/consortia that 
produced the record. RG line is mainly used in submission reference 
blocks, but could also be used in paper reference if the working group is 
cited as an author in the paper.
"""
  if(string == "ra" or string == "author" or string == "all"):
    output += string+""":
EMBL specific
3.4.10.6  The RA Line
The RA (Reference Author) lines list the authors of the paper (or other 
work) cited. All of the authors are included, and are listed in the order 
given in the paper. The names are listed surname first followed by a blank
followed by initial(s) with stops. Occasionally the initials may not 
be known, in which case the surname alone will be listed. The author names 
are separated by commas and terminated by a semicolon; they are not split 
between lines. The RA line in the example is:
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;    
As many RA lines as necessary are included for each reference.
"""
  if(string == "rt" or string == "all"):
    output += string+""":
EMBL specific
3.4.10.7  The RT Line
The RT (Reference Title) lines give the title of the paper (or other work) as
exactly as is possible given the limitations of computer character sets. Note
that the form used is that which would be used in a citation rather than that
displayed at the top of the published paper. For instance, where journals
capitalise major title words this is not preserved. The title is enclosed in
double quotes, and may be continued over several lines as necessary. The title
lines are terminated by a semicolon. The title lines from the example are:
RT   "Nucleotide and derived amino acid sequence of the cyanogenic
RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.)";
Greek letters in titles are spelled out; for example, a title in an entry 
would contain "kappa-immunoglobulin" even though the letter itself may be
present in the original title. Similar simplifications have been made in 
other cases (e.g. subscripts and superscripts). Note that the RT line of
a citation which has no title (such as a submission to the database) contains
only a semicolon.
"""
  if(string == "rl" or string == "all"):
    output += string+""":
EMBL specific
3.4.10.8  The RL Line
The RL (Reference Location) line contains the conventional citation 
information for the reference.  In general, the RL lines alone are 
sufficient to find the paper in question. They include the journal,
volume number, page range and year for each paper. 
Journal names are abbreviated according to existing ISO standards 
(International Standard Serial Number)
The format for the location lines is:
     RL   journal vol:pp-pp(year).
Thus, the reference location line in the example is:
     RL   Plant Mol. Biol. 17:209-219(1991).
Very occasionally a journal is encountered which does not consecutively 
number pages within a volume, but rather starts the numbering anew for
each issue number. In this case the issue number must be included, and the 
format becomes:
     RL   journal vol(no):pp-pp(year).
 
If a paper is in press, the RL line will appear with such information as 
we have available, the missing items appearing as zeros. For example:
     RL   Nucleic Acids Res. 0:0-0(2004).
This indicates a paper which will be published in Nucleic Acids Research at some
point in 2004, for which we have no volume or page information. Such references
are updated to include the missing information when it becomes available.
Another variation of the RL line is used for papers found in books 
or other similar publications, which are cited as shown below:
     RA   Birnstiel M., Portmann R., Busslinger M., Schaffner W.,
     RA   Probst E., Kressmeann A.;
     RT   "Functional organization of the histone genes in the
     RT   sea urchin Psammechinus:  A progress report";
     RL   (in) Engberg J., Klenow H., Leick V. (Eds.);
     RL   SPECIFIC EUKARYOTIC GENES:117-132;
     RL   Munksgaard, Copenhagen (1979).
Note specifically that the line where one would normally encounter the 
journal location is replaced with lines giving the bibliographic citation
of the book. The first RL line in this case contains the designation "(in)",
which indicates that this is a book reference.
The following examples illustrate RL line formats that are used for data
submissions:
     RL   Submitted (19-NOV-1990) to the INSDC.
     RL   M.A. Hughes, UNIVERSITY OF NEWCASTLE UPON TYNE, MEDICAL SCHOOL, NEW
     RL   CASTLE UPON TYNE, NE2  4HH, UK
Submitter address is always included in new entries, but some older 
submissions do not have this information. 
RL lines take another form for thesis references. 
For example:
     RL   Thesis (1999), Department of Genetics,
     RL   University of Cambridge, Cambridge, U.K.
For an unpublished reference, the RL line takes the following form:
     RL   Unpublished.
Patent references have the following form:
     RL   Patent number EP0238993-A/3, 30-SEP-1987.
     RL   BAYER AG.
The words "Patent number" are followed by the patent application number, the
patent type (separated by a hyphen), the sequence's serial number within the
patent (separated by a slash) and the patent application date. The subsequent RL
lines list the patent applicants, normally company names.
Finally, for journal publications where no ISSN number is available for the
journal (proceedings and abstracts, for example), the RL line contains the
designation "(misc)" as in the following example.
     RL   (misc) Proc. Vth Int. Symp. Biol. Terr. Isopods 2:365-380(2003).
"""
  if(string == "version" or string == "all"):
    output += string+""":
EMBL specific
This parameter is used to set up the Sequence version number which is the 2nd token of the ID line.
Default value = 1
"""
  if(string == "translate" or string == "all"):
    output += string+""":
EMBL specific
This parameter is used to add the translation of the CDS sequence into the <translation> qualifier of the CDS feature.
The option doesn't expect any value, use --translate to activate the translation.
"""
  if(string == "x" or string == "taxonomy" or string == "all"):   
    output += string+""":
EMBL specific
This option is used to set the taxonomic division within ID line (6th token).
Default value = XXX

3.2  Taxonomic Division
The entries which constitute the database are grouped into taxonomic divisions,
the object being to create subsets of the database which reflect areas of
interest for many users.
In addition to the division, each entry contains a full taxonomic
classification of the organism that was the source of the stored sequence,
from kingdom down to genus and species (see below).
Each entry belongs to exactly one taxonomic division. The ID line of each entry
indicates its taxonomic division (6th token), using the three letter codes shown below:

  Division     Code
  -----------------  ----
  BacteriophagePHG
  Environmental Sample     ENV
  Fungal FUN
  Human  HUM
  Invertebrate INV
  Other Mammal MAM
  Other Vertebrate   VRT
  Mus musculus MUS 
  Plant  PLN
  Prokaryote   PRO
  Other Rodent ROD
  Synthetic    SYN
  Transgenic   TGN
  Unclassified UNC
  Viral  VRL
"""
  if(string == "g" or string == "organelle" or string == "all"):
    output += string+""":
EMBL specific
This option is used to set uo the OG line.
There is no default value. Possible choices are: 
"chromatophore", "hydrogenosome", "mitochondrion", "nucleomorph", "plastid", 
"mitochondrion:kinetoplast", "plastid:chloroplast", "plastid:apicoplast", 
"plastid:chromoplast", "plastid:cyanelle", "plastid:leucoplast", "plastid:proplastid"

3.4.9  The OG Line
The OG (OrGanelle) linetype indicates the sub-cellular location of non-nuclear
sequences.  It is only present in entries containing non-nuclear sequences
and appears after the last OC line in such entries.
The OG line contains
a) one data item (title cased) from the controlled list detailed under the
/organelle qualifier definition in the Feature Table Definition document
that accompanies this release or
b) a plasmid name.
Examples include "Mitochondrion", "Plastid:Chloroplast" and "Plasmid pBR322".

For example, a chloroplast sequence from Euglena gracilis would appear as:
     OS   Euglena gracilis (green algae)
     OC   Eukaryota; Planta; Phycophyta; Euglenophyceae.
     OG   Plastid:Chloroplast
"""

  if(string == "o" or string == "output" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
This option allows to set the file name where the output will be written.
By default the output is written to STDOUT, which means it will be displayed within the shell you are executing the script, 
except if you are using the "> afile" to redirect it to the file <afile>.
"""
  if(string == "shame" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Suppress the shameless plug. A "shameless plug" is a term often used to refer when someone has included (or "plugged") 
some information that helps advance their own selfish interests.

i.e remove the following lines:
#######################################################################
# NBIS 2016 - Sweden  #
# Authors: Martin Norling, Jacques Dainat     #
# Please cite NBIS (www.nbis.se) when using this tool.    #
#######################################################################

Yes being silent is useful when you do piping and so on.
"""
  if(string == "v" or string == "verbose" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
This option allows to increase the verbosity level. V stands for verbose.

There is 5 verbosity level:  Critical > Error > Warning > Info > Debug
By default the verbosity level is set to Warning (Critical, Error, Warning messages will be displayed).
-v will increase the level to the Info level
-vv will increase the level to the Debug level
"""
  if(string == "q" or string == "quiet" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
This option allows to decrease the verbosity level. Q stands for quiet

There is 5 verbosity level:  Critical > Error > Warning > Info > Debug
By default the verbosity level is set to Warning (Critical, Error, Warning messages will be displayed).
-q will decrease the level to the Error level
-qq will decrease the level to the Critical level
"""
  if(string == "z" or string == "gzip" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
This option allows to compress the output file in gzip format. This option does not expect any value.
"""
  if(string == "uncompressed_log" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Some of the log could be compressed for a better lisibility, using this option they won't.
"""
  if(string == "email"  or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Email used to fetch information from NCBI taxonomy database. Default value 'EMBLmyGFF3@tool.org'.
/!\ This email could be blocked at any time if someone try to access the database too many time in a short period.
! It is advised to use its own email address.
""" 
  if(string == "interleave_genes" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Print gene features with interleaved mRNA and CDS features.
""" 
  if(string == "keep_duplicates" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Do not remove duplicate features during the process. 
/!\ Option not suitable for submission purpose. Features that have the same key (feature type) and location as another feature are considered as duplicates and aren't allowed by the EMBL database. So they are remove during the process. If you don't plan to submit the file to ENA and you wish to keep these features, use the --keep_duplicates option.
""" 
  if(string == "force_unknown_features" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Force to keep feature types not accepted by EMBL. 
/!\ Option not suitable for submission purpose. Indeed, an EMBL file with a feature type not accepted by ENA will not be valid.
"""
  if(string == "force_uncomplete_features" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Force to keep features whithout all the mandatory qualifiers. 
/!\ Option not suitable for submission purpose. Indeed, an EMBL file with a feature that have a mandatory qualifier missing will not be valid.
""" 
  if(string == "no_wrap_qualifier" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
Do not wrap the qualifier line output. One line by qualifier. 
By default the output is wrapped at 80 characters, and we cut at word levels.
"""
  if(string == "expose_translations" or string == "all"):
    output += string+""":
EMBLmyGFF3 tool specific
Bolean - Doesnt expect any value
When the option is called, the mapping json files will be copied into the current folder. A local mapping json file will always be used instead of the default internal one.
""" 
##Procaryote specific parameters
  if(string == "strain" or string == "all"):
    output += string+""":
EMBL specific - used to create the source feature keys.
At least one of the following qualifiers \"strain, environmental_sample, isolate\" must exist when organism belongs to Bacteria.
Strain from which sequence was obtained. Entries including /environmental_sample must not include the /strain qualifier.
""" 
  if(string == "environmental_sample" or string == "all"):
    output += string+""":
EMBL specific - used to create the source feature keys.
Bolean - Doesnt expect any value
At least one of the following qualifiers \"strain, environmental_sample, isolate\" must exist when organism belongs to Bacteria.
Definition:  identifies sequences derived by direct molecular isolation from a bulk environmental DNA sample (by PCR with or without subsequent cloning of the product, DGGE, or other anonymous methods) with no reliable identification of the source organism. Environmental samples include clinical samples, gut contents, and other sequences from anonymous organisms that may be associated with a particular host. They do not include endosymbionts that can be reliably recovered from a particular host, organisms from a readily identifiable but uncultured field sample (e.g., many cyanobacteria), or phytoplasmas that can be reliably recovered from diseased plants (even though these cannot be grown in axenic culture).
Comment: source feature keys containing the /environmental_sample qualifier should also contain the /isolation_source qualifier. Entries including /environmental_sample must not include the /strain qualifier.
""" 
  if(string == "isolation_source" or string == "all"):
    output += string+""":
EMBL specific - used to create the source feature keys.
Definition:  describes the physical, environmental and/or local geographical source of the biological sample from which the sequence was derived.
Comment: source feature keys containing an /environmental_sample qualifier should also contain an /isolation_source qualifier; the /country qualifier should be used to describe the country and major geographical sub-region.
""" 
  if(string == "isolate" or string == "all"):
    output += string+""":
EMBL specific - used to create the source feature keys.
At least one of the following qualifiers \"strain, environmental_sample, isolate\" must exist when organism belongs to Bacteria.
Individual isolate from which the sequence was obtained.
""" 

  if not output:
    output += string+""":Sorry no advanced help for this option: """+string+"\n"

  return output