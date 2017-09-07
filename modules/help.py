
def help(string):
   
   if(string == "a" or string == "accession"):   
      return """EMBL specific
      This option is used to set up the accession number of the AC line (see above) and the first token of the ID line as well as the prefix of the locus_tag qualifier.
      Default value = UNKNOWN

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

      Example:        AC   X56734; S46826;

      Similarly, if an existing entry is split into two or more entries (a rare 
      occurrence), the original accession number list is retained in all the derived
      entries.
      An accession number is dropped from the database only when the data to
      which it was assigned have been completely removed from the database.
      """
   elif(string == "c" or string == "created"):   
      return """EMBL specific
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
   elif(string == "d" or string == "data_class"):   
      return """EMBL specific
      This option is used to set up the 5th token of the ID line.
      Default value = STD

      3.1  Data Class
      The data class of each entry, representing a methodological approach to the
      generation of the data or a type of data, is indicated on the first (ID) line
      of the entry. Here an example of ID line:
      ID   Unknown; SV 1; linear; genomic DNA; STD; FUN; 3422859 BP. 
      Each entry belongs to exactly one data class.
        Class          Definition
        -----------    -----------------------------------------------------------
        CON     Entry constructed from segment entry sequences; if unannotated,
                       annotation may be drawn from segment entries
        PAT            Patent
        EST            Expressed Sequence Tag
        GSS            Genome Survey Sequence
        HTC            High Thoughput CDNA sequencing
        HTG            High Thoughput Genome sequencing
        MGA            Mass Genome Annotation
        WGS            Whole Genome Shotgun
        TSA            Transcriptome Shotgun AssEMBLy
        STS            Sequence Tagged Site
        STD            Standard (all entries not classified as above)
      """
   elif(string == "k" or string == "keyword"):
      return """EMBL specific
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
   elif(string == "l" or string == "classification"):
      return """EMBL specific
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
   elif(string == "m" or string == "molecule_type"):
      return """EMBL specific
      This option is used to set up the Molecule type that is the 4th token of the ID line.
      Default value = genomic DNA

      Note 1 - Molecule type: this represents the type of molecule as stored and can
      be any value from the list of current values for the mandatory mol_type source
      qualifier. This item should be the same as the value in the mol_type
      qualifier(s) in a given entry.
      possible value: ["genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", "other RNA", "other DNA", 
                                  "transcribed RNA", "viral cRNA", "unassigned DNA", "unassigned RNA"]
      """
   elif(string == "p" or string == "project_id"):
      return """EMBL specific
      This option is used to set up the PR line. The value should have been provided to you by EMBL.
      Default value = Unknown

      3.4.3  The PR Line
      The PR (PRoject) line shows the International Nucleotide Sequence Database
      Collaboration (INSDC) Project Identifier that has been assigned to the entry.
      Full details of INSDC Project are available at
      http://www.ebi.ac.uk/ena/about/page.php?page=project_guidelines.
      Example:        PR   Project:17285;
      """
   elif(string == "r" or string == "table"):
      return """EMBL specific
      This option is used to set up the translation table qualifier transl_table of the CDS features.
      Default value = 1
      """
   elif(string == "s" or string == "species"):
      return """EMBL specific
      This option is used to set up the OS line.
      Default value = Genus species (name)

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
   elif(string == "t" or string == "topology"):
      return """EMBL specific
      This option is used to set up the Topology that is the 3th token of the ID line.
      Default value = linear

      The possible choice are 'circular' or 'linear'.
      """
   elif(string == "rc"):
      return """EMBL specific
      3.4.10.2  The RC Line
      The RC (Reference Comment) linetype is an optional linetype which appears if 
      The reference has a comment. The comment is in English and as many RC lines as
      are required to display the comment will appear. They are formatted thus:
      RC   comment
      """
   elif(string == "rx"):
      return """EMBL specific
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
           PUBMED         PUBMED bibliographic database (NLM)
           DOI            Digital Object Identifier (International DOI Foundation)
           AGRICOLA       US National Agriculture Library (NAL) of the US Department
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
   elif(string == "rg"):
      return """EMBL specific
      3.4.10.5  The RG Line
      The RG (Reference Group) lines list the working groups/consortia that 
      produced the record. RG line is mainly used in submission reference 
      blocks, but could also be used in paper reference if the working group is 
      cited as an author in the paper.
      """
   elif(string == "ra" or string == "author"):
      return """EMBL specific
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
   elif(string == "rt"):
      return """EMBL specific
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
   elif(string == "rl"):
      return """EMBL specific
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
   elif(string == "version"):
      return """EMBL specific
      This parameter is used to set up the Sequence version number which is the 2nd token of the ID line.
      Default value = 1
      """
   elif(string == "translate"):
      return """EMBL specific
      This parameter is used to add the translation of the CDS sequence into the <translation> qualifier of the CDS feature.
      Default value = 1
      """
   elif(string == "x" or string == "taxonomy"):   
      return """EMBL specific
      This option is used to set the taxonomic division within ID line (6th token).
      This option is mandatory and do not have any default value.

      3.2  Taxonomic Division
      The entries which constitute the database are grouped into taxonomic divisions,
      the object being to create subsets of the database which reflect areas of
      interest for many users.
      In addition to the division, each entry contains a full taxonomic
      classification of the organism that was the source of the stored sequence,
      from kingdom down to genus and species (see below).
      Each entry belongs to exactly one taxonomic division. The ID line of each entry
      indicates its taxonomic division (6th token), using the three letter codes shown below:

                          Division                 Code
                          -----------------        ----
                          Bacteriophage            PHG
                          Environmental Sample     ENV
                          Fungal                   FUN
                          Human                    HUM
                          Invertebrate             INV
                          Other Mammal             MAM
                          Other Vertebrate         VRT
                          Mus musculus             MUS 
                          Plant                    PLN
                          Prokaryote               PRO
                          Other Rodent             ROD
                          Synthetic                SYN
                          Transgenic               TGN
                          Unclassified             UNC
                          Viral                    VRL
      """
   elif(string == "g" or string == "organelle"):
      return """EMBL specific
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

   elif(string == "o" or string == "output"):
      return """EMBLmyGFF3 tool specific
      This option allows to set the file name where the output will be written.
      By default the output is written to STDOUT, which means it will be displayed within the shell you are executing the script, 
      except if you are using the "> afile" to redirect it to the file <afile>.
      """
   elif(string == "shame"):
      return """EMBLmyGFF3 tool specific
      Suppress the shameless plug. A “shameless plug” is a term often used to refer when someone has included (or “plugged”) 
      some information that helps advance their own selfish interests.

      i.e remove the following lines:
      #######################################################################
      # NBIS 2016 - Sweden                                                  #
      # Authors: Martin Norling, Jacques Dainat                             #
      # Please cite NBIS (www.nbis.se) when using this tool.                #
      #######################################################################

      Yes being silent is useful when you do piping and so on.
      """
   elif(string == "v" or string == "verbose"):
      return """EMBLmyGFF3 tool specific
      This option allows to increase the verbosity level. V stands for verbose.

      There is 5 verbosity level:  Critical > Error > Warning > Info > Debug
      By default the verbosity level is set to Warning (Critical, Error, Warning messages will be displayed).
      -v will increase the level to the Info level
      -vv will increase the level to the Debug level
      """
   elif(string == "q" or string == "quiet"):
      return """EMBLmyGFF3 tool specific
      This option allows to decrease the verbosity level. Q stands for quiet

      There is 5 verbosity level:  Critical > Error > Warning > Info > Debug
      By default the verbosity level is set to Warning (Critical, Error, Warning messages will be displayed).
      -q will decrease the level to the Error level
      -qq will decrease the level to the Critical level
      """
   elif(string == "z" or string == "gzip"):
      return """EMBLmyGFF3 tool specific
      This option allows to compress the output file in gzip format. This option does not expect any value.
      """
   else:
      return """Sorry we don't recognize this option"""+string
