#!/usr/bin/env python3
"""
Formatting class for keep track of the values associated with an EMBL file
header, as described in ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.
"""

from datetime import datetime

from .embl_reference import EMBLReference
from .embl_utilities import embl_line, get_ena_release, taxid_to_species, \
                            classification_from_taxid, species_to_taxid

class EMBLHeader():
    """Formatting class for keep track of the values associated with an EMBL
    file header, as described in:
    ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt.
    """

    def __init__(self, args=None):
        """
        Takes the input arguments from the main gffemblconverter, as an
        argparse Namespace object. This class then formats and checks what goes
        into the EMBL file header.
        """

        self.settings = {"accession":None, "version":None, "topology":None,
                         "molecule_type":None, "data_class":None,
                         "taxonomic_division":None, "sequence_length":None,
                         "created":None}

        if args:
            for key, value in vars(args).items():
                self.settings[key] = value

        # format settings
        if self.settings.get('created', False):
            try:
                self.settings['created'] = \
                    datetime.strptime(self.settings['created'], "%Y-%m-%d")
            except ValueError:
                self.settings['created'] = \
                    datetime.strptime(self.settings['created'], "%d-%b-%Y")
            self.settings['created_release'] = \
                get_ena_release(self.settings['created'])

        # Add generated settings
        self.settings["updated"] = datetime.today()
        self.settings["updated_release"] = \
            get_ena_release(self.settings['updated'])

        # Add references
        self.references = [EMBLReference(args)]

    def __repr__(self):
        # ID - identification             (begins each entry; 1 per entry)
        output = self.id_line()
        # AC - accession number           (>=1 per entry)
        output += self.accession_line()
        # Submission AC
        if self.settings.get("record_id", False):
            output += self.accession_line(True, self.settings["record_id"])
        # PR - project identifier         (0 or 1 per entry)
        output += self.project_line()
        # DT - date                       (2 per entry)
        output += self.date_line()
        # DE - description                (>=1 per entry)
        output += self.description()
        # KW - keyword                    (>=1 per entry)
        output += self.keywords()
        # OS - organism species           (>=1 per entry)
        output += self.species()
        # OC - organism classification    (>=1 per entry)
        output += self.classification()
        # OG - organelle                  (0 or 1 per entry)
        output += self.organelle()
        # References, including RN, RC, RP, RX, RG, RA, RT and RL
        for reference in self.references:
            output += str(reference)
        # DR - database cross-reference   (>=0 per entry)
        output += self.xref()
        # CC - comments or notes          (>=0 per entry)
        output += self.comment()

        # output += self.third_party_annotation()

        # output += self.assembly_information()

        return output

    def id_line(self):
        """
        3.4.1  The ID Line
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

        Note 1 - Molecule type: this represents the type of molecule as stored
        and can be any value from the list of current values for the mandatory
        mol_type source qualifier. This item should be the same as the value
        in the mol_type qualifier(s) in a given entry.
        Note 2 - Sequence length: The last item on the ID line is the length
        of the sequence (the total number of bases in the sequence). This
        number includes base positions reported as present but undetermined
        (coded as "N").
        An example of a complete identification line is shown below:
        ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.
        """
        line_code = "ID"
        template = ("{accession}; SV {version}; {topology}; "
                    "{molecule_type}; {data_class}; "
                    "{taxonomy}; {sequence_length} BP.")

        return embl_line(line_code, template.format(**self.settings))

    def accession_line(self, submission=False, accession=False):
        """
        3.4.2  The AC Line
        The AC (ACcession number) line lists the accession numbers associated
        with the entry.

        Examples of accession number lines are shown below:
        AC   X56734; S46826;
        AC   Y00001; X00001-X00005; X00008; Z00001-Z00005;
        Each accession number, or range of accession numbers, is terminated by
        a semicolon. Where necessary, more than one AC line is used.
        Consecutive secondary accession numbers in ENA flatfiles are shown in
        the form of inclusive accession number ranges.
        Accession numbers are the primary means of identifying sequences
        providing a stable way of identifying entries from release to release.
        An accession number, however, always remains in the accession number
        list of the latest version of the entry in which it first appeared.
        Accession numbers allow unambiguous citation of database entries.
        Researchers who wish to cite entries in their publications should
        always cite the first accession number in the list (the "primary"
        accession number) to ensure that readers can find the relevant data in
        a subsequent release. Readers wishing to find the data thus cited must
        look at all the accession numbers in each entry's list. Secondary
        accession numbers: One reason for allowing the existence of several
        accession numbers is to allow tracking of data when entries are merged
        or split. For example, when two entries are merged into one, a
        "primary" accession number goes at the start of the list, and those
        from the merged entries are added after this one as "secondary"
        numbers.

        Example:        AC   X56734; S46826;

        Similarly, if an existing entry is split into two or more entries (a
        rare occurrence), the original accession number list is retained in all
        the derived entries.
        An accession number is dropped from the database only when the data to
        which it was assigned have been completely removed from the database.

        The submission case:
        The entry name is extracted from the AC * line . The entry name must be
        prefixed with a '_' when using the flat file format. No spaces or pipe
        character ('|') are allowed in the name.

        Example: AC * _contig1
        """
        line_code = "AC *" if submission else "AC"
        template = "{accession};"

        if accession:
            if not accession.startswith("_"):
                accession = f"_{accession}"
            value = template.format(accession=accession)
        else:
            value = template.format(**self.settings)

        return embl_line(line_code, value)

    def project_line(self):
        """
        3.4.3  The PR Line
        The PR (PRoject) line shows the International Nucleotide Sequence
        Database Collaboration (INSDC) Project Identifier that has been
        assigned to the entry.
        Full details of INSDC Project are available at
        http://www.ebi.ac.uk/ena/about/page.php?page=project_guidelines.
        Example:        PR   Project:17285;
        """
        if not self.settings["project_id"]:
            return ""

        line_code = "PR"
        template = "Project:{project_id};"

        return embl_line(line_code, template.format(**self.settings))

    def date_line(self):
        """
        3.4.4  The DT Line
        The DT (DaTe) line shows when an entry first appeared in the database
        and when it was last updated.  Each entry contains two DT lines,
        formatted as follows:
        DT   DD-MON-YYYY (Rel. #, Created)
        DT   DD-MON-YYYY (Rel. #, Last updated, Version #)
        The DT lines from the above example are:
        DT   12-SEP-1991 (Rel. 29, Created)
        DT   13-SEP-1993 (Rel. 37, Last updated, Version 8)
        The date supplied on each DT line indicates when the entry was created
        or Last updated; that will usually also be the date when the new or
        modified Entry became publicly visible via the EBI network servers.
        The release number indicates the first quarterly release made *after*
        the entry was created or last updated. The version number appears only
        on the "Last updated" DT line.
        The absolute value of the version number is of no particular
        significance; its purpose is to allow users to determine easily if the
        version of an entry which they already have is still the most up to
        date version. Version numbers are incremented by one every time an
        entry is updated; since an entry may be updated several times before
        its first appearance in a quarterly release, the version number at the
        time of its first release appearance may be greater than one. Note
        that because an entry may also be updated several times between two
        quarterly releases, there may be gaps in the sequence of version
        numbers which appear in consecutive releases.
        If an entry has not been updated since it was created, it will still
        have two DT lines and the "Last updated" line will have the same date
        (and release number) as the "Created" line.
        """
        line_code = "DT"
        created_template = ("{created:%d-%b-%Y}"
                            " (Rel. {created_release}, Created)")
        updated_template = ("{updated:%d-%b-%Y} (Rel. {updated_release},"
                            " Last updated, Version {version})")

        if not self.settings.get('created', False):
            self.settings['created'] = self.settings['updated']
            self.settings['created_release'] = self.settings['updated_release']

        return embl_line(line_code,
                         created_template.format(**self.settings),
                         add_spacer=False) \
             + embl_line(line_code, updated_template.format(**self.settings))

    def description(self):
        """
        3.4.5  The DE Line
        The DE (Description) lines contain general descriptive information
        about the sequence stored. This may include the designations of genes
        for which thesequence codes, the region of the genome from which it is
        derived, or other information which helps to identify the sequence. The
        format for a DE line is:
        DE   description
        The description is given in ordinary English and is free-format. Often,
        more than one DE line is required; when this is the case, the text is
        divided only between words. The description line from the example above
        is
        DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
        The first DE line generally contains a brief description, which can
        stand alone for cataloguing purposes.
        """
        line_code = "DE"
        template = "{description}"
        if isinstance(self.settings["description"], list):
            self.settings["description"] = \
                " ".join(self.settings["description"])

        return embl_line(line_code, template.format(**self.settings))

    def keywords(self):
        """
        3.4.6  The KW Line
        The KW (KeyWord) lines provide information which can be used to
        generate cross-reference indexes of the sequence entries based on
        functional, structural, or other categories deemed important.
        The format for a KW line is:
         KW   keyword[; keyword ...].
        More than one keyword may be listed on each KW line; the keywords are
        separated by semicolons, and the last keyword is followed by a full
        stop. Keywords may consist of more than one word, and they may contain
        embedded blanks and stops. A keyword is never split between lines.
        An example of a keyword line is:
         KW   beta-glucosidase.
        The keywords are ordered alphabetically; the ordering implies no
        hierarchy of importance or function.  If an entry has no keywords
        assigned to it, it will contain a single KW line like this:
         KW   .
        """
        line_code = "KW"

        return embl_line(line_code, "; ".join(self.settings["keywords"]) + ".",
                         split_on="; ")

    def species(self):
        """
        3.4.7  The OS Line
        The OS (Organism Species) line specifies the preferred scientific name
        of the organism which was the source of the stored sequence. In most
        cases this is done by giving the Latin genus and species designations,
        followed (in parentheses) by the preferred common name in English where
        known. The format is:
         OS   Genus species (name)
        In some cases, particularly for viruses and genetic elements, the only
        accepted designation is a simple name such as
        "Canine adenovirus type 2". In these cases only this designation is
        given. The species line from the example is:
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

        Where only a naturally occurring part of a plasmid is reported, the
        plasmid name will appear on the OG line and the OS/OC lines will
        describe the natural source.
        For example:
         OS   Escherichia coli
         OC   Prokaryota; ... Enterobacteriaceae.
         XX
         OG   Plasmid pUC8
        """
        line_code = "OS"
        template = "{species}"
        if self.settings["species"].isdigit():
            self.settings["species"] = \
                taxid_to_species(self.settings["species"])

        return embl_line(line_code, template.format(**self.settings))

    def classification(self):
        """
        3.4.8  The OC Line
        The OC (Organism Classification) lines contain the taxonomic
        classification of the source organism as described in Section 2.2
        above. The classification is listed top-down as nodes in a taxonomic
        tree in which the most general grouping is given first. The
        classification may be distributed over several OC lines, but nodes are
        not split or hyphenated between lines. The individual items are
        separated by semicolons and the list is terminated by a full stop. The
        format for the OC line is:
         OC   Node[; Node...].

        Example classification lines:
        OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
        OC   euphyllophytes; Spermatophyta; Magnoliophyta; eudicotyledons; Rosidae;
        OC   Fabales; Fabaceae; Papilionoideae; Trifolium.
        """
        line_code = "OC"
        template = "{classification}"

        if self.settings['classification'] is None:
            if self.settings["species"].isdigit():
                taxid = self.settings["species"]
            else:
                taxid = species_to_taxid(self.settings["species"])
            self.settings['classification'] = classification_from_taxid(taxid)

        return embl_line(line_code, template.format(**self.settings))

    def organelle(self):
        """
        3.4.9  The OG Line
        The OG (OrGanelle) linetype indicates the sub-cellular location of
        non-nuclear sequences.  It is only present in entries containing
        non-nuclear sequences and appears after the last OC line in such
        entries.
        The OG line contains
        a) one data item (title cased) from the controlled list detailed under
        the /organelle qualifier definition in the Feature Table Definition
        document that accompanies this release or
        b) a plasmid name.
        Examples include "Mitochondrion", "Plastid:Chloroplast" and
        "Plasmid pBR322".

        For example, a chloroplast sequence from Euglena gracilis would appear
        as:
            OS   Euglena gracilis (green algae)
            OC   Eukaryota; Planta; Phycophyta; Euglenophyceae.
            OG   Plastid:Chloroplast
        """
        if not self.settings["organelle"]:
            return ""

        line_code = "OG"
        template = "{organelle}"
        return embl_line(line_code, template.format(**self.settings))

    def xref(self):
        """
        3.4.11  The DR Line
        The DR (Database Cross-reference) line cross-references other databases
        which contain information related to the entry in which the DR line
        appears. For example, if an annotated/assembled sequence in ENA is cited
        in the IMGT/LIGM database there will be a DR line pointing to the
        relevant IMGT/LIGM entry. The format of the DR line is as follows:
            DR   database_identifier; primary_identifier; secondary_identifier.
        The first item on the DR line, the database identifier, is the
        abbreviated name of the data collection to which reference is made.
        The second item on the DR line, the primary identifier, is a pointer to
        the entry in the external database to which reference is being made.
        The third item on the DR line is the secondary identifier, if available,
        from the referenced database.
        An example of a DR line is shown below:
        DR   MGI; 98599; Tcrb-V4.
        """
        if not self.settings["reference_xref"]:
            return ""

        line_code = "DR"
        template = "{reference_xref};"
        return embl_line(line_code, template.format(**self.settings))

    def comment(self):
        """
        3.4.19  The CC Line
        CC lines are free text comments about the entry, and may be used to
        convey any sort of information thought to be useful that is unsuitable
        for inclusion in other line types.
        """
        if not self.settings["reference_comment"]:
            return ""

        line_code = "CC"
        template = "{reference_comment};"
        return embl_line(line_code, template.format(**self.settings))

    @staticmethod
    def spacer():
        """
        3.4.20 The XX Line
        The XX (spacer) line contains no data or comments. Its purpose is to
        make an entry easier to read on a page or terminal screen by setting
        off the various types of information in appropriate groupings. XX is
        used instead of blank lines to avoid confusion with the sequence data
        lines.
        The XX lines can always be ignored by computer programs.

        ----------
        Currently the spacer is inserted by the 'embl_line' function in
        embl_utilities.py, but this function is kept for completeness.
        """
        return "XX\n"
