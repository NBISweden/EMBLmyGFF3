#!/usr/bin/env python3
"""
Reference formatting class for EMBL.

3.4.10  The Reference (RN, RC, RP, RX, RG, RA, RT, RL) Lines
These lines comprise the literature citations within the database.
The citations provide access to the papers from which the data has been
abstracted. The reference lines for a given citation occur in a block, and
are always in the order RN, RC, RP, RX, RG, RA, RT, RL. Within each such
reference block the RN line occurs once, the RC, RP and RX lines occur zero
or more times, and the following lines must occur at least once: the RA (or RG),
RT, RL.  If several references are given, there will be a reference block for
each. Example of references :

RN   [5]
RP   1-1859
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
RT   "Nucleotide and derived amino acid sequence of the cyanogenic
RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.).";
RL   Plant Mol. Biol. 17:209-219(1991).

The formats of the individual lines are explained in the following
paragraphs.
RN   [2]
RP   1-1657990
RG   Prochlorococcus genome consortium
RA   Larimer F.;
RT   ;
RL   Submitted (03-JUL-2003) to the INSDC.
RL   Larimer F., DOE Joint Genome Institute, Production Genomics Facility,
RL   2800 Mitchell Drive, Walnut Creek, CA 94598, USA, and the Genome
RL   Analysis Group, Oak Ridge National Laboratory, 1060 Commerce Park Drive,
RL   Oak Ridge, TN 37831, USA;

"""

from datetime import datetime
from .embl_utilities import embl_line, ensure_quoted

class EMBLReference():
    """
    Reference handler class for EMBL.
    """

    def __init__(self, args=None, number=1):
        self.settings = {"reference_comment":"", "reference_xref":"",
                         "reference_group":"", "reference_author":"",
                         "reference_title":"", "reference_publisher":"",
                         "reference_position":""}
        self.reference_number = number

        if args:
            for key, value in vars(args).items():
                self.settings[key] = value

    def __repr__(self):
        output = ""
        output += self.number()
        output += self.comment()
        output += self.position()
        output += self.x_ref()
        output += self.group()
        output += self.author()
        output += self.title()
        output += self.location()
        output += "XX\n"
        return output

    def number(self):
        """
        3.4.10.1  The RN Line
        The RN (Reference Number) line gives a unique number to each reference
        Citation within an entry. This number is used to designate the
        reference in comments and in the feature table. The format of the RN
        line is:
            RN   [n]
        The reference number is always enclosed in square brackets. Note that
        the set of reference numbers which appear in an entry does not
        necessarily form a continuous sequence from 1 to n, where the entry
        contains "n" references. As references are added to and removed from an
        entry, gaps may be introduced into the sequence of numbers. The
        important point is that once an RN number has been assigned to a
        reference within an entry it never changes. The reference number line
        in the example above is:
            RN   [5]
        """
        line_code = "RN"
        template = f"[{self.reference_number}]"

        return embl_line(line_code, template, False)

    def comment(self):
        """
        3.4.10.2  The RC Line
        The RC (Reference Comment) linetype is an optional linetype which
        appears if The reference has a comment. The comment is in English and
        as many RC lines as are required to display the comment will appear.
        They are formatted thus:
            RC   comment
        """
        if not self.settings['reference_comment']:
            return ""
        line_code = "RC"
        template = "{reference_comment}"

        return embl_line(line_code, template.format(**self.settings), False)

    def position(self):
        """
        3.4.10.3  The RP Line
        The RP (Reference Position) linetype is an optional linetype which
        appears if one or more contiguous base spans of the presented sequence
        can be attributed to the reference in question. As many RP lines as are
        required to display the base span(s) will appear.
        The base span(s) indicate which part(s) of the sequence are covered by
        the reference.  Note that the numbering scheme is for the sequence as
        presented in the database entry (i.e. from 5' to 3' starting at 1), not
        the scheme used by the authors in the reference should the two differ.
        The RP line is formatted thus:
            RP   i-j[, k-l...]
        The RP line in the example above is:
            RP   1-1859
        """
        positions = self.settings['reference_position']
        if not positions:
            if 'sequence_length' in self.settings:
                positions = ["1-{}".format(self.settings['sequence_length'])]
        line_code = "RP"
        template = "{position}"
        value = ",".join([template.format(position=r) for r in positions])

        return embl_line(line_code, value, False)

    def x_ref(self):
        """
        3.4.10.4  The RX Line
        The RX (reference cross-reference) linetype is an optional linetype
        which contains a cross-reference to an external citation or abstract
        resource. For example, if a journal citation exists in the PUBMED
        database, there will be an RX line pointing to the relevant PUBMED
        identifier. The format of the RX line is as follows:
            RX  resource_identifier; identifier.
        The first item on the RX line, the resource identifier, is the
        abbreviated name of the data collection to which reference is made.
        The current set of cross-referenced resources is:
            Resource ID    Fullname
            -----------    ------------------------------------
            PUBMED         PUBMED bibliographic database (NLM)
            DOI            Digital Object Identifier (International DOI
                           Foundation)
            AGRICOLA       US National Agriculture Library (NAL) of the US
                           Department of Agriculture (USDA)
        The second item on the RX line, the identifier, is a pointer to the
        entry in the external resource to which reference is being made. The
        data item used as the primary identifier depends on the resource being
        referenced. For example:
        RX   DOI; 10.1016/0024-3205(83)90010-3.
        RX   PUBMED; 264242.
        Note that further details of DOI are available at http://www.doi.org/.
        URLs formulated in the following way are resolved to the correct full
        text URLs:
            http://dx.doi.org/<doi>
            eg. http:/dx.doi.org/10.1016/0024-3205(83)90010-3
        """
        if not self.settings['reference_xref']:
            return ""
        line_code = "RX"
        template = "{reference_xref}"

        return embl_line(line_code, template.format(**self.settings), False)

    def group(self):
        """
        3.4.10.5  The RG Line
        The RG (Reference Group) lines list the working groups/consortia that
        produced the record. RG line is mainly used in submission reference
        blocks, but could also be used in paper reference if the working group
        is cited as an author in the paper.
        """
        if not self.settings['reference_group']:
            return ""
        line_code = "RG"
        template = "{reference_group}"

        return embl_line(line_code, template.format(**self.settings), False)

    def author(self):
        """
        3.4.10.6  The RA Line
        The RA (Reference Author) lines list the authors of the paper (or other
        work) cited. All of the authors are included, and are listed in the
        order given in the paper. The names are listed surname first followed
        by a blank followed by initial(s) with stops. Occasionally the initials
        may not be known, in which case the surname alone will be listed. The
        author names are separated by commas and terminated by a semicolon;
        they are not split between lines. The RA line in the example is:
        RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
        As many RA lines as necessary are included for each reference.
        """
        if not self.settings['reference_author']:
            return ""
        line_code = "RA"
        template = "{authors}"
        authors = ", ".join(self.settings['reference_author'])

        return embl_line(line_code,
                         template.format(authors=authors),
                         False,
                         ",")

    def title(self):
        """
        3.4.10.7  The RT Line
        The RT (Reference Title) lines give the title of the paper (or other
        work) as exactly as is possible given the limitations of computer
        character sets. Note that the form used is that which would be used in
        a citation rather than that displayed at the top of the published
        paper. For instance, where journals capitalise major title words this
        is not preserved. The title is enclosed in double quotes, and may be
        continued over several lines as necessary. The title lines are
        terminated by a semicolon. The title lines from the example are:
        RT   "Nucleotide and derived amino acid sequence of the cyanogenic
        RT   beta-glucosidase (linamarase) from white clover (Trifolium
        RT   repens L.)";
        Greek letters in titles are spelled out; for example, a title in an
        entry would contain "kappa-immunoglobulin" even though the letter
        itself may be present in the original title. Similar simplifications
        have been made in other cases (e.g. subscripts and superscripts). Note
        that the RT line of a citation which has no title (such as a submission
        to the database) contains only a semicolon.
        """
        line_code = "RT"
        template = "{title};"
        if self.settings['reference_title'] is None:
            self.settings['reference_title'] = ""
        title = ensure_quoted(self.settings['reference_title'])

        return embl_line(line_code, template.format(title=title), False)

    def location(self):
        """
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
        each issue number. In this case the issue number must be included, and
        the format becomes:
            RL   journal vol(no):pp-pp(year).

        If a paper is in press, the RL line will appear with such information as
        we have available, the missing items appearing as zeros. For example:
            RL   Nucleic Acids Res. 0:0-0(2004).
        This indicates a paper which will be published in Nucleic Acids
        Research at some point in 2004, for which we have no volume or page
        information. Such references are updated to include the missing
        information when it becomes available. Another variation of the RL line
        is used for papers found in books or other similar publications, which
        are cited as shown below:
            RA   Birnstiel M., Portmann R., Busslinger M., Schaffner W.,
            RA   Probst E., Kressmeann A.;
            RT   "Functional organization of the histone genes in the
            RT   sea urchin Psammechinus:  A progress report";
            RL   (in) Engberg J., Klenow H., Leick V. (Eds.);
            RL   SPECIFIC EUKARYOTIC GENES:117-132;
            RL   Munksgaard, Copenhagen (1979).
        Note specifically that the line where one would normally encounter the
        journal location is replaced with lines giving the bibliographic
        citation of the book. The first RL line in this case contains the
        designation "(in)", which indicates that this is a book reference.
        The following examples illustrate RL line formats that are used for data
        submissions:
            RL   Submitted (19-NOV-1990) to the INSDC.
            RL   M.A. Hughes, UNIVERSITY OF NEWCASTLE UPON TYNE, MEDICAL SCHOOL,
            RL   NEW CASTLE UPON TYNE, NE2  4HH, UK
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
        The words "Patent number" are followed by the patent application number,
        the patent type (separated by a hyphen), the sequence's serial number
        within the patent (separated by a slash) and the patent application
        date. The subsequent RL lines list the patent applicants, normally
        company names. Finally, for journal publications where no ISSN number is
        available for the journal (proceedings and abstracts, for example), the
        RL line contains the designation "(misc)" as in the following example.
            RL   (misc)Proc. Vth Int. Symp. Biol. Terr. Isopods 2:365-380(2003).
        """
        if not self.settings['reference_publisher']:
            self.settings['reference_publisher'] = \
                f"Submitted ({datetime.today():%d-%b-%Y}) to the INSDC."
        line_code = "RL"
        template = "{reference_publisher}"

        return embl_line(line_code, template.format(**self.settings), False)
