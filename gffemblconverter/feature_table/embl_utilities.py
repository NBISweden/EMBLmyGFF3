#!/usr/bin/env python3
"""A small collection of EMBL writing utilities, that is used by the various
parts of the EMBL writer.
"""

import logging
from datetime import datetime

from Bio import Entrez

ENTREZ_EMAIL = "gffconverter@nbis.se"

def embl_line(line_code, information, add_spacer=True):
    """
    This function formats an EMBL output line.

    These lines are formatted as <line code>   <information>, but are limited
    to a maximum of 79 characters.
    ex. PR   Project:17285;
    """
    output = ensure_row_length(f"{line_code:5.4s}{information}\n", 79)

    if add_spacer:
        output += "XX\n"
    return output

def ensure_row_length(line, max_length=79):
    """Function that ensures that no line is longer than the maximum allowed
    line length.

    The function will attempt to split lines between words, but will split
    words if they cannot fit in a single line.
    """
    if len(line) < max_length:
        return line

    output = ""
    line_code = line[:2]
    row = line[:5] # line code and spaces
    words = line[5:].split()
    while words:
        # if a word must be split
        if len(words[0]) > max_length-6:
            split = len(row)
            row += words[0][:max_length-split]
            words[0] = words[0][max_length-split:]

        if len(row + " " + words[0]) > max_length:
            output += row
            output = output.rstrip() + "\n"
            row = f"{line_code}   " + words[0] + " "
        else:
            row += words[0] + " "
        words = words[1:]

    if row:
        output += row
    if output[-1] != "\n":
        output += "\n"

    return output

def get_ena_release(date):
    """
    Returns the correct ENA release version for a given date.

    This far I haven't found a way to get a complete list, or the current
    version from ENA, so for now we rely on the incomplete list.
    """
    release_dates = {136:datetime.strptime("2018-06-11", "%Y-%m-%d"),
                     135:datetime.strptime("2018-03-19", "%Y-%m-%d"),
                     134:datetime.strptime("2018-01-05", "%Y-%m-%d"),
                     133:datetime.strptime("2017-10-04", "%Y-%m-%d"),
                     132:datetime.strptime("2017-05-27", "%Y-%m-%d"),
                     131:datetime.strptime("2017-04-03", "%Y-%m-%d"),
                     130:datetime.strptime("2016-11-13", "%Y-%m-%d"),
                     125:datetime.strptime("2015-09-23", "%Y-%m-%d"),
                     124:datetime.strptime("2015-07-01", "%Y-%m-%d"),
                     123:datetime.strptime("2015-03-23", "%Y-%m-%d"),
                     122:datetime.strptime("2014-12-09", "%Y-%m-%d"),
                     121:datetime.strptime("2014-09-24", "%Y-%m-%d"),
                     120:datetime.strptime("2014-07-01", "%Y-%m-%d"),
                     119:datetime.strptime("2014-03-17", "%Y-%m-%d"),
                     118:datetime.strptime("2013-12-17", "%Y-%m-%d"),
                     117:datetime.strptime("2013-09-12", "%Y-%m-%d"),
                     116:datetime.strptime("2013-06-27", "%Y-%m-%d"),
                     114:datetime.strptime("2012-12-21", "%Y-%m-%d"),
                     }
    for release in sorted(release_dates.keys(), reverse=True):
        if date > release_dates[release]:
            return release
    return 0

def taxid_to_species(taxid):
    """
    Attempts to convert a taxid to a species name, by querying Entrez.
    """
    if taxid.isdigit():
        Entrez.email = ENTREZ_EMAIL
        try:
            search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
            data = Entrez.read(search)
            species = data[0]['ScientificName']
        except IOError as error:
            logging.error(error)
            logging.error("Could not get species from taxid '%s'", taxid)
            species = taxid

    return "%s%s" % (species[0].upper(), species[1:].lower())
