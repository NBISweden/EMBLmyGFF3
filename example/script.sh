#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#FASTA file used to produce the annotation
GENOME="dmel_chr4.fa"

# ANNOATION in gff3 FORMAT
ANNOTATION="dmel_chr4.gff3"

#PROJECT name registered on EMBL
PROJECT="my_EMBL_project"

#Sample accession registered on EMBL
ACCESSION="my_sample_accession"

# species name
SPECIES="Drosophila melanogaster"

# Taxonomy
TAXONOMY="INV"

# The Organism Classification (can be retrieved here http://www.ncbi.nlm.nih.gov/Taxonomy/)
CLASSIFICATION="cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Eremoneura; Cyclorrhapha; Schizophora; Acalyptratae; Ephydroidea; Drosophilidae; Drosophilinae; Drosophilini; Drosophila; Sophophora; melanogaster group; melanogaster subgroup"

# Converter script location
GFF2EMBL="../GFF2EMBL.py"

$GFF2EMBL -a $ACCESSION -p $PROJECT -l $CLASSIFICATION -t linear -s "$SPECIES" -x "$TAXONOMY" $ANNOTATION $GENOME $@ > result.embl
