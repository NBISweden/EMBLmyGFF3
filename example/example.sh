#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#FASTA file used to produce the annotation
GENOME="dmel_chr4.fa"

# ANNOATION in gff3 FORMAT
ANNOTATION="dmel_chr4.gff3"

#PROJECT name registered on EMBL
PROJECT="17285"

#Sample accession registered on EMBL
ACCESSION="MY_SAMPLE_ACCESSION"

# species name
SPECIES="Drosophila melanogaster"

# Taxonomy
TAXONOMY="INV"

# The Organism Classification (can be retrieved here http://www.ncbi.nlm.nih.gov/Taxonomy/)
CLASSIFICATION="cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Eremoneura; Cyclorrhapha; Schizophora; Acalyptratae; Ephydroidea; Drosophilidae; Drosophilinae; Drosophilini; Drosophila; Sophophora; melanogaster group; melanogaster subgroup"

#The working groups/consortia that produced the record. No default value
REFERENCE_GROUP="MyGroup"

#Translation table
TABLE="1"

#Data class of the sample.
CLASS="STD"

#Molecule type of the sample.
MOLECULE="genomic DNA"

# Converter script location
GFF3_to_EMBL="../GFF3_to_EMBL.py"

$GFF3_to_EMBL --rg $REFERENCE_GROUP -a $ACCESSION -p $PROJECT -l "$CLASSIFICATION" -m "$MOLECULE" -d $CLASS -r $TABLE -t linear -s "$SPECIES" -x $TAXONOMY $ANNOTATION $GENOME $@ > result.embl
